#include "tsp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void gen_random_triplet(int n, int* h, int* j, int* k) {
    int a, b, c;
    do {
        a = rand() % n;
        b = rand() % n;
        c = rand() % n;
    } while (a == b || b == c || a == c);
	*h = min(a, min(b, c));
	*k = max(a, max(b, c));
	*j = a + b + c - *h - *k;
}

void random_3opt(const instance* inst, int* tour, double* cost) {
    if (inst->num_nodes < 6) {
        printf("Error: instance has less than 6 nodes\n");
		exit(1);
    }
    int h, j, k;
	gen_random_triplet(inst->num_nodes, &h, &j, &k);
    
    int move = rand() % 4;
    debug(80, "3-opt kick type %d: %d %d %d\n", move, h, j, k);
    
    switch (move) {
        case 0:
            invert_subtour(tour, h, j);
            invert_subtour(tour, j, k);
            break;
        case 1:
            invert_subtour(tour, h, j);
            invert_subtour(tour, k, inst->num_nodes - 1);
            break;
        case 2:
            invert_subtour(tour, j, k);
            invert_subtour(tour, k, inst->num_nodes - 1);
            break;
        case 3:
            invert_subtour(tour, h, inst->num_nodes - 1);
            break;
    }
    
    *cost = compute_tour_cost(inst, tour);
}

void variable_neigh_search(instance* inst) {
	if (inst->time_limit <= 0) {
		debug(5, "WARNING: VNS called without time limit\n");
	}
	if (inst->sol_cost == INF_COST) {
		debug(20, "Initializing solution for VNS with nn\n");
		nearest_neighbor(inst);
	}
	int tour[inst->num_nodes];
	memcpy(tour, inst->sol, inst->num_nodes * sizeof(int));
	double cost = inst->sol_cost;
	if (inst->plot_cost) list_d_init(&inst->iter_costs);
	int iter = 1, num_3opts = 0;
	const int MIN_K = 2, MAX_K = 7;
	int k = MIN_K;
	double prev_cost = cost;
	while(!is_out_of_time(inst)) {
		// make a kick
		for (int i = 0; i < k; i++)
			random_3opt(inst, tour, &cost);
		num_3opts++, iter++;
		if (inst->plot_cost) list_d_push(&inst->iter_costs, cost);
		// local search with 2opt
		while (two_opt_once(inst, tour, &cost) && !is_out_of_time(inst)) {
			iter++;
			if (inst->plot_cost) list_d_push(&inst->iter_costs, cost);
		}
		update_sol(inst, tour, cost);
		// update k
		if (cost < prev_cost - EPS_COST) {
			k = MIN_K;
		} else if (cost > prev_cost + EPS_COST) {
			k = max(MIN_K, MIN_K + (k - MIN_K) / 2);
		} else {
			k = min(MAX_K, k + 1);
		}
		debug(60, "VNS iteration %d, cost = %f -> %f, k = %d\n", iter, prev_cost, cost, k);
		prev_cost = cost;
	}
	debug(40, "VNS: ran for %d iterations, did %d kicks, final cost: %f\n", iter, num_3opts, cost);
}
