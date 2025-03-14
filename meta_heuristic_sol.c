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
        printf("Error: instance has less than 4 nodes\n");
        return;
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
            invert_subtour(tour, k, inst->num_nodes - 2);
            break;
        case 2:
            invert_subtour(tour, j, k);
            invert_subtour(tour, k, inst->num_nodes - 2);
            break;
        case 3:
            invert_subtour(tour, h, inst->num_nodes - 2);
            break;
    }
    
    *cost = compute_tour_cost(inst, tour);
}

void variable_neigh_search(instance* inst) {
	if (inst->sol_cost == INF_COST) {
		debug(20, "Initializing solution for VNS with basic_sol\n");
		basic_sol(inst);
	}
	int tour[inst->num_nodes];
	memcpy(tour, inst->sol, inst->num_nodes * sizeof(int));
	double cost = inst->sol_cost;
	for (int iter = 1; iter < 1000; iter++) { // TODO time limit
		for (int i = 0; i < 10; i++)
			random_3opt(inst, tour, &cost);
		while (two_opt_once(inst, tour, &cost)) iter++;
		debug(60, "VNS iteration %d, cost = %.2f\n", iter, cost);
		update_sol(inst, tour, cost);
	}
}
