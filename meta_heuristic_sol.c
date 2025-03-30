#include "tsp.h"
#include "hash.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

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

void random_3opt(const Instance* inst, int* tour, double* cost) {
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

void variable_neigh_search(Instance* inst) {
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

// hash functions and comparison functions for cuckoo hash
size_t hash_int(int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    return (x >> 16) ^ x;
}
size_t move_hash1(const void* key, size_t table_size) {
    const Move* m = (const Move*)key;
	size_t i = hash_int(m->i), j = hash_int(m->j);
    return (i * 83492791ULL + j * 3) % table_size;
}
size_t move_hash2(const void* key, size_t table_size) {
    const Move* m = (const Move*)key;
	size_t i = hash_int(m->i), j = hash_int(m->j);
	return (i * 1103515245ULL + j * 7 + 1) % table_size;
}
bool move_compare(const void* pa, const void* pb) {
    const Move* a = (const Move*)pa, *b = (const Move*)pb;
    return (a->i == b->i && a->j == b->j) || (a->i == b->j && a->j == b->i);
}

void best_non_tabu(const Instance* inst, const int* tour, CuckooHash* tabu_list, Move* best_move, double* best_delta) {
	*best_delta = INF_COST;
	for (int i = 0; i < inst->num_nodes - 1; i++) {
		for (int j = i + 1; j < inst->num_nodes; j++) {
			Move move = {i, j};
			if (cuckoo_contains(tabu_list, &move)) {
				continue;
			}
			int a = tour[i];
			int b = tour[i + 1];
			int c = tour[j];
			int d = tour[(j + 1) % inst->num_nodes];
			double cost_next = get_cost(inst, a, c) + get_cost(inst, b, d);
			double cost_curr = get_cost(inst, a, b) + get_cost(inst, c, d);
			double delta = cost_next - cost_curr;
			if (delta < *best_delta) {
				*best_delta = delta;
				best_move->i = i;
				best_move->j = j;
			}
		}
	}
}

void tabu_search(Instance* inst) {
    int n = inst->num_nodes;
	if (inst->time_limit <= 0) {
		debug(5, "WARNING: TABU called without time limit\n");
	}
	if (inst->sol_cost == INF_COST) {
		debug(20, "Initializing solution for TABU with nn\n");
		nearest_neighbor(inst);
	}
	if (inst->plot_cost) list_d_init(&inst->iter_costs);

    // the main parameter of tabue searh is the tenure -> size of the tabu list
    double tenure_min = max(n / 16.0, 2.0);
    double tenure_max = max(n / 4.0, 10.0);
    double tenure_frequency = 200;
	debug(40, "Tabu parameters: %f -> %f, freq = %f\n", tenure_min, tenure_max, tenure_frequency);
	assert(tenure_min <= tenure_max);
	double tenure_amplitude = tenure_max - tenure_min, tenure_center = tenure_min + (tenure_max + tenure_min) / 2;
	double tenure_omega = 2 * M_PI / tenure_frequency;
	debug(80, "Tabu parameters: %f +- %f sin(%f t)\n", tenure_center, tenure_amplitude, tenure_omega);
    
    int current_tour[n];
    memcpy(current_tour, inst->sol, n * sizeof(int));
    double current_cost = inst->sol_cost;
    CuckooHash tabu_list;
    cuckoo_create(&tabu_list, n, 0.2, sizeof(Move), move_hash1, move_hash2, move_compare);
    List_mv move_history;
    list_mv_init(&move_history);
    
    int iteration = 0;
    while (!is_out_of_time(inst)) {
        iteration++;
        double angle = iteration * tenure_omega;
        int tenure = (int)(tenure_center + tenure_amplitude * sin(angle));
        
		// find and apply the best non-tabu move, even if it is not improving
        Move best_move;
        double best_delta;
		best_non_tabu(inst, current_tour, &tabu_list, &best_move, &best_delta);
		two_opt_apply(inst, current_tour, &current_cost, best_move, best_delta);
		update_sol(inst, current_tour, current_cost);
        
		// forbid the move it is not undone in the next iterations
		if (best_delta > -EPS_COST) {
			cuckoo_insert(&tabu_list, &best_move);
			list_mv_push(&move_history, best_move);
			debug(70, "TABU: iteration %d, forbid move: %d <-> %d, delta = %f\n", iteration, best_move.i, best_move.j, best_delta);
		}
        
        while (move_history.len > tenure) {
            Move old_move = move_history.buf[0];
			list_mv_pop_front(&move_history);
            cuckoo_remove(&tabu_list, &old_move);
			debug(70, "TABU: iteration %d, allowing move: %d <-> %d\n", iteration, old_move.i, old_move.j);
        }
		assert(move_history.len == cuckoo_size(&tabu_list));

		if (inst->plot_cost) list_d_push(&inst->iter_costs, current_cost);
    }

	debug(40, "TABU: ran for %d iterations, final cost: %f\n", iteration, current_cost);
    
    cuckoo_free(&tabu_list);
	free(move_history.buf);
}
