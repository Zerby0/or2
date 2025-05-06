#include "list.h"
#include "tsp.h"
#include "hash.h"
#include "max_heap.h"


#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double rand_double(void) {
    return (double)rand() / (double)RAND_MAX;
}

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
    
	//tolgono i segmenti (j,h), (j+1,j+2), (k+1,k+2)
    switch (move) {
        case 0:
            invert_subtour(tour, h, j);
            invert_subtour(tour, j + 1, k);
            break;
        case 1:
            invert_subtour(tour, k, j);
            invert_subtour(tour, k-1, h);
            break;
        case 2:
            invert_subtour(tour, h, j);
            invert_subtour(tour, h, k);
            break;
        case 3:
            invert_subtour(tour, j+1, k);
			invert_subtour(tour, h, k);
            break;
    }
    
    *cost = compute_tour_cost(inst, tour);
}

void variable_neigh_search_iteration(Instance* inst, int k, bool incremental) {
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
	inst_init_plot(inst);
	int iter = 1, num_3opts = 0;
	const int MIN_K = 1, MAX_K = 5;
	if (incremental) {
		k = MIN_K;
	}
	double prev_cost = cost;
	while(!is_out_of_time(inst)) {
		// make a kick
		for (int i = 0; i < k; i++)
			random_3opt(inst, tour, &cost);
		num_3opts++, iter++;
		inst_plot_cost(inst, cost);

		// local search with 2opt
		while (two_opt_once(inst, tour, &cost) && !is_out_of_time(inst)) {
			iter++;
			inst_plot_cost(inst, cost);
		}
		update_sol(inst, tour, cost);

		if (incremental) {
			// update k
			if (cost < prev_cost - EPS_COST) {
				k = MIN_K;
			} else if (cost > prev_cost + EPS_COST) {
				k = min(MAX_K, k + 1);
			} else {
				k = max(MIN_K, k - 1);
			}
		} 

		debug(60, "VNS iteration %d, cost = %f -> %f, k = %d\n", iter, prev_cost, cost, k);
		prev_cost = cost;
	}
	debug(40, "VNS: ran for %d iterations, did %d kicks, final cost: %f\n", iter, num_3opts, cost);
}

void variable_neigh_search(Instance* inst) {
	variable_neigh_search_iteration(inst, 1, true);
}

// hash functions and comparison functions for cuckoo hash
size_t hash_int(int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    return (x >> 16) ^ x;
}
size_t move_hash1(const void* key, size_t table_size) {
    const TabuMove* m = (const TabuMove*)key;
	size_t a = hash_int(m->a), b = hash_int(m->b), c = hash_int(m->c), d = hash_int(m->d);
    return (a * 75436877ULL + b * 95094731ull + c * 80067287ull + d) % table_size;
}
size_t move_hash2(const void* key, size_t table_size) {
    const TabuMove* m = (const TabuMove*)key;
	size_t a = hash_int(m->a), b = hash_int(m->b), c = hash_int(m->c), d = hash_int(m->d);
    return (a * 62470469ull + b * 41616881ull + c * 91229777ull + d) % table_size;
}
bool move_compare(const void* pa, const void* pb) {
    const TabuMove* x = (const TabuMove*)pa, *y = (const TabuMove*)pb;
    return x->a == y->a && x->b == y->b && x->c == y->c && x->d == y->d;
}

int compare_int(const void* a, const void* b) {
	return (*(int*)a - *(int*)b);
}
void norm_move(TabuMove* tm) {
	qsort((void*) tm, 4, sizeof(int), compare_int);
	assert(tm->a <= tm->b && tm->b <= tm->c && tm->c <= tm->d);
}

void best_non_tabu(const Instance* inst, const int* tour, CuckooHash* tabu_list, Move* best_move, TabuMove* best_tabu, double* best_delta) {
	*best_delta = INF_COST;
	for (int i = 0; i < inst->num_nodes - 1; i++) {
		for (int j = i + 1; j < inst->num_nodes; j++) {
			int a = tour[i];
			int b = tour[i + 1];
			int c = tour[j];
			int d = tour[(j + 1) % inst->num_nodes];
			double cost_next = get_cost(inst, a, c) + get_cost(inst, b, d);
			double cost_curr = get_cost(inst, a, b) + get_cost(inst, c, d);
			double delta = cost_next - cost_curr;
			if (fabs(delta) < EPS_COST) { // skip moves that do not change the cost
				continue;
			}
			if (delta + EPS_COST >= *best_delta) {
				continue;
			}
			TabuMove tabu = { a, b, c, d };
			norm_move(&tabu);
			if (cuckoo_contains(tabu_list, &tabu)) {
				continue;
			}
			*best_delta = delta;
			*best_move = (Move){ i, j };
			*best_tabu = tabu;
		}
	}
}

void tabu_search_iteration(Instance* inst, double min_factor, double max_factor, double freq) {
    int n = inst->num_nodes;
	if (inst->time_limit <= 0) {
		debug(5, "WARNING: TABU called without time limit\n");
	}
	if (inst->sol_cost == INF_COST) {
		debug(20, "Initializing solution for TABU with nn\n");
		nearest_neighbor(inst);
	}
	inst_init_plot(inst);

    // the main parameter of tabue searh is the tenure -> size of the tabu list
    double tenure_min = max(n * min_factor, 2.0);
    double tenure_max = max(n * max_factor, 10.0);
    double tenure_frequency = freq;
	debug(40, "Tabu parameters: %f -> %f, freq = %f\n", tenure_min, tenure_max, tenure_frequency);
	assert(tenure_min <= tenure_max);
    
    int current_tour[n];
    memcpy(current_tour, inst->sol, n * sizeof(int));
    double current_cost = inst->sol_cost;
    CuckooHash tabu_list;
    cuckoo_create(&tabu_list, n, 0.4, sizeof(TabuMove), move_hash1, move_hash2, move_compare);
    List_tm move_history;
    list_tm_init(&move_history);
    
    int iteration = 0;
    while (!is_out_of_time(inst)) {
        iteration++;
		tenure_max += 0.0005 * tenure_min;
		tenure_frequency += 0.02;
		double tenure_amplitude = (tenure_max - tenure_min) / 2, tenure_center = tenure_min + tenure_amplitude;
		double tenure_omega = 2 * M_PI / tenure_frequency;
        double angle = iteration * tenure_omega;
        int tenure = (int)(tenure_center + tenure_amplitude * sin(angle));
        
		// find and apply the best non-tabu move, even if it is not improving
        Move best_move;
		TabuMove best_tabu;
        double best_delta;
		best_non_tabu(inst, current_tour, &tabu_list, &best_move, &best_tabu, &best_delta);
		two_opt_apply(inst, current_tour, &current_cost, best_move, best_delta);
		bool new_global = update_sol(inst, current_tour, current_cost);
        
		// forbid the move it is not undone in the next iterations
		if (best_delta > -EPS_COST) {
			cuckoo_insert(&tabu_list, &best_tabu);
			list_tm_push(&move_history, best_tabu);
			debug(70, "TABU: iteration %d, forbid move: %d - %d - %d - %d, delta = %f, tenure = %d\n", iteration, best_tabu.a, best_tabu.b, best_tabu.c, best_tabu.d, best_delta, tenure);
		}

		if (new_global) {
			debug(70, "TABU: iteration %d, found new global minima, allowing %d moves\n", iteration, move_history.len);
			for (int i = 0; i < move_history.len; i++)
				cuckoo_remove(&tabu_list, move_history.buf + i);
			move_history.len = 0;
		}
        while (move_history.len > tenure) {
            TabuMove old_move = move_history.buf[0];
			list_tm_pop_front(&move_history);
            cuckoo_remove(&tabu_list, &old_move);
			debug(70, "TABU: iteration %d, allowing move: %d - %d - %d - %d, tenure = %d\n", iteration, old_move.a, old_move.b, old_move.c, old_move.d, tenure);
        }
		assert(move_history.len == cuckoo_size(&tabu_list));

		inst_plot_cost(inst, current_cost);
    }

	debug(20, "TABU: ran for %d iterations, final cost: %f\n", iteration, current_cost);
    
    cuckoo_free(&tabu_list);
	free(move_history.buf);
}

void tabu_search(Instance* inst) {
	tabu_search_iteration(inst, 1./16, 1./4, 200);
}

double grasp_iteration(Instance* inst, int* tour, int start, int max_num_edges, int prob_time) {
	for(int i = 0; i < inst->num_nodes; i++) {
		tour[i] = i;
	}
	swap(tour, 0, start);
	if(max_num_edges < 2) {
		max_num_edges = 2;
	}
	double prob[max_num_edges];
	double tot_cost = 0;
	for (int i = 0; i < inst->num_nodes - 1; i++) {
		Edge heap[max_num_edges];
		int heap_size = 0;
		
		//fill the heap with k minimum edges
		for (int j = i + 1; j < inst->num_nodes; j++) {
			int candidate_node = tour[j];
			double cost = get_cost(inst, tour[i], candidate_node);
			Edge e = { .node = j, .cost = cost };
			insert_edge(heap, &heap_size, max_num_edges, e);
		}

		//calculate the probability of each edge
		double inv_sum = 0;
		for (int k = 0; k < heap_size; k++) {
			inv_sum += 1.0 / heap[k].cost;
		}
		if (inv_sum < 1e-6) inv_sum = 1e-6;
		for (int k = 0; k < heap_size; k++) {
			prob[k] = (1.0 / heap[k].cost) / inv_sum;
		}

		int selected;
		if (rand() % prob_time == 0) {
			double r = rand_double();
			double cumulative = 0;
			for (int k = 0; k < heap_size; k++) {
				cumulative += prob[k];
				if (r <= cumulative) {
					selected = k;
					break;
				}
			}
		} else {
			selected = heap_size - 1;
		}

		debug(80, "[%d] extending %d -> %d [%f]\n", i, tour[i], tour[heap[selected].node], heap[selected].cost);

		swap(tour, i + 1, heap[selected].node);
		tot_cost += heap[selected].cost;
		if (inst->verbose >= 90)
			plot_partial_sol(inst, tour, i + 2);
	}
	int last_node = tour[inst->num_nodes - 1];
	tot_cost += get_cost(inst, last_node, start);
	return tot_cost;
}

void grasp_parameter(Instance* inst, int max_num_edges, int prob_time) {
	int iter_count = 0;
	for (; !is_out_of_time(inst); iter_count++) {
		int tour[inst->num_nodes];
		int start = rand() % inst->num_nodes;
		double cost = grasp_iteration(inst, tour, start, max_num_edges, prob_time);
		if (inst->two_opt) {
			two_opt_from(inst, tour, &cost, true);
		}
		debug(40, "GRASP: iteration %d found cost %f\n", iter_count, cost);
		update_sol(inst, tour, cost);
	}
	debug(30, "GRASP: ran for %d iterations\n", iter_count);
}

void grasp(Instance* inst) {
	grasp_parameter(inst, 2, 10);
}
