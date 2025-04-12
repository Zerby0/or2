#include "tsp.h"

#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

bool two_opt_best(const Instance* inst, int* tour, bool only_improving, Move* out_move, double* out_delta) {
	int best_i, best_j;
	double best_delta = INF_COST;
	for (int i = 0; i < inst->num_nodes - 1; i++) {
		for (int j = i + 1; j < inst->num_nodes; j++) {
			int a = tour[i];
			int b = tour[i + 1];
			int c = tour[j];
			int d = tour[(j + 1) % inst->num_nodes];
			double cost_next = get_cost(inst, a, c) + get_cost(inst, b, d);
			double cost_curr = get_cost(inst, a, b) + get_cost(inst, c, d);
			double delta = cost_next - cost_curr;
			if (delta < best_delta) {
				best_delta = delta;
				best_i = i;
				best_j = j;
			}
		}
	}
	if (only_improving && best_delta >= -EPS_COST) {
		return false;
	}
	out_move->i = best_i;
	out_move->j = best_j;
	*out_delta = best_delta;
	return true;
}

void two_opt_apply(const Instance* inst, int* tour, double* cost, Move move, double delta) {
	debug(80, "2-opt move: %d <-> %d swapping (%d -- %d, %d -- %d) [%f]\n", move.i, move.j, tour[move.i], tour[move.i + 1], tour[move.j], tour[(move.j + 1) % inst->num_nodes], delta);
	invert_subtour(tour, move.i + 1, move.j);
	*cost += delta;
}

bool two_opt_once(const Instance* inst, int* tour, double* cost) {
	Move move;
	double delta;
	if (two_opt_best(inst, tour, true, &move, &delta)) {
		two_opt_apply(inst, tour, cost, move, delta);
		return true;
	} else {
		return false;
	}
}

void two_opt_from(const Instance* inst, int* tour, double* cost, bool check_time) {
	double cost0 = *cost;
	int iter = 0;
	while (two_opt_once(inst, tour, cost) && !(check_time && is_out_of_time(inst))) {
		iter++;
	}
	debug(40, "2-opt: %d iterations, cost: %f  [%f]\n", iter, *cost, *cost - cost0);
}

void two_opt(Instance* inst) {
    if (inst->num_nodes < 3) {
        printf("Error: instance has less than 3 nodes\n");
        return;
    }
	if (inst->sol_cost == INF_COST) {
		fprintf(stderr, "Error: 2-opt: no solution yet\n");
		exit(1);
	}
    int tour[inst->num_nodes];
	memcpy(tour, inst->sol, inst->num_nodes * sizeof(int));
	double cost = inst->sol_cost;
	two_opt_from(inst, tour, &cost, true);
    update_sol(inst, tour, cost);
}
