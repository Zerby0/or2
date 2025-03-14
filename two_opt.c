#include "tsp.h"

#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

bool two_opt_once(const instance* inst, int* tour, double* cost) {
	int best_i, best_j;
	double best_delta = 0;
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
	if (best_delta < -EPS_COST) {
		debug(80, "2-opt move: %d <-> %d swapping (%d -- %d, %d -- %d) [%f]\n", best_i, best_j, tour[best_i], tour[best_i + 1], tour[best_j], tour[(best_j + 1) % inst->num_nodes], best_delta);
		invert_subtour(tour, best_i + 1, best_j);
		*cost += best_delta;
		return true;
	} else {
		return false;
	}
}

void two_opt(instance* inst) {
	// TODO time limit
    if (inst->num_nodes < 3) {
        printf("Error: instance has less than 3 nodes\n");
        return;
    }
    int tour[inst->num_nodes];
	memcpy(tour, inst->sol, inst->num_nodes * sizeof(int));
	double cost = inst->sol_cost;
	while (two_opt_once(inst, tour, &cost)) {
		debug(80, "2-opt improved cost = %.2f\n", cost);
	}
    update_sol(inst, tour, cost);
}
