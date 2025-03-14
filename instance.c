#include "tsp.h"
#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>

bool check_sol(const instance* inst, int* tour, double cost) {
	int count[inst->num_nodes];
	memset(count, 0, sizeof(int) * inst->num_nodes);
	for (int i = 0; i < inst->num_nodes; i++) {
		int a = tour[i];
		if (a < 0 || a >= inst->num_nodes) {
			fprintf(stderr, "Invalid tour: node %d is out of range: %d\n", i, a);
			return false;
		}
		count[a]++;
	}
	for (int i = 0; i < inst->num_nodes; i++) {
		if (count[i] != 1) {
			fprintf(stderr, "Invalid tour: node %d appears %d times\n", i, count[i]);
			return false;
		}
	}
	double tot = compute_tour_cost(inst, tour);
	if (fabs(tot - cost) > EPS_COST) {
		fprintf(stderr, "Invalid tour: provided cost %f != recomputed cost %f\n", cost, tot);
		return false;
	}
	return true;
}

void update_sol(instance* inst, int* tour, double cost) {
	assert(tour != NULL);
	assert(check_sol(inst, tour, cost));
	if (inst->sol_cost <= cost) return;
	memcpy(inst->sol, tour, sizeof(int) * inst->num_nodes);
	inst->sol[inst->num_nodes] = inst->sol[0]; // loop back
	inst->sol_cost = cost;
	debug(50, "New sol cost: %f\n", cost);
}

double get_cost(const instance* inst, int i, int j) {
	assert(0 <= i && i < inst->num_nodes);
	assert(0 <= j && j < inst->num_nodes);
	return inst->costs_array[i * inst->num_nodes + j];
}
