#include "tsp.h"
#include <assert.h>
#include <memory.h>

void update_sol(instance* inst, int* tour, double cost) {
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
