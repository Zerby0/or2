#include "tsp.h"
#include <memory.h>
#include <stdio.h>

void update_sol(instance* inst, int* tour, double cost) {
	memcpy(inst->sol, tour, sizeof(int) * inst->num_nodes);
	inst->sol[inst->num_nodes] = inst->sol[0]; // loop back
	inst->sol_cost = cost;
	debug(50, "New sol cost: %f\n", cost);
}
