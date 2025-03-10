#include "tsp.h"

void basic_sol(instance* inst) {
	int tour[inst->num_nodes];
	double cost = 0;
    for (int i = 0; i < inst->num_nodes; i++) {
        tour[i] = i;
		cost += get_cost(inst, i, (i + 1) % inst->num_nodes);
    }
	update_sol(inst, tour, cost);
}
