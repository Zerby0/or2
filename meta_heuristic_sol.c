#include "tsp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void random_3opt(const instance* inst, int* tour, double* cost) {
    if (inst->num_nodes < 6) {
        printf("Error: instance has less than 4 nodes\n");
        return;
    }
    int h, j, k;
    do {
        h = rand() % (inst->num_nodes - 2);
        j = h + 1 + (rand() % (inst->num_nodes - h - 3)); // TODO this can be zero
        k = j + 1 + (rand() % (inst->num_nodes - j - 2));
    } while (k >= inst->num_nodes - 1);
    
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
