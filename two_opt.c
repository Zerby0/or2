#include "tsp.h"
#include <math.h>
#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

static void swap(int* a, int pos1, int pos2) {
    int temp = a[pos1];
    a[pos1] = a[pos2];
    a[pos2] = temp;
}

void two_opt(instance* inst) {
    if (inst->num_nodes < 3) {
        printf("Error: instance has less than 3 nodes\n");
        return;
    }
    int* tour = (int*) malloc((inst->num_nodes + 1) * sizeof(int));
    for (int i = 0; i < inst->num_nodes; i++) {
        tour[i] = inst->sol[i];
    }
    tour[inst->num_nodes] = inst->sol[0];
    double cost = inst->sol_cost;
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < inst->num_nodes - 1; i++) {
            for (int j = i + 1; j < inst->num_nodes; j++) {
                int a = tour[i];
                int b = tour[i + 1];
                int c = tour[j];
                int d = tour[(j + 1) % inst->num_nodes];
                double delta = inst->costs_array[a * inst->num_nodes + c] + inst->costs_array[b * inst->num_nodes + d] - inst->costs_array[a * inst->num_nodes + b] - inst->costs_array[c * inst->num_nodes + d];
                if (delta < 0) {
                    improved = true;
                    for (int k = 0; k < (j - i) / 2; k++) {
                        swap(tour, i + 1 + k, j - k);
                    }
                    cost += delta;
                }
            }
        }
    }
    update_sol(inst, tour, cost);
    free(tour);
}