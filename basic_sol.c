#include "tsp.h"
#include <stdio.h>

void fill_connections(instance* inst) {
    for (int i = 0; i < inst->num_nodes; i++) {
        inst->connections[i] = i;
    }
    inst->connections[inst->num_nodes] = 0;
}