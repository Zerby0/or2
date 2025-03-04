#include "tsp.h"
#include <stdio.h>

void basic_sol(instance* inst) {
    for (int i = 0; i < inst->num_nodes; i++) {
        inst->sol[i] = i;
    }
    inst->sol[inst->num_nodes] = 0;
}