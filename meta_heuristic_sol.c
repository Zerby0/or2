#include "tsp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static void swap(int* a, int pos1, int pos2) {
    int temp = a[pos1];
    a[pos1] = a[pos2];
    a[pos2] = temp;
}

//forse da spostare in un file a parte come 2-opt? filename: opt.c?
void kick_rand_3_opt(instance* inst, int times){
    if (inst->num_nodes < 6) {
        printf("Error: instance has less than 4 nodes\n");
        return;
    }
    int h, j, k;
    for (int i = 0; i < times; i++) {
        do {
            h = rand() % (inst->num_nodes - 3);
            j = h + 1 + (rand() % (inst->num_nodes - h - 3)); 
            k = j + 1 + (rand() % (inst->num_nodes - j - 2));
            debug(80, "3-opt rand value: %d %d %d\n", h, j, k);
        } while (k >= inst->num_nodes - 1);

        int move = rand() % 4;
        switch (move) {
            case 0:
                for (int a = h, b = j; a < b; a++, b--) swap(inst->sol, a, b);
                for (int a = j, b = k; a < b; a++, b--) swap(inst->sol, a, b);
                debug(80, "3-opt kick: %d %d %d\n", h, j, k);
                break;
            case 1:
                for (int a = h, b = j; a < b; a++, b--) swap(inst->sol, a, b);
                for (int a = k, b = inst->num_nodes - 2; a < b; a++, b--) swap(inst->sol, a, b);
                debug(80, "3-opt kick: %d %d %d\n", h, j, k);
                break;
            case 2:
                for (int a = j, b = k; a < b; a++, b--) swap(inst->sol, a, b);
                for (int a = k, b = inst->num_nodes - 2; a < b; a++, b--) swap(inst->sol, a, b);
                debug(80, "3-opt kick: %d %d %d\n", h, j, k);
                break;
            case 3:
                for (int a = h, b = inst->num_nodes - 2; a < b; a++, b--) swap(inst->sol, a, b);
                debug(80, "3-opt kick: %d %d %d\n", h, j, k);
                break;
        }
    }
    //TODO: aggiornare il costo della soluzione in modo sensato
    double tot = 0.0;
    for (int i = 0; i < inst->num_nodes; i++) {
		int a = inst->sol[i], b = inst->sol[(i + 1) % inst->num_nodes];
		tot += get_cost(inst, a, b);
	}
    inst->sol_cost = tot;
}