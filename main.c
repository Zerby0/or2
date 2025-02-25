#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsp.h"

void parse_arguments(int argc, char *argv[], instance *inst) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {inst->num_nodes = atoi(argv[++i]);}
        if (strcmp(argv[i], "-seed") == 0) {inst->seed = atoi(argv[++i]);}
        if (strcmp(argv[i], "-time") == 0) {inst->time_limit = atoi(argv[++i]);}
    }
}



int main(int argc, char *argv[]) {
    if (argc > 1) {
        printf("Argument received: %s\n", argv[1]);
    } else {
        printf("No arguments received.\n");
    }
    instance inst;
    parse_tsp_file("data/rat99.tsp", &inst);
    printf("Coordinate x of the first one %d\n", inst.x_coords[0]);
    printf("Coordinate y of the first one %d\n", inst.y_coords[0]);
    return 0;
}