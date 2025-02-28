#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "tsp.h"

void parse_arguments(int argc, char *argv[], instance *inst) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {inst->num_nodes = atoi(argv[++i]);}
        if (strcmp(argv[i], "-seed") == 0) {inst->seed = atoi(argv[++i]);}
        if (strcmp(argv[i], "-time") == 0) {inst->time_limit = atoi(argv[++i]);}
        if (strcmp(argv[i], "-file") == 0) {inst->file = argv[++i];}
    }
}



int main(int argc, char *argv[]) {
    instance inst;
    clock_t t1,t2;
    double time;
    t1 = clock();
    if (argc > 1) {
        printf("Arguments received\n");
        parse_arguments(argc, argv, &inst);
        if (parse_tsp_file(inst.file, &inst) == -1) return 1;
    }
    else {
        printf("No arguments received.\n");
        if (parse_tsp_file("data/d198.tsp", &inst) == -1) return 1;
    }
    printf("Data collected\n");
    fill_connections(&inst);
    printf("Connections filled\n");
    plot_instance(&inst);
    printf("Data plotted\n");
    t2 = clock();
    time = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("Time: %fs\n", time);
    return 0;
}
