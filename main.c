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
    clock_t t1,t2; // we must use the world time instead of the cpu time because if the cpu is busy the time will be slower (or parallelize the code)
    double time;
    t1 = clock();
    if (argc > 1) {
        printf("Arguments received\n");
        parse_arguments(argc, argv, &inst);
        if (parse_tsp_file(inst.file, &inst) == -1) return 1;
    }
    else {
        printf("No arguments received.\nDefault values will be used.\n");
        if (parse_tsp_file("data/d198.tsp", &inst) == -1) return 1;
    }
    printf("Data collected\n");
    basic_sol(&inst);
    printf("Connections filled\n");
    plot_instance(&inst);
    printf("Data plotted\n");
    t2 = clock();
    time = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("Time: %fs\n", time);
    int rand0_1 = (rand() + 0.0)/RAND_MAX; //random number between 0 and 1
    printf("Random number: %d\n", rand0_1); 
    return 0;
}

//First thing to implement: check feasibility of a solution using a counter array
/*Have to check if a permutation of 1 to n is in the array sol 
    - check if 0 <= sol[h] <= n-1 (very like to not appear)
    - if the first is satisfied we increased a counter ++count[sol[h]]
    - using check.sol(sol,cost,inst)
    - check if the cost is correct (re calculate the cost of the solution) -> if cost int require ==, but if double require a range of values (epsilon defined in 
    the header file)
    
Another function that we hae to implement is the update_best_sol(sol, cost, inst) (adding a parameter best_sol):
    - whenever you have a new sol, check if is feasible, check if the cost is better than the previous one (better), then copy the new solution in best_sol
    
Adding a double best_value which is the cost of the best sol at time t*/