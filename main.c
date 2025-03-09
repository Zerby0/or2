#include "tsp.h"

#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void parse_arguments(int argc, char *argv[], instance *inst) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0) inst->verbose = atoi(argv[++i]);
        if (strcmp(argv[i], "-n") == 0) inst->num_nodes = atoi(argv[++i]);
        if (strcmp(argv[i], "--seed") == 0) inst->seed = atoi(argv[++i]);
        if (strcmp(argv[i], "--time") == 0) inst->time_limit = atoi(argv[++i]);
        if (strcmp(argv[i], "--file") == 0) inst->file = argv[++i];
    }
	if (inst->file == NULL) inst->file = "data/d198.tsp";
}

int main(int argc, char *argv[]) {
    instance inst_data = {0};
	instance* inst = &inst_data;
    clock_t t1,t2; // we must use the world time instead of the cpu time because if the cpu is busy the time will be slower (or parallelize the code)
    double time;
    t1 = clock();
	parse_arguments(argc, argv, inst);
	if (parse_tsp_file(inst, inst->file) == -1) return -1;
    debug(10, "Data collected, instance size: %d\n", inst->num_nodes);
	//nearest_neighbor(inst);
	extra_milage(inst);
    debug(10, "Connections filled\n");
    plot_instance(inst);
    debug(10, "Data plotted\n");
    t2 = clock();
    time = (double)(t2 - t1) / CLOCKS_PER_SEC;
    debug(5, "Time: %fs\n", time);

	printf("%f\n", inst->sol_cost);

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
