#include "tsp.h"

#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define _argi(name, field) \
	else if (strcmp(argv[i], name) == 0) inst->field = atoi(argv[++i])
#define _args(name, field) \
	else if (strcmp(argv[i], name) == 0) inst->field = argv[++i]

void parse_arguments(int argc, char *argv[], instance *inst) {
	inst->seed = time(NULL);
	inst->file = "data/d198.tsp";
	inst->solver = "basic";
    for (int i = 1; i < argc; i++) {
		if (0);
		_argi("-v", verbose);
		_argi("--verb", verbose);
		_argi("-n", num_nodes);
		_argi("--num-nodes", num_nodes);
		_argi("--seed", seed);
		_argi("-t", time_limit);
		_argi("--time", time_limit);
		_args("-f", file);
		_args("--file", file);
		_args("-s", solver);
		_args("--solver", solver);
		else {
			fprintf(stderr, "Unknown option: %s\n", argv[i]);
			exit(1);
		}
    }
}

void solve_instance(instance *inst) {
	if (strcmp(inst->solver, "basic") == 0) basic_sol(inst);
	else if (strcmp(inst->solver, "nn") == 0) nearest_neighbor(inst);
	else if (strcmp(inst->solver, "em") == 0) extra_milage(inst);
	else {
		fprintf(stderr, "Unknown solver: %s\n", inst->solver);
		exit(1);
	}
}

int main(int argc, char *argv[]) {
    instance inst_data = {0};
	instance* inst = &inst_data;
    clock_t t1,t2; // we must use the world time instead of the cpu time because if the cpu is busy the time will be slower (or parallelize the code)
	parse_arguments(argc, argv, inst);
	srand(inst->seed);
	if (parse_tsp_file(inst, inst->file) == -1) return -1;
    debug(10, "Data collected, instance size: %d\n", inst->num_nodes);
    t1 = clock();
	solve_instance(inst);
    t2 = clock();
    debug(10, "Connections filled\n");
    plot_instance(inst);
    debug(10, "Data plotted\n");
    double took = (double)(t2 - t1) / CLOCKS_PER_SEC;
    debug(5, "Time: %fs\n", took);

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
