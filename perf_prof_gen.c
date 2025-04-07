#include "tsp.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

static void solve(Instance* inst, bool is_last, char* solver, bool two_opt) {
	inst->sol_cost = INF_COST;
	inst->solver = solver;
	inst->two_opt = two_opt;
	inst->start_time = get_time();
	if (solve_Instance(inst) == -1) return;
	double took = get_time() - inst->start_time;
	debug(20, "Solver: %s, Time: %fs, Cost: %f\n", inst->solver, took, inst->sol_cost);
	printf("%f", inst->sol_cost);
	if (!is_last) printf(", ");
}

void perf_prof_gen(Instance* inst) {
	debug(10, "Running performance profile\n");
	printf("7, nearest_neighbor, nearest_neighbor_two_opt, nearest_neighbor_all_starts, nearest_neighbor_all_starts_two_opt, extra_milage, extra_milage_two_opt, variable_neigh_search, tabu_search\n"); 
	int base_seed = inst->seed;
	for(int i = base_seed; i < base_seed + 10; i++) {
		debug(15, "Running with seed %d\n", i);
		inst->seed = i;
		random_inst(inst);
		printf("Instance_%d, ", i-base_seed);

		solve(inst, 0, "nn",   0);
        solve(inst, 0, "nn",   1);
		solve(inst, 0, "nna",  0);
        solve(inst, 0, "nna",  1);
		solve(inst, 0, "em",   0);
		solve(inst, 0, "em",   1);
		solve(inst, 0, "vns",  0);
		solve(inst, 1, "tabu", 0);

		printf("\n");
	}
}
