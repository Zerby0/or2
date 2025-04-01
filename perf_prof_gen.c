#include "tsp.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

static void solve(Instance* inst, char* solver) {
	inst->solver = solver;
	debug(10, "Solver: %s\n", inst->solver);
	inst->start_time = get_time();
	if (solve_Instance(inst) == -1) return;
	double took = get_time() - inst->start_time;
	debug(5, "Time: %fs\n", took);
	printf("%f\n", inst->sol_cost);
    inst->two_opt = false;
}

void perf_prof_gen(Instance* inst) {
	FILE *file = fopen("perf_prof.txt", "w");
    
    if (file == NULL) {
        printf("error opening file\n");
        return;
    }
	fputs("7, nearest_neighbor, nearest_neighbor_two_opt, nearest_neighbor_all_starts, nearest_neighbor_all_starts_two_opt, extra_milage, variable_neigh_search, tabu_search", file); 
	int a = inst->seed;
	for(int i = a; i < a + 10; i++) {
		inst->seed = i;
		random_inst(inst);
		fprintf(file, "\nInstance_%d, ", i-a);

		solve(inst, "nn");
		fprintf(file, "%f, ", inst->sol_cost);
		inst->sol_cost = INF_COST;

        inst->two_opt = true;
        solve(inst, "nn");
		fprintf(file, "%f, ", inst->sol_cost);
		inst->sol_cost = INF_COST;

		solve(inst, "nna");
		fprintf(file, "%f, ", inst->sol_cost);
		inst->sol_cost = INF_COST;
        
        inst->two_opt = true;
        solve(inst, "nna");
		fprintf(file, "%f, ", inst->sol_cost);
		inst->sol_cost = INF_COST;

		solve(inst, "em");
		fprintf(file, "%f, ", inst->sol_cost);
		inst->sol_cost = INF_COST;

		solve(inst, "vns");
		fprintf(file, "%f, ", inst->sol_cost);
		inst->sol_cost = INF_COST;

		solve(inst, "tabu");
		fprintf(file, "%f, ", inst->sol_cost);
		inst->sol_cost = INF_COST;
	}
}