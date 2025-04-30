#include "tsp.h"

#include <stdio.h>

static void solve(Instance* inst, bool is_last, char* solver, bool two_opt) {
	inst->sol_cost = INF_COST;
	inst->solver = solver;
	inst->two_opt = two_opt;
	inst->start_time = get_time();
	if (solve_instance(inst) == -1) return;
	double took = get_time() - inst->start_time;
	debug(20, "Solver: %s, Time: %fs, Cost: %f\n", inst->solver, took, inst->sol_cost);
	printf("%f", inst->sol_cost);
	if (!is_last) printf(", ");
}

void run_perf_profile(Instance* inst) {
	debug(10, "Running performance profile\n");
	init_instance_data(inst);
	printf("7, nearest_neighbor, nearest_neighbor_two_opt, nearest_neighbor_all_starts, nearest_neighbor_all_starts_two_opt, extra_milage, extra_milage_two_opt, variable_neigh_search, tabu_search\n"); 
	int base_seed = inst->seed;
	for(int i = base_seed; i < base_seed + 10; i++) {
		debug(15, "Running with seed %d\n", i);
		inst->seed = i;
		random_inst_data(inst);
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
	free_instance_data(inst);
}

void solve_grasp(Instance* inst, int k, int t, bool is_last){
	inst->sol_cost = INF_COST;
	inst->start_time = get_time();
	grasp_parameter(inst, k, t);
	double took = get_time() - inst->start_time;
	debug(20, "Solver: grasp_%d, Time: %fs, Cost: %f\n", k, took, inst->sol_cost);
	printf("%f", inst->sol_cost);
	if (!is_last) printf(", ");
}

void perf_profile_tuning_grasp(Instance* inst) {
	debug(10, "Running performance profile tuning\n");
	init_instance_data(inst);
	printf("12, grasp_2_5, grasp_3_5, grasp_4_5, grasp_5_5, grasp_2_10, grasp_3_10, grasp_4_10, grasp_5_10, grasp_2_20, grasp_3_20, grasp_4_20, grasp_5_20\n"); 
	int base_seed = inst->seed;
	for(int i = base_seed; i < base_seed + 10; i++) {
		debug(15, "Running with seed %d\n", i);
		inst->seed = i;
		random_inst_data(inst);
		printf("Instance_%d, ", i-base_seed);

		solve_grasp(inst, 2, 5, 0);
		solve_grasp(inst, 3, 5,	0);
		solve_grasp(inst, 4, 5,	0);
		solve_grasp(inst, 5, 5,	0);
		solve_grasp(inst, 2, 10, 0);
		solve_grasp(inst, 3, 10, 0);
		solve_grasp(inst, 4, 10, 0);
		solve_grasp(inst, 5, 10, 0);
		solve_grasp(inst, 2, 20, 0);
		solve_grasp(inst, 3, 20, 0);
		solve_grasp(inst, 4, 20, 0);
		solve_grasp(inst, 5, 20, 1);
		

		printf("\n");
	}
	free_instance_data(inst);
}

void run_perf_profile_tuning(Instance* inst){
	//perf_profile_tuning_vns(inst);
	//perf_profile_tuning_tabu(inst);
	perf_profile_tuning_grasp(inst);
}