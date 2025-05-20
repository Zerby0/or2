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
	printf("10, nearest_neighbor, nearest_neighbor_two_opt, nearest_neighbor_all_starts, nearest_neighbor_all_starts_two_opt, extra_milage, extra_milage_two_opt, variable_neigh_search, tabu_search, grasp, grasp_two_opt\n"); 
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
		solve(inst, 0, "tabu", 0);
		solve(inst, 0, "grasp", 0);
		solve(inst, 1, "grasp", 1);

		printf("\n");
	}
	free_instance_data(inst);
}

static void solve_grasp(Instance* inst, int k, int t, bool is_last){
	inst->sol_cost = INF_COST;
	inst->start_time = get_time();
	grasp_parametrized(inst, k, t);
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

static void solve_tabu(Instance* inst, double min_div, double max_div, double freq, bool is_last){
	inst->sol_cost = INF_COST;
	inst->start_time = get_time();
	tabu_search_parametrized(inst, min_div, max_div, freq);
	double took = get_time() - inst->start_time;
	debug(20, "Solver: tabu_%f_%f_%f, Time: %fs, Cost: %f\n", min_div, max_div, freq, took, inst->sol_cost);
	printf("%f", inst->sol_cost);
	if (!is_last) printf(", ");
}

void perf_profile_tuning_tabu(Instance* inst) {
	debug(10, "Running performance profile tuning\n");
	init_instance_data(inst);
	printf("10, ts_0.1_0.4_150, ts_0.2_0.6_200, ts_0.05_0.2_100, ts_0.1_0.4_50, ts_0.15_0.5_120, ts_0.08_0.3_180, ts_0.12_0.35_90, ts_0.2_0.5_160, ts_0.05_0.25_140, ts_0.1_0.45_200\n");
	int base_seed = inst->seed;
	for(int i = base_seed; i < base_seed + 10; i++) {
		debug(15, "Running with seed %d\n", i);
		inst->seed = i;
		random_inst_data(inst);
		printf("Instance_%d, ", i-base_seed);

		solve_tabu(inst, 0.1, 0.4, 150, 0);
		solve_tabu(inst, 0.2, 0.6, 200, 0);
		solve_tabu(inst, 0.05, 0.2, 100, 0);
		solve_tabu(inst, 0.1, 0.4, 50, 0);
		solve_tabu(inst, 0.15, 0.5, 120, 0);
		solve_tabu(inst, 0.08, 0.3, 180, 0);
		solve_tabu(inst, 0.12, 0.35, 90, 0);
		solve_tabu(inst, 0.2, 0.5, 160, 0);
		solve_tabu(inst, 0.05, 0.25, 140, 0);
		solve_tabu(inst, 0.1, 0.45, 200, 1);

		
		printf("\n");
	}
	free_instance_data(inst);
}

static void solve_vns(Instance* inst, int k, bool incremental, bool is_last){
	inst->sol_cost = INF_COST;
	inst->start_time = get_time();
	variable_neigh_search_parametrized(inst, k, incremental);
	double took = get_time() - inst->start_time;
	debug(20, "Solver: vns_%d, Time: %fs, Cost: %f\n", k, took, inst->sol_cost);
	printf("%f", inst->sol_cost);
	if (!is_last) printf(", ");
}
void perf_profile_tuning_vns(Instance* inst) {
	debug(10, "Running performance profile tuning\n");
	init_instance_data(inst);
	printf("10, vns_1, vns_2, vns_3, vns_4, vns_5, vns_6, vns_7, vns_8, vns_9, vns_incr\n");
	int base_seed = inst->seed;
	for(int i = base_seed; i < base_seed + 10; i++) {
		debug(15, "Running with seed %d\n", i);
		inst->seed = i;
		random_inst_data(inst);
		printf("Instance_%d, ", i-base_seed);

		solve_vns(inst, 1, 0, 0);
		solve_vns(inst, 2, 0, 0);
		solve_vns(inst, 3, 0, 0);
		solve_vns(inst, 4, 0, 0);
		solve_vns(inst, 5, 0, 0);
		solve_vns(inst, 6, 0, 0);
		solve_vns(inst, 7, 0, 0);
		solve_vns(inst, 8, 0, 0);
		solve_vns(inst, 9, 0, 0);
		solve_vns(inst, 1, 1, 1);

		
		printf("\n");
	}
	free_instance_data(inst);
}

static void solve_hard_fixing(Instance* inst, bool seqence_fixings, double p0, double p_decay, double iter_nl, bool is_last){
	inst->sol_cost = INF_COST;
	inst->start_time = get_time();
	hard_fixing_parametrized(inst, seqence_fixings, p0, p_decay, 0.1, iter_nl);
	double took = get_time() - inst->start_time;
	debug(20, "Solver: hard_fixing_%f_%f, Time: %fs, Cost: %f\n", p0, p_decay, took, inst->sol_cost);
	printf("%f", inst->sol_cost);
	if (!is_last) printf(", ");
}

void perf_profile_tuning_hard_fixing(Instance* inst) {
	debug(10, "Running performance profile tuning\n");
	init_instance_data(inst);
	printf("18, hard_fixing_0, hard_fixing_1, hard_fixing_2, hard_fixing_3, hard_fixing_4, hard_fixing_5, hard_fixing_6, hard_fixing_7, hard_fixing_8, hard_fixing_9, hard_fixing_10, hard_fixing_12, hard_fixing_13, hard_fixing_14, hard_fixing_15, hard_fixing_16, hard_fixing_17\n");
	int base_seed = inst->seed;
	double p0[3] = {0.5,0.6,0.7};
	double dk[3] = {0.95, 0.97, 0.985};
	double iter_nl[2] = {0.1, 2};
	for(int i = base_seed; i < base_seed + 10; i++) {
		debug(15, "Running with seed %d\n", i);
		inst->seed = i;
		random_inst_data(inst);
		printf("Instance_%d, ", i-base_seed);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 2; l++) {
					bool is_last = (j == 2 && k == 2 && l == 1);
					solve_hard_fixing(inst, true, p0[j], dk[k], iter_nl[l], is_last);
				}
			}
		}
		printf("\n");
	}
	free_instance_data(inst);
}

static void solve_exact(Instance* inst, char* solver, bool posting, bool fcuts, bool warm, bool is_last) {
	inst->sol_cost = INF_COST;
	inst->solver = solver;
	inst->bc_posting = posting;
	inst->bc_fcuts = fcuts;
	inst->bc_warm = warm;
	inst->start_time = get_time();
	if (solve_instance(inst) == -1) return;
	double took = get_time() - inst->start_time;
	debug(20, "Solver: %s, Time: %fs, Cost: %f\n", inst->solver, took, inst->sol_cost);
	printf("%f", took);
	if (!is_last) printf(", ");
}

void run_perf_profile_exact_method(Instance* inst) {
	debug(10, "Running performance profile tuning\n");
	init_instance_data(inst);
	printf("9, benders, bc0_0_0, bc1_0_0, bc0_1_0, bc0_0_1, bc1_1_0, bc1_0_1, bc0_1_1, bc1_1_1\n"); 
	int base_seed = inst->seed;
	for(int i = base_seed; i < base_seed + 10; i++) {
		debug(15, "Running with seed %d\n", i);
		inst->seed = i;
		random_inst_data(inst);
		printf("Instance_%d, ", i-base_seed);

		solve_exact(inst, "benders", 0, 0, 0, 0);
		solve_exact(inst, "bc", 0, 0, 0, 0);
		solve_exact(inst, "bc", 1, 0, 0, 0);
		solve_exact(inst, "bc", 0, 1, 0, 0);
		solve_exact(inst, "bc", 0, 0, 1, 0);
		solve_exact(inst, "bc", 1, 1, 0, 0);
		solve_exact(inst, "bc", 1, 0, 1, 0);
		solve_exact(inst, "bc", 0, 1, 1, 0);
		solve_exact(inst, "bc", 1, 1, 1, 1);
		
		printf("\n");
	}
	free_instance_data(inst);
}

void run_perf_profile_tuning(Instance* inst){
	//perf_profile_tuning_vns(inst);
	//perf_profile_tuning_tabu(inst);
	//perf_profile_tuning_grasp(inst);
	perf_profile_tuning_hard_fixing(inst);
}
