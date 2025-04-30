#include "tsp.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define _argi(name, field) \
	else if (strcmp(argv[i], name) == 0) inst->field = atoi(argv[++i])
#define _args(name, field) \
	else if (strcmp(argv[i], name) == 0) inst->field = argv[++i]
#define _argb(name, field) \
	else if (strcmp(argv[i], name) == 0) inst->field = true
#define _argnb(name, field) \
	else if (strcmp(argv[i], name) == 0) inst->field = false
#define _argf(name, field) \
	else if (strcmp(argv[i], name) == 0) inst->field = atof(argv[++i])

int parse_arguments(int argc, char *argv[], Instance *inst) {
	inst->seed = 42;
	inst->file = "data/d198.tsp";
	inst->solver = "basic";
	inst->bc_posting = true;
	inst->bc_fcuts = true;
	
    for (int i = 1; i < argc; i++) {
		if (0);
		_argi("-v", verbose);
		_argi("--verb", verbose);
		_argi("-n", num_nodes);
		_argi("--num-nodes", num_nodes);
		_argi("--seed", seed);
		_args("-f", file);
		_args("--file", file);
		_argb("-r", random_inst);
		_argb("--random", random_inst);
		_args("-s", solver);
		_args("--solver", solver);
		_argb("-2", two_opt);
		_argb("--two-opt", two_opt);
		_argb("--plot-cost", plot_cost);
		_argf("-t", time_limit);
		_argf("--time-limit", time_limit);
		_argb("--perf-profile", perf_profile);
		_argb("--perf-profile-tuning", perf_profile_tuning);
		_args("--write-prob", write_prob);
		_argnb("--no-posting", bc_posting);
		_argnb("--no-fcuts", bc_fcuts);
		else {
			fprintf(stderr, "Unknown option: %s\n", argv[i]);
			return -1;
		}
    }
	return 0;
}

int solve_instance(Instance *inst) {
	if (strcmp(inst->solver, "basic") == 0) basic_sol(inst);
	else if (strcmp(inst->solver, "nn") == 0) nearest_neighbor(inst);
	else if (strcmp(inst->solver, "nna") == 0) nearest_neighbor_all_starts(inst);
	else if (strcmp(inst->solver, "em") == 0) extra_milage(inst);
	else if (strcmp(inst->solver, "vns") == 0) variable_neigh_search(inst);
	else if (strcmp(inst->solver, "tabu") == 0) tabu_search(inst);
	else if (strcmp(inst->solver, "grasp") == 0) grasp(inst);
	else if (strcmp(inst->solver, "benders") == 0) benders_method(inst);
	else if (strcmp(inst->solver, "bc") == 0) branch_and_cut(inst);
	else {
		fprintf(stderr, "Unknown solver: %s\n", inst->solver);
		return -1;
	}
	if (inst->two_opt) two_opt(inst);
	return 0;
}

int main(int argc, char *argv[]) {
    Instance inst_data = {0};
	Instance* inst = &inst_data;

	if (parse_arguments(argc, argv, inst) == -1)
		return -1;

	if (inst->perf_profile) {
		run_perf_profile(inst);
		return 0;
	}

	if(inst->perf_profile_tuning) {
		run_perf_profile_tuning(inst);
		return 0;
	}

	if (inst->random_inst) {
		init_instance_data(inst);
		random_inst_data(inst);
	} else if (inst->file != NULL) {
		if (parse_tsp_file(inst, inst->file) == -1) {
			free_instance_data(inst);
			return -1;
		}
	} else {
		fatal_error("No input specified\n");
	}
    debug(10, "Data collected, Instance size: %d\n", inst->num_nodes);
	srand(inst->seed);
	
	inst->start_time = get_time();
	if (solve_instance(inst) == -1) return -1;
    double took = get_time() - inst->start_time;
    debug(5, "Time: %fs\n", took);
	if (inst->time_patching > 0) debug(10, "Time patching: %fs\n", inst->time_patching);
	if (inst->time_twoopt > 0) debug(10, "Time 2-opt: %fs\n", inst->time_twoopt);
	if (inst->time_fcuts > 0) debug(10, "Time generating fractional cuts: %fs\n", inst->time_fcuts);

	if (inst->sol_cost == INF_COST) {
		fatal_error("No solution found\n");
	}
	printf("%f\n", inst->sol_cost);

	if (inst->plot_cost && inst->iter_data.len > 0) {
		plot_cost_iteration(inst->iter_data.buf, inst->iter_data.len);
		free(inst->iter_data.buf);
	}
    plot_instance(inst);
    debug(20, "Data plotted\n");
	
	free_instance_data(inst);
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
