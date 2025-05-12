#include "tsp.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ilcplex/cplex.h>
#include <string.h>
#include <strings.h>

#define _c(what) if ((error = what)) { \
	fatal_error("CPLEX error in %s: %d\n", #what, error); \
}

void cplex_to_tour(const Instance* inst, const double* xstar, int* tour) {
	int n = inst->num_nodes;
	tour[0] = 0;
	for (int i = 1; i < n; i++) {
		int prev = i > 1 ? tour[i - 2] : -1;
		int curr = tour[i - 1];
		int next = -1;
		for (int j = 0; j < n; j++) {
			if (j == curr) continue; // don't access x[i, i]
			if (j == prev) continue; // don't go backwards
			if (xstar[xpos(inst, curr, j)] > 0.5) {
				if (i > 1) assert(next == -1);
				next = j;
			}
		}
		assert(next != -1);
		tour[i] = next;
	}
}

// estimate the distance between edge numbered i and j in the tour
double edge_dist(const Instance* inst, int i, int j) {
	int a = inst->sol[i], b = inst->sol[(i+1) % inst->num_nodes];
	int c = inst->sol[j], d = inst->sol[(j+1) % inst->num_nodes];
	return (get_cost(inst, a, c) + get_cost(inst, a, d) + get_cost(inst, b, c) + get_cost(inst, b, d)) / 4;
}

// fox edges independently one from the other
void fix_independent(const Instance* inst, double p, bool* fix) {
	int n = inst->num_nodes;
	for (int i = 0; i < n; i++)
		fix[i] = (double) rand() / RAND_MAX < p;
}

// fix a countinuous segment of edges
void fix_segment(const Instance* inst, double p, bool* fix) {
	int n = inst->num_nodes;
	memset(fix, 0, n * sizeof(bool));
	int start = rand() % n, len = round(n * p);
	for (int i = 0; i < len; i++) {
		fix[start] = true;
		start = (start + 1) % n;
	}
}

// leave free a patch of edges close together in space around a center
void fix_hole(const Instance* inst, double p, bool* fix) {
	int n = inst->num_nodes;
	for (int i = 0; i < n; i++) fix[i] = true;
	int center = rand() % n;
	fix[center] = false;
	int nfree = 1, want_free = n - round(n * p);
	while (nfree < want_free) {
		double min_d = INFINITY;
		int min_i = -1;
		for (int i = 0; i < n; i++) {
			if (!fix[i]) continue;
			double d = edge_dist(inst, center, i);
			if (d < min_d) {
				min_d = d;
				min_i = i;
			}
		}
		assert(min_i != -1);
		fix[min_i] = false;
		nfree++;
	}
}

/** sequence_fixings: if false only use independent fixings, otherwise also use segment and hole fixings
  * p0: initial probability of fixing an edge
  * p_decay: decay factor for the probability of fixing an edge
  * iter_tl: time limit for each iteration of the fixing loop, realtive to totale time limit
  */
void hard_fixing_parametrized(Instance *inst, bool seqence_fixings, double p0, double p_decay, double iter_tl) {
	if (inst->time_limit == 0) {
		fatal_error("hard fixing called without a time limit");
	}

	// seed a solution using an heuristic
	double total_tl = inst->time_limit;
	inst->time_limit *= 0.15;
	inst->verbose -= 10;
	tabu_search(inst);
	inst->time_limit = total_tl;
	inst->verbose += 10;

	inst_init_plot(inst);
	inst_plot_cost(inst, inst->sol_cost);

	// setup Cplex
	int n = inst->num_nodes;
	int error;
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	open_cplex(inst, &env, &lp);
	build_base_model(inst, env, lp);
	int num_cols = inst->num_cols = CPXgetnumcols(env, lp);
	install_cplex_callbacks(inst, env, lp);
	mip_warm_start(inst, env, lp); // will use the tabu search result

	// hard fixing loop
	int* indicies = malloc(num_cols * sizeof(int));
	char* bound_t = malloc(num_cols * sizeof(char));
	double* buf_d = malloc(num_cols * sizeof(double));
	int iteration = 0;
	double p = p0;
	while (!is_out_of_time(inst)) {
		iteration++;

		// fix variables
		double* bound_v = buf_d;
		for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) {
			int h = xpos(inst, i, j);
			indicies[h] = h;
			bound_t[h] = 'L';
			bound_v[h] = 0;
		}
		bool fix[n];
		if (iteration % 3 == 0) {
			fix_independent(inst, p, fix);
		} else if (iteration % 3 == 1) {
			fix_segment(inst, p, fix);
		} else {
			fix_hole(inst, p, fix);
		}
		if (inst->verbose >= 90) plot_solution_subset(inst, inst->sol, fix);
		for (int i = 0; i < n; i++) {
			if (fix[i]) {
				int a = inst->sol[i];
				int b = inst->sol[(i+1) % n];
				int h = xpos(inst, a, b);
				bound_v[h] = 1;
			}
		}
		_c(CPXchgbds(env, lp, num_cols, indicies, bound_t, bound_v));

		// solve the problem
		_c(CPXsetdblparam(env, CPX_PARAM_TILIM, min(inst->time_limit * iter_tl, get_remaining_time(inst))));
		double start = get_time();
		_c(CPXmipopt(env, lp));
		double elapsed = get_time() - start;
		int status;
		double* xstar = buf_d;
		double cost;
		_c(CPXsolution(env, lp, &status, &cost, xstar, NULL, NULL, NULL));
		bool timeout = false;
		if (status == CPXMIP_TIME_LIM_FEAS || status == CPXMIP_TIME_LIM_INFEAS) {
			timeout = true;
		} else if (!(status == CPXMIP_OPTIMAL || status == CPXMIP_OPTIMAL_TOL)) {
			fatal_error("MIP optimization failed with status %d\n", status);
		}
		inst_plot_cost(inst, cost);

		// check the new solution
		bool improved = cost < inst->sol_cost;
		int tour[n];
		cplex_to_tour(inst, xstar, tour);
		if (inst->verbose >= 90) 
			plot_solution(inst, tour);
		update_sol(inst, tour, cost);

		debug(30, "iter = %d, \tp = %f, \ttime = %f%c, \tcost = %f%c\n", iteration, p, elapsed, timeout ? '*' : 's', cost, improved ? '*' : ' ');

		// update p
		p = p * p_decay;
		if (improved) p = 1 - (1 - p) * pow(p_decay, 5);
		else p *= pow(p_decay, 1.5);
		if (timeout) p = 1 - (1 - p) * pow(p_decay, 3);
		else p *= pow(p_decay, 0.3);
		p = clamp(p, 0.2, 0.9);
	}

	debug(20, "Hard fixing run for %d iterations\n", iteration);

	free(indicies);
	free(bound_t);
	free(buf_d);
	close_cplex(&env, &lp);
}

void hard_fixing(Instance *inst) {
	hard_fixing_parametrized(inst, true, 0.75, 0.985, 0.05);
}
