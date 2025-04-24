#include "tsp.h"

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

#include <ilcplex/cplex.h>

typedef struct {
	int* succ;   // successor in the cycle
	int* comp;   // component number of the element 
	int ncomp;  // number of connected components
} Cycles;

typedef struct {
	int rcnt, nzcnt;
	double *rhs;
	char *sense;
	int *rmatbeg, *rmatind;
	double *rmatval;
	char **rowname;
} SecData;

#define _c(what) if ((error = what)) { \
	fatal_error("CPLEX error in %s: %d\n", #what, error); \
}

int xpos(const Instance *inst, int i, int j) {
	assert(i >= 0 && i < inst->num_nodes);
	assert(j >= 0 && j < inst->num_nodes);
	assert(i != j);
	// ensure i < j
	if (i > j) {
		int temp = i;
		i = j;
		j = temp;
	}
	// maps pairs (i, j) with i < j to a unique index 0..N*(N-1)/2 - 1
	// for row i, we skip (i+1)*(i+2)/2 elements from previous rows full of unused pairs (k, l) where k < i
	// then within row i, the element is j
	return i * inst->num_nodes - ((i+1) * (i+2))/2 + j;
}

void write_problem(const Instance *inst, CPXENVptr env, CPXLPptr lp) {
	int error;
	if (inst->write_prob != NULL) {
		debug(50, "Writing problem to %s\n", inst->write_prob);
		error = CPXwriteprob(env, lp, inst->write_prob, "LP");
		if (error) fprintf(stderr, "Warning: CPXwriteprob failed for problem\n");
	}
}

void build_base_model(const Instance *inst, CPXENVptr env, CPXLPptr lp) {
	int n = inst->num_nodes;

	int error;
	char binary = 'B';
	char sense_equal;
	double rhs;
	int izero = 0;
	char varname[128];
	char constrname[128];

	// create variables
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			sprintf(varname, "x(%d,%d)", i + 1, j + 1); // start names from 1
			double cost = get_cost(inst, i, j);
			double lb = 0.0, ub = 1.0;
			char* cname[1] = {varname};
			_c(CPXnewcols(env, lp, 1, &cost, &lb, &ub, &binary, cname));
			// check consistency of `xpos` function
			int current_cols = CPXgetnumcols(env, lp);
			if (current_cols - 1 != xpos(inst, i, j)) {
				fatal_error("Variable index mismatch for (%d, %d): expected %d, got %d\n",
					i, j, xpos(inst, i, j), current_cols - 1);
			}
		}
	}

	// create degree constraints
	int index[n];
	double value[n];
	sense_equal = 'E';
	rhs = 2.0;
	for (int i = 0; i < n; i++) {
		sprintf(constrname, "degree(%d)", i + 1); // start names from 1
		int nnz = 0;
		for (int j = 0; j < n; j++) {
			if (i == j) continue;
			index[nnz] = xpos(inst, i, j);
			value[nnz] = 1.0;
			nnz++;
		}
		char* cname[1] = {constrname};
		_c(CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense_equal, &izero, index, value, NULL, cname));
	}
	write_problem(inst, env, lp);
}

void find_cycles(const Instance *inst, const double *xstar, Cycles* cycles) {
	int nnodes = inst->num_nodes;
	int* succ = cycles->succ;
	int* comp = cycles->comp;
	int* ncomp = &cycles->ncomp;

	// reset output structure
	*ncomp = 0;
	for (int i = 0; i < nnodes; i++) {
		succ[i] = -1;
		comp[i] = -1; // -1 indicates node not visited yet [cite: 61]
	}

	// source a DFS at each node
	for (int start_node = 0; start_node < nnodes; start_node++) {
		if (comp[start_node] != -1) continue;

		// start a new components
		(*ncomp)++; // 1-based
		comp[start_node] = *ncomp;
		int current_node = start_node;

		while (1) {
			bool found_next = false;
			for (int next_node = 0; next_node < nnodes; next_node++) {
				if (comp[next_node] == *ncomp) continue; // don't go backwards
				int var_idx = xpos(inst, current_node, next_node);
				if (xstar[var_idx] < 0.5) continue;
					found_next = true;
					succ[current_node] = next_node;
					comp[next_node] = *ncomp;
					current_node = next_node;
					break; // move to next_node & don't look for another successor
				}
			// all neighbours of `current_node` are already in the component
			// -> close the cycle
			if (!found_next) {
				assert(xstar[xpos(inst, current_node, start_node)] > 0.5);
				succ[current_node] = start_node;
				break; // move to the next component
			}
		}
	}
}

void sec_alloc(const Instance* inst, SecData* s) {
	int maxrows = inst->num_nodes;
	int cols = inst->num_cols;
	s->rmatbeg = malloc(maxrows * sizeof(int));
	s->rmatind = malloc(cols * sizeof(int));
	s->rmatval = malloc(cols * sizeof(double));
	s->rhs = malloc(maxrows * sizeof(double));
	s->sense = malloc(maxrows * sizeof(char));
	if (inst->write_prob) {
		s->rowname = malloc(maxrows * sizeof(char*));
		for (int i = 0; i < maxrows; i++) {
			s->rowname[i] = malloc(64 * sizeof(char));
		}
	} else {
		s->rowname = NULL;
	}
}

void sec_free(const Instance* inst, SecData* s) {
	int maxrows = inst->num_nodes;
	free(s->rmatbeg);
	free(s->rmatind);
	free(s->rmatval);
	free(s->rhs);
	free(s->sense);
	if (s->rowname) {
		for (int i = 0; i < maxrows; i++) {
			free(s->rowname[i]);
		}
		free(s->rowname);
	}
}

void sec_generate(const Instance *inst, const Cycles* cycles, SecData* s) {
	int error;
	int n = inst->num_nodes;
	int ncomp = cycles->ncomp;
	assert(ncomp > 1);

	int cols = inst->num_cols;
	s->rcnt = ncomp;
	s->nzcnt = 0;

	for (int c = 1; c <= ncomp; c++) {
		int row = c - 1;
		s->rmatbeg[row] = s->nzcnt;
		int comp_size = 0;
		for (int i = 0; i < n; i++) {
			if (cycles->comp[i] != c) continue;
			comp_size++;
			for (int j = i + 1; j < n; j++) {
				if (cycles->comp[j] != c) continue;
				// i in S, j in S -> e_ij \in E[S]
				s->rmatind[s->nzcnt] = xpos(inst, i, j);
				s->rmatval[s->nzcnt] = 1.0;
				s->nzcnt++;
			}
		}
		assert(comp_size >= 3);
		assert(s->nzcnt - s->rmatbeg[row] == comp_size * (comp_size - 1) / 2);
		s->sense[row] = 'L'; // less than or equal
		s->rhs[row] = comp_size - 1.0; // right hand side is |S| - 1
		if (s->rowname)
			sprintf(s->rowname[row], "SEC_comp%d_size%d", c, comp_size);
	}
}

void add_sec_constraints(const Instance *inst, CPXENVptr env, CPXLPptr lp, const Cycles* cycles) {
	int error;
	SecData s;
	sec_alloc(inst, &s);
	sec_generate(inst, cycles, &s);
	_c(CPXaddrows(env, lp, 0, s.rcnt, s.nzcnt, s.rhs, s.sense, s.rmatbeg, s.rmatind, s.rmatval, NULL, s.rowname));
	sec_free(inst, &s);
	write_problem(inst, env, lp);
}

// transforsm the cycle representation of the solution into a tour representation
void cycle_to_tour(const Instance *inst, const Cycles* cycles, int* tour) {
	assert(cycles->ncomp == 1);
	int n = inst->num_nodes;
	int current = 0;
	for (int i = 0; i < n; ++i) {
		if (i) assert(current != 0);
		tour[i] = current;
		current = cycles->succ[current];
	}
	assert(current == 0);
}

// reconstruct the tour permutation from `cycles`, runs 2-top, updates the incumbent
void build_incumbent_sol(Instance *inst, const Cycles* cycles, double* objval) {
	int n = inst->num_nodes;
	int tour[n];
	cycle_to_tour(inst, cycles, tour);
	*objval = compute_tour_cost(inst, tour); // recompute cost to overcome CPLEX numeric tolerance
	if (inst->two_opt) {
		two_opt_from(inst, tour, objval, false);
	}
	if (inst->verbose >= 90) {
		plot_solution(inst, tour);
	}
	update_sol(inst, tour, *objval);
}

// transforms a feasible solution represented as tour into a CPLEX solution
void tour_to_cplex(const Instance *inst, const int* tour, int* ind, double* x) {
	int n = inst->num_nodes;
	int nc = inst->num_cols;
	for (int j=0; j<nc; j++) ind[j] = j;
	for (int j=0; j<nc; j++) x[j] = 0;
	for (int i=0; i<n; i++) {
		int j = (i + 1) % n;
		int a = tour[i], b = tour[j];
		x[xpos(inst, a, b)] = 1.0;
	}
}

void invert_cycle(const Instance* inst, int *succ, int start) {
	int n = inst->num_nodes;
    int prev = start;
    int current = succ[start];
    int next;

    if (current == start) return;
    int first = start;
    succ[first] = first;
	while (current != start) {
        next = succ[current];
        succ[current] = prev;
        prev = current;
        current = next;
    }
	succ[first] = prev;
}

double get_delta1(const Instance *inst, int i, int j, int *succ) {
	double delta1 = (get_cost(inst, i, j) + get_cost(inst, succ[i], succ[j])) 
		- (get_cost(inst, i, succ[i]) + get_cost(inst, j, succ[j]));
	return delta1;
}

double get_delta2(const Instance *inst, int i, int j, int *succ) {
	double delta2 = (get_cost(inst, i, succ[j]) + get_cost(inst, j, succ[i]))
		- (get_cost(inst, succ[i], i) + get_cost(inst, j, succ[j]));
	return delta2;
}

// modifies `cycles` to only have a single compoent
double patch_heuristic(Instance *inst, Cycles *cycles) {
	if (cycles->ncomp <= 1) return -1;
	double start = get_time();
    int n = inst->num_nodes;
	int *succ = cycles->succ, *comp = cycles->comp;

	// find the best pair of components to merge and merge them
    while (cycles->ncomp > 1) {
        int best_i = -1, best_j = -1;
		bool best_inversion;
        double best_cost = INF_COST;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (comp[i] == comp[j]) continue;
				double delta1 = get_delta1(inst, i, j, succ);
				double delta2 = get_delta2(inst, i, j, succ);
				double delta = min(delta1, delta2);
				if (delta < best_cost) {
					best_cost = delta;
					best_inversion = delta1 < delta2;
					best_i = i;
					best_j = j;
                }
            }
        }
        
		assert(best_i != -1 && best_j != -1);
		debug(80, "Best pair: (%d, %d) on components (%d, %d), cost %f\n", best_i, best_j, comp[best_i], comp[best_j], best_cost);

		int after_i = succ[best_i];
		int after_j = succ[best_j];
		if (best_inversion) {
			debug(80, "Patching: %d -> %d, %d -> %d with invert_cycle\n", best_j, after_i, best_i, after_j);
			invert_cycle(inst, succ, best_j);
			succ[best_i] = best_j;
			succ[after_j] = after_i;
		} else{
			debug(80, "Patching: %d -> %d, %d -> %d\n", best_i, after_j, best_j, after_i);
			succ[best_i] = after_j;
			succ[best_j] = after_i;
		}
		
		int old_comp = comp[best_j];
		int new_comp = comp[best_i];
		for (int h = 0; h < n; ++h) {
			if (comp[h] == old_comp) {
				comp[h] = new_comp;
			}
		}
		cycles->ncomp--;

        // TODO: check if the `cycle` structure is still consistent
	}
	
	// reconstruct the tour
	double objval = 0;
	for (int i = 0; i < n; ++i) {
		objval += get_cost(inst, i, succ[i]);
	}
	debug(30, "Patched tour cost: %f\n", objval);
	inst->time_patching += get_time() - start;
	return objval;
}

static int cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
	Instance *inst = (Instance*) userhandle;
	int n = inst->num_nodes, ncols = inst->num_cols;
	int error;
	double xstar[ncols];
	double objval;
	_c(CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval));
	int succ[n], comp[n];
	Cycles cycles = { succ, comp, -1 };
	find_cycles(inst, xstar, &cycles);
	debug(30, "Found %d component(s) in callback\n", cycles.ncomp);
	if (cycles.ncomp > 1) {
		SecData s;
		sec_alloc(inst, &s);
		sec_generate(inst, &cycles, &s);
		_c(CPXcallbackrejectcandidate(context, s.rcnt, s.nzcnt, s.rhs, s.sense, s.rmatbeg, s.rmatind, s.rmatval));
		sec_free(inst, &s);
		if (inst->bc_posting) {
			double hc = patch_heuristic(inst, &cycles);
			int tour[n], ind[ncols];
			cycle_to_tour(inst, &cycles, tour);
			if (inst->two_opt) two_opt_from(inst, tour, &hc, false);
			debug(60, "Posting heuristic solution with cost %f\n", hc);
			tour_to_cplex(inst, tour, ind, xstar);
			_c(CPXcallbackpostheursoln(context, ncols, ind, xstar, hc, CPXCALLBACKSOLUTION_NOCHECK));
		}
	}
	return 0;
}

void open_cplex(const Instance* inst, CPXENVptr* env, CPXLPptr* lp) {
	int error = 0;
	int n = inst->num_nodes;

	debug(40, "Initializing CPLEX environment...\n");
	*env = CPXopenCPLEX(&error);
	if (*env == NULL) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(*env, error, errmsg);
		fatal_error("Could not open CPLEX environment: %s\n", errmsg);
	}

	CPXsetintparam(*env, CPX_PARAM_RANDOMSEED, inst->seed);
	CPXsetintparam(*env, CPX_PARAM_THREADS, 1); // use a single thread for now
	CPXsetintparam(*env, CPX_PARAM_SCRIND, (inst->verbose >= 50) ? CPX_ON : CPX_OFF);
	CPXsetintparam(*env, CPX_PARAM_MIPDISPLAY, clamp((inst->verbose - 50) / 5, 0, 4));
	if (inst->time_limit > 0)
		CPXsetdblparam(*env, CPX_PARAM_TILIM, get_remaining_time(inst));

	*lp = CPXcreateprob(*env, &error, "TSP_Benders");
	if (*lp == NULL) {
		fatal_error("Failed to create CPLEX problem.\n");
	}
}

void close_cplex(CPXENVptr* env, CPXLPptr* lp) {
	if (*lp != NULL) CPXfreeprob(*env, lp);
	if (*env != NULL) CPXcloseCPLEX(env);
}

void benders_method(Instance *inst) {
	int n = inst->num_nodes;
	int error;
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;

	open_cplex(inst, &env, &lp);

	// start with the base model without any SEC
	build_base_model(inst, env, lp);

	int num_cols = inst->num_cols = CPXgetnumcols(env, lp);
	int iteration = 0;
	double xstar[num_cols];
	inst_init_plot(inst);

	int succ[n], comp[n];
	Cycles cycles = { succ, comp, -1 };

	while (!is_out_of_time(inst)) {
		iteration++;
		debug(35, "Benders Iteration %d...\n", iteration);

		if (inst->time_limit > 0)
			CPXsetdblparam(env, CPX_PARAM_TILIM, get_remaining_time(inst));

		double start = get_time();
		_c(CPXmipopt(env, lp));
		double elapsed = get_time() - start;
		debug(30, "CPLEX took %f seconds\n", elapsed); // TODO: table-like prints

		int status;
		double lowerbound; // when ncomp>1 this is a lowerbound to the TSP cost, otherwise it is the optimal cost up to tolerance
		_c(CPXsolution(env, lp, &status, &lowerbound, xstar, NULL, NULL, NULL));
		_c(CPXgetbestobjval(env, lp, &lowerbound)); // this will be a lowerbound even if mipopt ran out of time
		if (status == CPXMIP_TIME_LIM_FEAS || status == CPXMIP_TIME_LIM_INFEAS) {
			debug(10, "MIP optimization timed out, building a feasible sol with patch heuristic\n");
			find_cycles(inst, xstar, &cycles);
			double cost = patch_heuristic(inst, &cycles);
			build_incumbent_sol(inst, &cycles, &cost);
			inst_plot_iter_data(inst, lowerbound, cost); // this lowerbound could be worse then earlier since cstar is not optimal
			break;
		}
		if (!(status == CPXMIP_OPTIMAL || status == CPXMIP_OPTIMAL_TOL)) {
			fatal_error("MIP optimization failed with status %d\n", status);
		}

		find_cycles(inst, xstar, &cycles);
		debug(30, "Found %d component(s), lowerbound = %f\n", cycles.ncomp, lowerbound);

		if (cycles.ncomp == 1) {
			// only one cycle -> valid tour
			build_incumbent_sol(inst, &cycles, &lowerbound);
			inst_plot_iter_data(inst, lowerbound, lowerbound); // optimal -> lowerbound == cost
			break;
		} else {
			// multiple cycles -> unfeasible sol -> add one SEC for each cycle
			if (inst->verbose >= 90) {
				plot_infeasible_solution(inst, xstar);
			}
			add_sec_constraints(inst, env, lp, &cycles);
			double cost = patch_heuristic(inst, &cycles); // after SEC since this modifies `cycles`
			build_incumbent_sol(inst, &cycles, &cost);
			inst_plot_iter_data(inst, lowerbound, cost);
		}
	}

	debug(20, "Benders required %d iterations\n", iteration);

	close_cplex(&env, &lp);
}

void branch_and_cut(Instance* inst) {
	int n = inst->num_nodes;
	int error;
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	open_cplex(inst, &env, &lp);

	build_base_model(inst, env, lp);
	int num_cols = inst->num_cols = CPXgetnumcols(env, lp);
	
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	_c(CPXcallbacksetfunc(env, lp, contextid, cplex_callback, inst));

	double start = get_time();
	_c(CPXmipopt(env, lp));
	double elapsed = get_time() - start;
	debug(30, "CPLEX took %f seconds\n", elapsed);

	int status;
	double cost;
	int iteration = 0;
	double xstar[num_cols];
	_c(CPXsolution(env, lp, &status, &cost, xstar, NULL, NULL, NULL));
	if (status == CPXMIP_TIME_LIM_FEAS || status == CPXMIP_TIME_LIM_INFEAS) {
		debug(10, "MIP optimization timed out, solution in not proven optimal\n");
	} else if (!(status == CPXMIP_OPTIMAL || status == CPXMIP_OPTIMAL_TOL)) {
		fatal_error("MIP optimization failed with status %d\n", status);
	}

	int succ[n], comp[n];
	Cycles cycles = { succ, comp, -1 };
	find_cycles(inst, xstar, &cycles);
	assert(cycles.ncomp == 1);
	build_incumbent_sol(inst, &cycles, &cost);

	close_cplex(&env, &lp);
}
