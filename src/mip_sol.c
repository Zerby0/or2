#include "tsp.h"

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

#include <ilcplex/cplex.h>
#include <mincut.h>
#include <stdlib.h>

typedef struct {
	int* succ;   // successor in the cycle
	int* comp;   // component number of the element 
	int ncomp;  // number of connected components
} Cycles;

typedef struct {
	int maxrows, maxcols;
	int rcnt, nzcnt;
	double *rhs;
	char *sense;
	int *rmatbeg, *rmatind;
	double *rmatval;
	char **rowname;
	int* purgeable;
	int* local;
} SecData;

#define _c(what) if ((error = what)) { \
	fatal_error("CPLEX error in %s: %d\n", #what, error); \
}

// maps pairs (i, j) with i < j to a unique index 0..N*(N-1)/2 - 1
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
	// for row i, we skip (i+1)*(i+2)/2 elements from previous rows full of unused pairs (k, l) where k < i
	// then within row i, the element is j
	return i * inst->num_nodes - ((i+1) * (i+2))/2 + j;
}

// saves problem description to disk if the user requested it
void write_problem(const Instance *inst, CPXENVptr env, CPXLPptr lp) {
	int error;
	if (inst->write_prob != NULL) {
		debug(50, "Writing problem to %s\n", inst->write_prob);
		error = CPXwriteprob(env, lp, inst->write_prob, "LP");
		if (error) fprintf(stderr, "Warning: CPXwriteprob failed for problem\n");
	}
}

// build the base model for the problem wihtout any SECs
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

// finds all cycles in the integer feasible solution `xstar` and stores them in `cycles`
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

// allocates heap memory for the SEC data structure
void sec_alloc(const Instance* inst, SecData* s) {
	s->maxrows = inst->num_nodes;
	s->maxcols = 2 * inst->num_cols; // sometimes fractional cuts produce more addends then `ncols`
	s->rcnt = s->nzcnt = 0;
	s->rmatbeg = malloc(s->maxrows * sizeof(int));
	s->rmatind = malloc(s->maxcols * sizeof(int));
	s->rmatval = malloc(s->maxcols * sizeof(double));
	s->rhs = malloc(s->maxrows * sizeof(double));
	s->sense = malloc(s->maxrows * sizeof(char));
	s->purgeable = malloc(s->maxrows * sizeof(char));
	s->local = malloc(s->maxrows * sizeof(char));
	if (inst->write_prob) {
		s->rowname = malloc(s->maxrows * sizeof(char*));
		for (int i = 0; i < s->maxrows; i++) {
			s->rowname[i] = malloc(64 * sizeof(char));
		}
	} else {
		s->rowname = NULL;
	}
}

// frees heap memory for the SEC data structure
void sec_free(const Instance* inst, SecData* s) {
	int maxrows = inst->num_nodes;
	free(s->rmatbeg);
	free(s->rmatind);
	free(s->rmatval);
	free(s->rhs);
	free(s->sense);
	free(s->purgeable);
	free(s->local);
	if (s->rowname) {
		for (int i = 0; i < maxrows; i++) {
			free(s->rowname[i]);
		}
		free(s->rowname);
	}
}

// uses cycles information to generate all the SEC constraints, saves them in `s`
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

// uses cycles information to generate & add all SEC constraints to the problem
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
void update_candidate_sol(Instance *inst, const Cycles* cycles, double* objval) {
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

// inverts a cycle in the tour starting at `start`
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

// computes the extra milage for merging cycles at (i,j)
double get_delta1(const Instance *inst, int i, int j, int *succ) {
	return (get_cost(inst, i, j) + get_cost(inst, succ[i], succ[j])) 
		- (get_cost(inst, i, succ[i]) + get_cost(inst, j, succ[j]));
}
double get_delta2(const Instance *inst, int i, int j, int *succ) {
	return (get_cost(inst, i, succ[j]) + get_cost(inst, j, succ[i]))
		- (get_cost(inst, succ[i], i) + get_cost(inst, j, succ[j]));
}

// patches solution in `cycles` in place to only have a single compoent
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
	inst->time_patching += get_time() - start;
	return objval;
}

// uses cycles information to generate all SEC constraints and reject the candidate solution
void callback_sec(const Instance* inst, CPXCALLBACKCONTEXTptr context, const Cycles* cycles) {
	int error;
	SecData s;
	sec_alloc(inst, &s);
	sec_generate(inst, cycles, &s);
	_c(CPXcallbackrejectcandidate(context, s.rcnt, s.nzcnt, s.rhs, s.sense, s.rmatbeg, s.rmatind, s.rmatval));
	sec_free(inst, &s);
}

// uses cycles information to patch the solution and post it to CPLEX
void callback_posting(Instance* inst, CPXCALLBACKCONTEXTptr context, Cycles* cycles) {
	if (!inst->bc_posting) return;
	int n = inst->num_nodes, ncols = inst->num_cols;
	double hc = patch_heuristic(inst, cycles);
	int tour[n];
	cycle_to_tour(inst, cycles, tour);
	if (inst->two_opt) two_opt_from(inst, tour, &hc, false);
	debug(60, "Posting heuristic solution with cost %f\n", hc);
	int* ind = malloc(ncols * sizeof(int));
	double* xheu = malloc(ncols * sizeof(double));
	tour_to_cplex(inst, tour, ind, xheu);
	int error;
	_c(CPXcallbackpostheursoln(context, ncols, ind, xheu, hc, CPXCALLBACKSOLUTION_NOCHECK));
	free(ind);
	free(xheu);
}

// callback at every incumbent candidate solution
void candidate_callback(Instance* inst, CPXCALLBACKCONTEXTptr context) {
	int n = inst->num_nodes, ncols = inst->num_cols;
	int error;
	double* xstar = malloc(ncols * sizeof(double));
	double objval;
	_c(CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval));
	int succ[n], comp[n];
	Cycles cycles = { succ, comp, -1 };
	find_cycles(inst, xstar, &cycles);
	free(xstar);
	debug(90, "Found %d component(s) in callback\n", cycles.ncomp);
	if (cycles.ncomp > 1) {
		callback_posting(inst, context, &cycles);
		callback_sec(inst, context, &cycles);
	}
}

// build the residual graph from the fractional solution `xstar` and store it in `elist` & `ecost`
int build_residual_graph(const Instance* inst, const double* xstar, int* elist, double* ecost) {
	int n = inst->num_nodes;
	int ecount = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (xstar[xpos(inst, i, j)] < 1e-4) continue;
			elist[2*ecount+0] = i;
			elist[2*ecount+1] = j;
			ecost[ecount] = xstar[xpos(inst, i, j)];
			ecount++;
		}
	}
	return ecount;
}

// adds a SEC constraint for the component `comp` to the SEC data structure `s`
void add_sec_for_component(const Instance* inst, SecData* s, int compcount, const int* comp) {
	int n = inst->num_nodes;
	int row = s->rcnt;
	if (row >= s->maxrows) {
		debug(10, "SEC: too many rows\n");
		return;
	}
	if (s->nzcnt + compcount * (compcount - 1) / 2 > s->maxcols) {
		debug(10, "SEC: too many addends (%d >= %d)\n", s->nzcnt + compcount * (compcount - 1) / 2, s->maxcols);
		return;
	}
	s->rcnt++;
	s->rmatbeg[row] = s->nzcnt;
	for (int i = 0; i < compcount; i++) {
		for (int j = i + 1; j < compcount; j++) {
			int a = comp[i], b = comp[j];
			s->rmatind[s->nzcnt] = xpos(inst, a, b);
			s->rmatval[s->nzcnt] = 1.0;
			s->nzcnt++;
		}
	}
	assert(s->nzcnt - s->rmatbeg[row] == compcount * (compcount - 1) / 2);
	s->sense[row] = 'L'; // less than or equal
	s->rhs[row] = compcount - 1.0; // right hand side is |S| - 1
	s->purgeable[row] = CPX_USECUT_FILTER; // can be removed from the cut pool later on
	s->local[row] = 0; // global cut
	if (s->rowname)
		sprintf(s->rowname[0], "SEC_comp%d_size%d", compcount, compcount);
}

// callback when the flow algo finds a violated cut
int flow_callback(double cut_value, int compsize, int* comp, void* handle) {
	void** handle_arr = (void**) handle;
	const Instance* inst = (const Instance*) handle_arr[0];
	SecData* s = (SecData*) handle_arr[1];
	debug(85, "Flow algo found a cut with value %f, size %d\n", cut_value, compsize);
	add_sec_for_component(inst, s, compsize, comp);
	return 0;
}

// CPLEX callback at every relaxation point
void relaxation_callback(Instance* inst, CPXCALLBACKCONTEXTptr context) {
	int error;
	int depth;
	_c(CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEDEPTH, &depth));
	bool is_root = depth == 0;
	if (!is_root) { // alwyas separate at root
		if ((double) rand() / RAND_MAX > inst->bc_theta) // separate at non-root nodes with probability theta
			return;
	}

	double start = get_time();
	int n = inst->num_nodes, ncols = inst->num_cols;
	double* xstar = malloc(ncols * sizeof(double));
	double objval;
	_c(CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols - 1, &objval));
	int* elist = malloc(2 * ncols * sizeof(int));
	double* ecost = malloc(ncols * sizeof(double));
	int ecount = build_residual_graph(inst, xstar, elist, ecost);
	int ncomp;
	int* compscount = NULL, *comps = NULL;
	_c(CCcut_connect_components(n, ecount, elist, ecost, &ncomp, &compscount, &comps));
	debug(80, "Residual graph has %d/%d edges, %d components\n", ecount, ncols, ncomp);
	SecData s;
	sec_alloc(inst, &s);
	if (ncomp > 1) { // more then one compoenent: add a SEC for each component
		int start = 0;
		for (int i = 0; i < ncomp; i++) {
			add_sec_for_component(inst, &s, compscount[i], comps + start);
			start += compscount[i];
		}
	} else { // only one component: run flow algo to find violated cuts
		void* handle[2] = { inst, &s };
		double t = is_root ? 1.9 : 1.7;
		_c(CCcut_violated_cuts(n, ecount, elist, ecost, t, flow_callback, handle));
	}
	_c(CPXcallbackaddusercuts(context, s.rcnt, s.nzcnt, s.rhs, s.sense, s.rmatbeg, s.rmatind, s.rmatval, s.purgeable, s.local));
	sec_free(inst, &s);
	inst->time_fcuts += get_time() - start;
	CC_IFFREE(compscount, int);
	CC_IFFREE(comps, int);
	free(elist);
	free(ecost);
	free(xstar);
}

// installed CPLEX callback, calls the right procedure depending on the context
static int cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
	Instance *inst = (Instance*) userhandle;
	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) candidate_callback(inst, context);
	else if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) relaxation_callback(inst, context);
	else fatal_error("Unknown CPLEX callback context: %lld\n", contextid);
	return 0;
}

// adds the instane incumbent as a warm start to CPLEX; first generates a solution if there is not one yet
void mip_warm_start(Instance* inst, CPXENVptr env, CPXLPptr lp) {
	if (inst->sol_cost == INF_COST) {
		if (!inst->bc_warm) return;
		double start = get_time();
		nearest_neighbor(inst);
		two_opt(inst);
		double elapsed = get_time() - start;
		debug(30, "No solution yet, using NN+2opt for a warm start took %f seconds\n", elapsed);
	}
	debug(30, "Adding a warm start with cost %f\n", inst->sol_cost);
	int ncols = inst->num_cols;
	int* ind = malloc(ncols * sizeof(int));
	double* x = malloc(ncols * sizeof(double));
	tour_to_cplex(inst, inst->sol, ind, x);
	int error;
	int beg = 0;
	int effortlevel = CPX_MIPSTART_NOCHECK;
	_c(CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, x, &effortlevel, NULL));
	free(ind);
	free(x);
}

// opens cplex, creates an empty problem, setups some common parameters
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

	*lp = CPXcreateprob(*env, &error, "TSP");
	if (*lp == NULL) {
		fatal_error("Failed to create CPLEX problem.\n");
	}
}

// installs our callbacks to CPLEX
void install_cplex_callbacks(Instance* inst, CPXENVptr env, CPXLPptr lp) {
	int error;
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (inst->bc_fcuts) contextid |= CPX_CALLBACKCONTEXT_RELAXATION;
	_c(CPXcallbacksetfunc(env, lp, contextid, cplex_callback, inst));

}

// frees CPLEX resources
void close_cplex(CPXENVptr* env, CPXLPptr* lp) {
	if (*lp != NULL) CPXfreeprob(*env, lp);
	if (*env != NULL) CPXcloseCPLEX(env);
}

// solves the instance to optimality using the Benders loop method
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
	double* xstar = malloc(num_cols * sizeof(double));
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
			debug(30, "Patched tour cost: %f\n", cost);
			update_candidate_sol(inst, &cycles, &cost);
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
			update_candidate_sol(inst, &cycles, &lowerbound);
			inst_plot_iter_data(inst, lowerbound, lowerbound); // optimal -> lowerbound == cost
			break;
		} else {
			// multiple cycles -> unfeasible sol -> add one SEC for each cycle
			if (inst->verbose >= 90) {
				plot_infeasible_solution(inst, xstar);
			}
			add_sec_constraints(inst, env, lp, &cycles);
			double cost = patch_heuristic(inst, &cycles); // after SEC since this modifies `cycles`
			debug(30, "Patched tour cost: %f\n", cost);
			update_candidate_sol(inst, &cycles, &cost);
			inst_plot_iter_data(inst, lowerbound, cost);
		}
	}

	debug(20, "Benders required %d iterations\n", iteration);

	close_cplex(&env, &lp);
	free(xstar);
}

// solves the instance to optimality using the branch and cut method
void branch_and_cut(Instance* inst) {
	int n = inst->num_nodes;
	int error;
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	open_cplex(inst, &env, &lp);

	build_base_model(inst, env, lp);
	int num_cols = inst->num_cols = CPXgetnumcols(env, lp);
	install_cplex_callbacks(inst, env, lp);

	mip_warm_start(inst, env, lp);
	double start = get_time();
	_c(CPXmipopt(env, lp));
	double elapsed = get_time() - start;
	debug(30, "CPLEX took %f seconds\n", elapsed);

	int status;
	double cost;
	int iteration = 0;
	double* xstar = malloc(num_cols * sizeof(double));
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
	update_candidate_sol(inst, &cycles, &cost);

	close_cplex(&env, &lp);
	free(xstar);
}
