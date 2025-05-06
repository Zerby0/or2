#include <stdio.h>
#include <stdbool.h>
#include <ilcplex/cplex.h>

#include "list.h"

#ifndef TSP_H
#define TSP_H

//define a epsilon for the double comparison EPS_COST 1e-5
//define infinite for the cost INF_COST 1e38
static const double EPS_COST = 1e-5;
static const double INF_COST = 1e38;

typedef struct {
	double bound; // bound at the current iteration
	double cost;  // cost of the solution at the current iteration
} IterData;
_LIST_DEF(id, IterData);

typedef struct {
	// parameters
	int verbose;       // Verbosity level 1-100
    char* file;        // File name
    int seed;          // Seed for random number generation (if used)
    double time_limit; // Time limit for solving the problem
	char* solver;      // Solver algorithm to use
    bool random_inst;  // Generate random instance
	bool two_opt;      // Apply two-opt algorithm after the solver
	bool plot_cost;    // plot the cost of the solution per iteration
	bool perf_profile; // run the performance profiler
	bool perf_profile_tuning; // run the performance profiler tuning
	char* write_prob;  // write the CPLEX problem to a file
	bool bc_posting;   // do solution posting during B&C using heuristics
	bool bc_fcuts;     // apply fractional cuts during B&C
	bool bc_warm;      // warm start the B&C with an heuristic solution
	double bc_theta;   // B&C: fraction of non-root nodes at which we separate fractional cuts
	// input data
    int num_nodes;     // Number of nodes
	int num_cols;      // Number of columns in the MIP problem
    double* x_coords;  // Array for x coordinates
    double* y_coords;  // Array for y coordinates
    double* costs_array; //array for the costs of the edges
	// working memory
    int* sol;          // current best solution; array for the connections O(n) should not be used for working memory, only store the best solution
	double sol_cost;   // cost of the best solution
	double start_time; // time when the program started, to check OOT
	// debugging data
	List_id iter_data; // list of costs/bounds per iteration
	double time_patching; // time spent in patching
	double time_twoopt;   // time spent in two-opt
	double time_fcuts;   // time spent computing fractional cuts

    //using a single data structure to store the coordinates of the pts can
    //be more efficient so the ram can do only one operation (the data separated are very far on the ram)
    //we should try to be very close to maximise the efficency 
} Instance;

typedef struct {
	// indexes in the tour array, sorted
	int i, j;
} Move;
_LIST_DEF(mv, Move)

typedef struct {
	// verticies ids, sorted
	int a, b, c, d;
} TabuMove;
_LIST_DEF(tm, TabuMove);

double get_time();
bool is_out_of_time(const Instance* inst);
double get_remaining_time(const Instance* inst);
void swap(int* a, int pos1, int pos2);
void invert_subtour(int* tour, int i, int j);
double compute_tour_cost(const Instance* inst, const int* tour); // O(n)

int parse_tsp_file(Instance *inst, const char* filename);
void random_inst_data(Instance* inst);
int plot_instance(Instance* inst);
int plot_solution(const Instance* inst, const int* sol);
int plot_partial_sol(const Instance* inst, const int* sol, int len);
void plot_cost_iteration(const IterData* data, int len);
void save_cost_to_file(const char* filename, int iteration, double cost);
int plot_infeasible_solution(const Instance* inst, const double* xstar);
void plot_solution_subset(const Instance* inst, const int* tour, const bool* subset);

int init_instance_data(Instance* inst);
void free_instance_data(Instance* inst);
int solve_instance(Instance* inst);
bool update_sol(Instance* inst, int* tour, double cost);
double get_cost(const Instance* inst, int i, int j);
void inst_init_plot(Instance* inst);
void inst_plot_cost(Instance* inst, double cost);
void inst_plot_iter_data(Instance* inst, double bound, double cost);

void basic_sol(Instance* inst);
void nearest_neighbor(Instance* inst);
void nearest_neighbor_all_starts(Instance* inst);
void extra_milage(Instance* inst);
bool two_opt_best(const Instance* inst, int* tour, bool only_improving, Move* out_move, double* out_delta);
void two_opt_apply(const Instance* inst, int* tour, double* cost, Move move, double delta);
bool two_opt_once(const Instance* inst, int* tour, double* cost);
void two_opt_from(Instance* inst, int* tour, double* cost, bool check_time);
void two_opt(Instance* inst);
void variable_neigh_search(Instance* inst);
void variable_neigh_search_iteration(Instance* inst, int k, bool incremental);
void tabu_search(Instance* inst);
void tabu_search_iteration(Instance* inst, double min_factor, double max_factor, double freq);
void grasp(Instance* inst);
void grasp_parameter(Instance* inst, int k, int t);
void benders_method(Instance* inst);
void branch_and_cut(Instance* inst);
void hard_fixing(Instance* inst);

void run_perf_profile(Instance* inst);
void run_perf_profile_tuning(Instance* inst);


// cplex stuff

int xpos(const Instance *inst, int i, int j);
void open_cplex(const Instance* inst, CPXENVptr* env, CPXLPptr* lp);
void close_cplex(CPXENVptr* env, CPXLPptr* lp);
void build_base_model(const Instance *inst, CPXENVptr env, CPXLPptr lp);
void install_cplex_callbacks(Instance* inst, CPXENVptr env, CPXLPptr lp);
void mip_warm_start(Instance* inst, CPXENVptr env, CPXLPptr lp);


// utility macros

#define debug(level, ...) \
    do { \
        if (inst->verbose >= level) fprintf(stderr, __VA_ARGS__); \
    } while (0)

#define fatal_error(...) \
	do { \
		fprintf(stderr, __VA_ARGS__); \
		fflush(NULL); \
		exit(1); \
	} while (0)

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

#define clamp(x, l, u) \
({                           \
	__typeof__ (x) _x = (x); \
	__typeof__ (l) _l = (l); \
	__typeof__ (u) _u = (u); \
	_x < _l ? _l : (_x > _u ? _u : _x); \
})

#endif // TSP_H
