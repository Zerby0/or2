#include <stdio.h>
#include <stdbool.h>
#include "list.h"

#ifndef TSP_H
#define TSP_H

//define a epsilon for the double comparison EPS_COST 1e-5
//define infinite for the cost INF_COST 1e38
static const double EPS_COST = 1e-5;
static const double INF_COST = 1e38;

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
	// input data
    int num_nodes;     // Number of nodes
    double* x_coords;  // Array for x coordinates
    double* y_coords;  // Array for y coordinates
    double* costs_array; //array for the costs of the edges
	// working memory
    int* sol;          // current best solution; array for the connections O(n) should not be used for working memory, only store the best solution
	double sol_cost;   // cost of the best solution
	double start_time; // time when the program started, to check OOT
	// debugging data
	List_d iter_costs; // list of costs per iteration

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

void swap(int* a, int pos1, int pos2);
void invert_subtour(int* tour, int i, int j);
double compute_tour_cost(const Instance* inst, const int* tour); // O(n)

int parse_tsp_file(Instance *inst, const char* filename);
void random_inst(Instance* inst);
int plot_Instance(Instance* inst);
int plot_solution(const Instance* inst, const int* sol);
int plot_partial_sol(const Instance* inst, const int* sol, int len);
void plot_cost_iteration(double* cost, int len);
void save_cost_to_file(const char* filename, int iteration, double cost);

int solve_Instance(Instance* inst);
bool update_sol(Instance* inst, int* tour, double cost);
double get_cost(const Instance* inst, int i, int j);

void basic_sol(Instance* inst);
void nearest_neighbor(Instance* inst);
void nearest_neighbor_all_starts(Instance* inst);
void extra_milage(Instance* inst);
bool two_opt_best(const Instance* inst, int* tour, bool only_improving, Move* out_move, double* out_delta);
void two_opt_apply(const Instance* inst, int* tour, double* cost, Move move, double delta);
bool two_opt_once(const Instance* inst, int* tour, double* cost);
void two_opt_from(const Instance* inst, int* tour, double* cost);
void two_opt(Instance* inst);
void variable_neigh_search(Instance* inst);
void tabu_search(Instance* inst);

void perf_prof_gen(Instance* inst);
// utility macros

#define debug(level, ...) \
    do { \
        if (inst->verbose >= level) printf(__VA_ARGS__); \
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

#endif // TSP_H
