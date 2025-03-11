#include <stdio.h>
#include <stdbool.h>

#ifndef TSP_H
#define TSP_H

//define a epsilon for the double comparison EPS_COST 1e-5
//define infinite for the cost INF_COST 1e38
static const double EPS_COST = 1e-5;
static const double INF_COST = 1e38;
static const char* FILE_COST_ITER = "cost_iter.dat";


typedef struct {
	// parameters
	int verbose;       // Verbosity level 1-100
    char* file;        // File name
    int seed;          // Seed for random number generation (if used)
    double time_limit; // Time limit for solving the problem
	char* solver;      // Solver algorithm to use
	bool two_opt;      // Apply two-opt algorithm after the solver
	// input data
    int num_nodes;     // Number of nodes
    double* x_coords;  // Array for x coordinates
    double* y_coords;  // Array for y coordinates
    double* costs_array; //array for the costs of the edges
	// working memory
    int* sol;          // current best solution; array for the connections O(n) should not be used for working memory, only store the best solution
	double sol_cost;   // cost of the best solution

    //using a single data structure to store the coordinates of the pts can
    //be more efficient so the ram can do only one operation (the data separated are very far on the ram)
    //we should try to be very close to maximise the efficency 
} instance;

#define debug(level, ...) \
    do { \
        if (inst->verbose >= level) printf(__VA_ARGS__); \
    } while (0)

int parse_tsp_file(instance *inst, const char* filename); 
int plot_instance(instance* inst);
int plot_solution(const instance* inst, const int* sol);
int plot_partial_sol(const instance* inst, const int* sol, int len);
void plot_cost_iteration(const char* filename);
void save_cost_to_file(const char* filename, int iteration, double cost);

void update_sol(instance* inst, int* tour, double cost);
double get_cost(const instance* inst, int i, int j);

void basic_sol(instance* inst);
void nearest_neighbor(instance* inst);
void extra_milage(instance* inst);
void two_opt(instance* inst);
void kick_rand_3_opt(instance* inst, int times);

#endif // TSP_H
