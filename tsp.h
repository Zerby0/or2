#ifndef TSP_H
#define TSP_H

//define a epsilon for the double comparison EPS_COST 1e-5
//define infinite for the cost INF_COST 1e38
#define EPS_COST 1e-5
#define INF_COST 1e38

typedef struct {
	// parameters
	int verbose;       // Verbosity level 1-100
    char* file;        // File name
    int seed;          // Seed for random number generation (if used)
    double time_limit; // Time limit for solving the problem
	char* solver;      // Solver algorithm to use
	// input data
    int num_nodes;     // Number of nodes
    double* x_coords;  // Array for x coordinates
    double* y_coords;  // Array for y coordinates
	double** costs;    // Cost matrix; TODO: maybe use a 1D array?
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
void update_sol(instance* inst, int* tour, double cost);

void basic_sol(instance* inst);
void nearest_neighbor(instance* inst);
void extra_milage(instance* inst);


#endif // TSP_H

