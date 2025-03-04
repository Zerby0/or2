#ifndef TSP_H
#define TSP_H

//define a epsilon for the double comparison EPS_COST 1e-5
//define infinite for the cost INF_COST 1e38


typedef struct{
    int num_nodes;       // Number of nodes
    double* x_coords;    // Array for x coordinates
    double* y_coords;    // Array for y coordinates
    int seed;            // Seed for random number generation (if used)
    double time_limit;   // Time limit for solving the problem
    int* sol;    // Array for the connections O(n) should not be used for working memory, only store the best solution
    char* file;          // File name

    //using a single data structure to store the coordinates of the pts can
    //be more efficient so the ram can do only one operation (the data separated are very far on the ram)
    //we should try to be very close to maximise the efficency 
} instance;

//use of #define VERBOSE 100 (or use it in the structure to change it at runtime)

int parse_tsp_file(const char* filename, instance *inst); 
void basic_sol(instance* inst);
int plot_instance(instance* inst);


#endif // TSP_H

