#ifndef TSP_H
#define TSP_H

typedef struct{
    int num_nodes;       // Number of nodes
    int* x_coords;    // Array for x coordinates
    int* y_coords;    // Array for y coordinates
    int seed;            // Seed for random number generation
    double time_limit;   // Time limit for solving the problem
} instance;


short parse_tsp_file(const char* filename, instance *inst); 
void plot_instance(instance* inst);


#endif // TSP_H
