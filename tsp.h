#ifndef TSP_H
#define TSP_H

typedef struct{
    int num_nodes;       // Number of nodes
    double* x_coords;    // Array for x coordinates
    double* y_coords;    // Array for y coordinates
    int seed;            // Seed for random number generation
    double time_limit;   // Time limit for solving the problem
    int* connections;    // Array for the connections
    char* file;          // File name

} instance;


int parse_tsp_file(const char* filename, instance *inst); 
void fill_connections(instance* inst);
int plot_instance(instance* inst);

#endif // TSP_H
