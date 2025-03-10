#include "tsp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//cost function (distance from node to node)
double calc_cost(const instance* inst, int i, int j) {
	double dx = inst->x_coords[i] - inst->x_coords[j];
	double dy = inst->y_coords[i] - inst->y_coords[j];
	return sqrt(dx*dx + dy*dy);
}

int parse_tsp_file(instance *inst, const char* filename) {
    if (!filename) {
        fprintf(stderr, "Error: filename is NULL\n");
        return -1;
    }
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return -1;
    }

    char line[256];
    short header = 1;
    char weight_type[256];
    while (fgets(line, sizeof(line), file)) {
        if (header) {
            if (sscanf(line, "DIMENSION : %d", &inst->num_nodes) == 1) {
                inst->x_coords = (double*) malloc(inst->num_nodes * sizeof(double));  //alternative to calloc (everything is set to 0)
                inst->y_coords = (double*) malloc(inst->num_nodes * sizeof(double));
                //I can make an if() to check if the malloc was successful
            }
            else if (sscanf(line, "EDGE_WEIGHT_TYPE : %255s", weight_type) == 1) { 
                if (strcmp(weight_type, "EUC_2D") != 0) {
                    fprintf(stderr, "Error: EDGE_WEIGHT_TYPE is not EUC_2D\n");
                    return -1;
                }
            }
            else if (strcmp(line, "NODE_COORD_SECTION\n") == 0) {
                header = 0;
                if(inst->x_coords == NULL || inst->y_coords == NULL) {
                    fprintf(stderr, "Error: DIMENSION not found\n");
                    return -1;
                }
            }
        }
        else {
            double x, y;
            int i;
            if (sscanf(line, "%d %lf %lf", &i, &x, &y) == 3) {
                inst->x_coords[i-1] = x;
                inst->y_coords[i-1] = y;
            } 
            else if (strcmp(line, "EOF\n") == 0) break;
            else fprintf(stderr, "Error parsing line: %s", line);
        };
    }
    fclose(file);

	// space for the solution
	inst->sol = (int*) malloc((inst->num_nodes + 1) * sizeof(int));
	inst->sol_cost = INF_COST;

    // Precompute the cost array
    inst->costs_array = (double*) malloc(inst->num_nodes * inst->num_nodes * sizeof(double));
    for (int i = 0; i < inst->num_nodes; i++) {
        for (int j = 0; j < inst->num_nodes; j++) {
            inst->costs_array[i * inst->num_nodes + j] = calc_cost(inst, i, j);
        }
    }

    return 0;
}


//I can do a unit test for this function building a main function that calls this function.
