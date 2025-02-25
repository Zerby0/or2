#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsp.h"

void parse_tsp_file(const char* filename, instance *inst) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    char line[256];
    short i = 0;
    int header = 1;
    char* weight_type;
    while (fgets(line, sizeof(line), file)) {
        if (header) {
            if (sscanf(line, "DIMENSION : %d", &inst->num_nodes) == 1) {
                inst->x_coords = (int*) malloc(inst->num_nodes * sizeof(int));
                inst->y_coords = (int*) malloc(inst->num_nodes * sizeof(int));
            }
            else if (sscanf(line, "EDGE_WEIGHT_TYPE : %s", weight_type) == 1) {
                if (strcmp(weight_type, "EUC_2D") != 0) {
                    fprintf(stderr, "Error: EDGE_WEIGHT_TYPE is not EUC_2D\n");
                    exit(1);
                }
            }
            else if (strcmp(line, "NODE_COORD_SECTION\n") == 0) {
                header = 0;
            }
        }
        else {
            int i, x, y;
            if (sscanf(line, "%d %d %d", &i, &x, &y) == 3) {
                inst->x_coords[i-1] = x;
                inst->y_coords[i-1] = y;
            } 
            else if (strcmp(line, "EOF\n") == 0) {
                break;
            }
            else {
                fprintf(stderr, "Error parsing line: %s", line);
            }
        };
    }

    fclose(file);
}
