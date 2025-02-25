#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsp.h"

short parse_tsp_file(const char* filename, instance *inst) {
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
                inst->x_coords = (int*) malloc(inst->num_nodes * sizeof(int));
                inst->y_coords = (int*) malloc(inst->num_nodes * sizeof(int));
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
            int i, x, y;
            if (sscanf(line, "%d %d %d", &i, &x, &y) == 3) {
                inst->x_coords[i-1] = x;
                inst->y_coords[i-1] = y;
            } 
            else if (strcmp(line, "EOF\n") == 0) break;
            else fprintf(stderr, "Error parsing line: %s", line);
        };
    }
    fclose(file);
    return 0;
}

void plot_instance(instance* inst) {

    //idea per efficienza -> sposta il calcolo di max_x e max_y in parse_tsp_file
    //così basta una singola lettura di tsp, anche se logicamente meglio dividere 
    //pars e find max/min in due funzioni diverse
    int max_x = -1; 
    int max_y = -1;
    short GAP = 3;
    for (int i = 0; i < inst->num_nodes; i++) {
        if (inst->x_coords[i] > max_x) max_x = inst->x_coords[i];
        if (inst->y_coords[i] > max_y) max_y = inst->y_coords[i];
    }

    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'Grafico dei Nodi'\n");
        fprintf(gnuplotPipe, "set xlabel 'X'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set xrange [0:%d+%d]\n", max_x, GAP);
        fprintf(gnuplotPipe, "set yrange [0:%d+%d]\n", max_y, GAP);

        // Plot dei punti con le coordinate dalla struttura
        fprintf(gnuplotPipe, "plot '-' with points pointtype 7 pointsize 1.5 lc rgb 'blue' title 'Nodi'\n");
        for (int i = 0; i < inst->num_nodes; i++) {
            fprintf(gnuplotPipe, "%d %d\n", inst->x_coords[i], inst->y_coords[i]);
        }
        fprintf(gnuplotPipe, "e\n");

        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        printf("Errore nell'apertura di Gnuplot.\n");
    }
}
