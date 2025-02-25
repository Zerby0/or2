#include "tsp.h"
#include <stdio.h>

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
