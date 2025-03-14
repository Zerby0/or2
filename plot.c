#include "tsp.h"

#include <stdio.h>

int plot_partial_sol(const instance* inst, const int* sol, int len) {
    double min_x = inst->x_coords[0], max_x = inst->x_coords[0];
	double min_y = inst->y_coords[0], max_y = inst->y_coords[0];
    for (int i = 0; i < inst->num_nodes; i++) {
		min_x = min(min_x, inst->x_coords[i]);
		max_x = max(max_x, inst->x_coords[i]);
		min_y = min(min_y, inst->y_coords[i]);
		max_y = max(max_y, inst->y_coords[i]);
    }

    FILE* pipe = popen("gnuplot", "w");  //file pipe is a pointer to a file stream, faster than a file
	if (!pipe) {
		fprintf(stderr, "Error opening pipe to gnuplot\n");
		return -1;
	}
	fprintf(pipe, "set title 'Nodes plot'\n");
	fprintf(pipe, "set xlabel 'X'\n");
	fprintf(pipe, "set ylabel 'Y'\n");
	fprintf(pipe, "set grid\n");
    const double GAP = 3;
	fprintf(pipe, "set xrange [%f:%f+%f]\n", min_x, max_x, GAP);
	fprintf(pipe, "set yrange [%f:%f+%f]\n", min_y, max_y, GAP);

	// plot points with coordinates arrays
	fprintf(pipe, "plot '-' with points pointtype 7 pointsize 1.5 lc rgb 'blue' title 'Nodes', '-' with lines lc rgb 'red' title 'Connections'\n");
	for (int i = 0; i < inst->num_nodes; i++) {
		fprintf(pipe, "%f %f\n", inst->x_coords[i], inst->y_coords[i]);
	}
	fprintf(pipe, "e\n");

	// plot the connections
    for (int i = 0; i < len; i++) {
        int idx1 = sol[i];
        int idx2 = sol[(i + 1) % len];
		fprintf(pipe, "%f %f\n", inst->x_coords[idx1], inst->y_coords[idx1]);
        fprintf(pipe, "%f %f\n", inst->x_coords[idx2], inst->y_coords[idx2]);
        fprintf(pipe, "\n"); 
    }
    fprintf(pipe, "e\n");

	// keep gnuplot open until the user closes it
	fprintf(pipe, "pause mouse keypress\n");

	fflush(pipe);
	pclose(pipe);
	return 0;
}

int plot_solution(const instance* inst, const int* sol) {
	return plot_partial_sol(inst, sol, inst->num_nodes);
}

int plot_instance(instance* inst) {
	return plot_partial_sol(inst, inst->sol, inst->num_nodes);
}

//to save the cost of the solution i in a file
void save_cost_to_file(const char *filename, int iteration, double cost) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        perror("Errore apertura file");
        return;
    }
    fprintf(file, "%d %.2f\n", iteration, cost);
    fclose(file);
}

//plot the cost of the solution in a file
void plot_cost_iteration(const char *filename) {
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot == NULL) {
        perror("Errore apertura Gnuplot");
        return;
    }

    fprintf(gnuplot, "set title 'Costo della soluzione - 2-opt'\n");
    fprintf(gnuplot, "set xlabel 'Iterazioni'\n");
    fprintf(gnuplot, "set ylabel 'Costo della soluzione'\n");
    fprintf(gnuplot, "plot '%s' with lines title 'Costo'\n", filename);

    pclose(gnuplot);
}
