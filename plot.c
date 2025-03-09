#include "tsp.h"

#include <stdio.h>

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

int plot_instance(instance* inst) {			// we need to also pass the solution for plotting all the variables solution (example for debugging)
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
    for (int i = 0; i < inst->num_nodes; i++) {
        int idx1 = inst->sol[i];
        int idx2 = inst->sol[(i + 1) % inst->num_nodes]; 
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
