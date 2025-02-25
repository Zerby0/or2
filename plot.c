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

int plot_instance(instance* inst) {
    int min_x = inst->x_coords[0], max_x = inst->x_coords[0];
	int min_y = inst->y_coords[0], max_y = inst->y_coords[0];
    for (int i = 0; i < inst->num_nodes; i++) {
		min_x = min(min_x, inst->x_coords[i]);
		max_x = max(max_x, inst->x_coords[i]);
		min_y = min(min_y, inst->y_coords[i]);
		max_y = max(max_y, inst->y_coords[i]);
    }

    FILE* pipe = popen("gnuplot -persistent", "w");
	if (!pipe) {
		fprintf(stderr, "Error opening pipe to gnuplot\n");
		return -1;
	}
	fprintf(pipe, "set title 'Nodes plot'\n");
	fprintf(pipe, "set xlabel 'X'\n");
	fprintf(pipe, "set ylabel 'Y'\n");
	fprintf(pipe, "set grid\n");
    const short GAP = 3;
	fprintf(pipe, "set xrange [%d:%d+%d]\n", min_x, max_x, GAP);
	fprintf(pipe, "set yrange [%d:%d+%d]\n", min_y, max_y, GAP);

	// Plot dei punti con le coordinate dalla struttura
	fprintf(pipe, "plot '-' with points pointtype 7 pointsize 1.5 lc rgb 'blue' title 'Nodi'\n");
	for (int i = 0; i < inst->num_nodes; i++) {
		fprintf(pipe, "%d %d\n", inst->x_coords[i], inst->y_coords[i]);
	}
	fprintf(pipe, "e\n");

	fflush(pipe);
	pclose(pipe);
	return 0;
}
