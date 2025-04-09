#include "tsp.h"

#include <assert.h>
#include <stdio.h>

int plot_partial_sol(const Instance* inst, const int* sol, int len) {
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
    double gap = (max_x - min_x + max_y - min_y) * 0.02;
    fprintf(pipe, "set xrange [%f:%f]\n", min_x - gap, max_x + gap);
    fprintf(pipe, "set yrange [%f:%f]\n", min_y - gap, max_y + gap);

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

int plot_solution(const Instance* inst, const int* sol) {
	return plot_partial_sol(inst, sol, inst->num_nodes);
}

int plot_instance(Instance* inst) {
	return plot_partial_sol(inst, inst->sol, inst->num_nodes);
}

void plot_cost_iteration(double* costs, int num_costs) {
	if (num_costs == 0) {
		return;
	}
    FILE *pipe = popen("gnuplot", "w");
    if (pipe == NULL) {
        perror("Error opening Gnuplot");
        return;
    }

    fprintf(pipe, "set title 'Cost of the solution'\n");
    fprintf(pipe, "set xlabel 'Iteration'\n");
    fprintf(pipe, "set ylabel 'Cost'\n");

    fprintf(pipe, "plot '-' with lines linecolor rgb 'blue' title 'Cost', '-' with lines linecolor rgb 'red' title 'Best Cost'\n");
    for (int i = 0; i < num_costs; i++) {
        fprintf(pipe, "%d %f\n", i, costs[i]);
    }
    fprintf(pipe, "e\n");
	double best = costs[0];
    for (int i = 0; i < num_costs; i++) {
		best = min(best, costs[i]);
        fprintf(pipe, "%d %f\n", i, best);
    }
    fprintf(pipe, "e\n");

	fprintf(pipe, "pause mouse keypress\n");
	fflush(pipe);
    pclose(pipe);
}

int plot_infeasible_solution(const Instance* inst, const double* xstar) {
	assert(inst != NULL && xstar != NULL);

    double min_x = inst->x_coords[0], max_x = inst->x_coords[0];
    double min_y = inst->y_coords[0], max_y = inst->y_coords[0];
    for (int i = 1; i < inst->num_nodes; i++) { // Start from 1 as 0 is already used
        min_x = min(min_x, inst->x_coords[i]);
        max_x = max(max_x, inst->x_coords[i]);
        min_y = min(min_y, inst->y_coords[i]);
        max_y = max(max_y, inst->y_coords[i]);
    }

    FILE* pipe = popen("gnuplot", "w");
    if (!pipe) {
        perror("Error opening pipe to gnuplot");
        return -1;
    }

    fprintf(pipe, "set title 'Infeasible Solution Plot (Selected Edges)'\n");
    fprintf(pipe, "set xlabel 'X Coordinate'\n");
    fprintf(pipe, "set ylabel 'Y Coordinate'\n");
    fprintf(pipe, "set grid\n");
    double gap = (max_x - min_x + max_y - min_y) * 0.02;
    fprintf(pipe, "set xrange [%f:%f]\n", min_x - gap, max_x + gap);
    fprintf(pipe, "set yrange [%f:%f]\n", min_y - gap, max_y + gap);

    fprintf(pipe, "plot '-' with points pointtype 7 pointsize 1.2 lc rgb 'blue' title 'Nodes', ");
    fprintf(pipe, "'-' with lines lc rgb 'red' title 'Selected Edges'\n");

	// verticies
    for (int i = 0; i < inst->num_nodes; i++) {
        fprintf(pipe, "%f %f\n", inst->x_coords[i], inst->y_coords[i]);
    }
    fprintf(pipe, "e\n"); // End of data for nodes

	// edges
	int xpos = 0;
    for (int i = 0; i < inst->num_nodes; i++) {
        for (int j = i + 1; j < inst->num_nodes; j++) {
            int k = xpos++;
            if (xstar[k] > 0.5) {
                fprintf(pipe, "%f %f\n", inst->x_coords[i], inst->y_coords[i]);
                fprintf(pipe, "%f %f\n", inst->x_coords[j], inst->y_coords[j]);
                fprintf(pipe, "\n");
            }
        }
    }
    fprintf(pipe, "e\n");

    fprintf(pipe, "pause mouse keypress 'Close window or press any key to exit plot...'\n");
    fflush(pipe);
    if (pclose(pipe) == -1) {
         perror("Error closing gnuplot pipe");
         return -1;
    }

    return 0; 
}
