#include "tsp.h"

#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

void free_instance_data(Instance* inst) {
	debug(60, "Freeing instance data\n");
	if (inst->x_coords) free(inst->x_coords);
	if (inst->y_coords) free(inst->y_coords);
	if (inst->costs_array) free(inst->costs_array);
	if (inst->sol) free(inst->sol);
	inst->x_coords = NULL;
	inst->y_coords = NULL;
	inst->costs_array = NULL;
	inst->sol = NULL;
}

int init_instance_data(Instance* inst) {
	if (inst->num_nodes <= 0) {
		fprintf(stderr, "Error: num_nodes is not set\n");
		return -1;
	}
	debug(60, "Allocating instance data for %d nodes\n", inst->num_nodes);
	int n = inst->num_nodes;
	inst->sol = (int*) malloc((n + 1) * sizeof(int));
    inst->costs_array = (double*) malloc(n * n * sizeof(double));
	inst->x_coords = (double*) malloc(n * sizeof(double));
	inst->y_coords = (double*) malloc(n * sizeof(double));
	if (!inst->x_coords || !inst->y_coords || !inst->costs_array || !inst->sol) {
		fprintf(stderr, "Error: malloc failed\n");
		free_instance_data(inst);
		return -1;
	}
	return 0;
}

bool check_sol(const Instance* inst, int* tour, double cost) {
	int count[inst->num_nodes];
	memset(count, 0, sizeof(int) * inst->num_nodes);
	for (int i = 0; i < inst->num_nodes; i++) {
		int a = tour[i];
		if (a < 0 || a >= inst->num_nodes) {
			fprintf(stderr, "Invalid tour: node %d is out of range: %d\n", i, a);
			return false;
		}
		count[a]++;
	}
	for (int i = 0; i < inst->num_nodes; i++) {
		if (count[i] != 1) {
			fprintf(stderr, "Invalid tour: node %d appears %d times\n", i, count[i]);
			return false;
		}
	}
	double tot = compute_tour_cost(inst, tour);
	if (fabs(tot - cost) > EPS_COST) {
		fprintf(stderr, "Invalid tour: provided cost %f != recomputed cost %f\n", cost, tot);
		return false;
	}
	return true;
}

bool update_sol(Instance* inst, int* tour, double cost) {
	assert(tour != NULL);
	assert(check_sol(inst, tour, cost));
	if (inst->sol_cost <= cost) return false; // no update
	memcpy(inst->sol, tour, sizeof(int) * inst->num_nodes);
	inst->sol[inst->num_nodes] = inst->sol[0]; // loop back
	inst->sol_cost = cost;
	debug(50, "New sol cost: %f\n", cost);
	return true;
}

double get_cost(const Instance* inst, int i, int j) {
	assert(0 <= i && i < inst->num_nodes);
	assert(0 <= j && j < inst->num_nodes);
	return inst->costs_array[i * inst->num_nodes + j];
}

double get_remaining_time(const Instance* inst) {
	if (inst->time_limit <= 0) return INFINITY;
	double elapsed = get_time() - inst->start_time;
	return max(0,  inst->time_limit - elapsed);
}

bool is_out_of_time(const Instance* inst) {
	if (inst->time_limit <= 0) return false;
	bool r = (get_time() - inst->start_time) > inst->time_limit;
	if (r) debug(50, "Out of time!\n");
	return r;
}


void inst_init_plot(Instance* inst) {
	if (!inst->plot_cost) return;
	if (inst->iter_data.len > 0) {
		debug(10, "Warning: erasing %d previous cost data\n", inst->iter_data.len);
	}
	if (inst->iter_data.buf) {
		free(inst->iter_data.buf);
	}
	list_id_init(&inst->iter_data);
}

void inst_plot_iter_data(Instance* inst, double bound, double cost) {
	if (!inst->plot_cost) return;
	IterData data = { bound, cost };
	list_id_push(&inst->iter_data, data);
}

void inst_plot_cost(Instance* inst, double cost) {
	inst_plot_iter_data(inst, NAN, cost);
}
