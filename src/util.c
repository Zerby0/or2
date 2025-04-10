#include "tsp.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

void swap(int* a, int pos1, int pos2) {
    int temp = a[pos1];
    a[pos1] = a[pos2];
    a[pos2] = temp;
}

void invert_subtour(int* tour, int i, int j) {
	while (i < j) {
		swap(tour, i, j);
		i++;
		j--;
	}
}

double compute_tour_cost(const Instance* inst, const int* tour) {
    double cost = 0.0;
    for (int i = 0; i < inst->num_nodes; i++) {
        int a = tour[i];
        int b = tour[(i + 1) % inst->num_nodes];
        cost += get_cost(inst, a, b);
    }
    return cost;
}

double get_time() {
	return clock() / (double) CLOCKS_PER_SEC;
}

#define _LIST_IMPL(name, type) \
	void list_ ## name ## _init(List_ ## name *l) { \
		l->len = 0; \
		l->capacity = 16; \
		l->buf = (type*) malloc(l->capacity * sizeof(type)); \
	} \
	void list_ ## name ## _push(List_ ## name *l, type val) { \
		if (l->len == l->capacity) { \
			l->capacity *= 2; \
			l->buf = (type*) realloc(l->buf, l->capacity * sizeof(type)); \
		} \
		l->buf[l->len++] = val; \
	} \
	void list_ ## name ## _pop_front(List_ ## name *l) { \
		if (l->len == 0) return; \
		memmove(l->buf, l->buf + 1, (l->len - 1) * sizeof(type)); \
		l->len--; \
	}

_LIST_IMPL(d, double)
_LIST_IMPL(mv, Move)
_LIST_IMPL(tm, TabuMove)
_LIST_IMPL(id, IterData)

void random_points(Instance* inst) {
    srand(inst->seed);
    
    for (int i = 0; i < inst->num_nodes; i++) {
		inst->x_coords[i] = (double)rand() / RAND_MAX * 1000;
		debug(40, "number on %d coord x is %f\n", i, inst->x_coords[i]);
		inst->y_coords[i] = (double)rand() / RAND_MAX * 1000;
		debug(40, "number on %d coord y is %f\n", i, inst->y_coords[i]);
    }
}

// fills instance data with random values, assuming it has already been initialized
void random_inst_data(Instance* inst) {
	random_points(inst);
	
	for (int i = 0; i < inst->num_nodes; i++) {
		for (int j = 0; j < inst->num_nodes; j++) {
			if (i != j) {
				double dx = inst->x_coords[i] - inst->x_coords[j];
				double dy = inst->y_coords[i] - inst->y_coords[j];
				inst->costs_array[i * inst->num_nodes + j] = sqrt(dx * dx + dy * dy);
			} else {
				inst->costs_array[i * inst->num_nodes + j] = INF_COST;
			}
		}
	}
	
	inst->sol_cost = INF_COST;
}
