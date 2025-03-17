#include "tsp.h"
#include <stdlib.h>

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

double compute_tour_cost(const instance* inst, const int* tour) {
    double cost = 0.0;
    for (int i = 0; i < inst->num_nodes; i++) {
        int a = tour[i];
        int b = tour[(i + 1) % inst->num_nodes];
        cost += get_cost(inst, a, b);
    }
    return cost;
}

void list_d_init(list_d *l) {
	l->buf = (double*) malloc(16 * sizeof(double));
	l->len = 0;
	l->capacity = 16;
}

void list_d_push(list_d *l, double val) {
	if (l->len == l->capacity) {
		l->capacity *= 2;
		l->buf = (double*) realloc(l->buf, l->capacity * sizeof(double));
	}
	l->buf[l->len++] = val;
}
