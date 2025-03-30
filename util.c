#include "tsp.h"

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
