#include "tsp.h"

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
