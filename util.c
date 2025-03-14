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
