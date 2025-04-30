#include <stdio.h>
#include "max_heap.h"


void swap_edge(Edge* a, Edge* b) {
    Edge temp = *a;
    *a = *b;
    *b = temp;
}

void print_heap(const Edge* heap, int size) {
    printf("Heap contents (size = %d):\n", size);
    for (int i = 0; i < size; i++) {
        printf("  [%d] node = %d, cost = %.4f\n", i, heap[i].node, heap[i].cost);
    }
    printf("\n");
}

void heapify_up(Edge* heap, int index) {
    while (index > 0) {
        int parent = (index - 1) / 2;
        if (heap[index].cost > heap[parent].cost) {
            swap_edge(&heap[index], &heap[parent]);
            index = parent;
        } else {
            break;
        }
    }
}

void heapify_down(Edge* heap, int size, int index) {
    while (2 * index + 1 < size) {
        int left = 2 * index + 1;
        int right = 2 * index + 2;
        int largest = index;
        if (left < size && heap[left].cost > heap[largest].cost)
            largest = left;
        if (right < size && heap[right].cost > heap[largest].cost)
            largest = right;
        if (largest != index) {
            swap_edge(&heap[index], &heap[largest]);
            index = largest;
        } else {
            break;
        }
    }
}

void insert_edge(Edge* heap, int* size, int k, Edge new_edge) {
    if (*size < k) {
        heap[*size] = new_edge;
        heapify_up(heap, *size);
        (*size)++;
    } else if (new_edge.cost < heap[0].cost) {
        heap[0] = new_edge;
        heapify_down(heap, *size, 0);
    }
}
