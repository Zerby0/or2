#ifndef MAX_HEAP_H
#define MAX_HEAP_H

typedef struct {
    int node;       //destination node
    double cost;    //cost of the edge
} Edge;

void heapify_up(Edge* heap, int index);
void heapify_down(Edge* heap, int size, int index);
void insert_edge(Edge* heap, int* size, int k, Edge new_edge);
void swap_edge(Edge* a, Edge* b);
void print_heap(const Edge* heap, int size);
#endif /* MAX_HEAP_H */