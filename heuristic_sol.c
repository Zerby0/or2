#include "tsp.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>

static void swap(int* a, int pos1, int pos2) {
    int temp = a[pos1];
    a[pos1] = a[pos2];
    a[pos2] = temp;
}

void nearest_neighbor_from(instance* inst, int start) {
    int tour[inst->num_nodes];
    for(int i = 0; i < inst->num_nodes; i++) {
        tour[i] = i;
    }
    swap(tour, 0, start);
    double tot_cost = 0;
    for (int i = 0; i < inst->num_nodes - 1; i++) {
        int min_index = -1;
        double min_cost = INF_COST;
        for (int j = i + 1; j < inst->num_nodes; j++) {
            double cost = inst->costs_array[tour[i] * inst->num_nodes + tour[j]];
            if (cost < min_cost) {
                min_cost = cost;
                min_index = j;
            }
        }
		debug(80, "[%d] extending %d -> %d [%f]\n", i, tour[i], tour[min_index], min_cost);
        swap(tour, i + 1, min_index);
        tot_cost += min_cost;
		if (inst->verbose >= 90)
			plot_partial_sol(inst, tour, i + 2);
    }
	int last_node = tour[inst->num_nodes - 1];
    tot_cost += inst->costs_array[last_node * inst->num_nodes + start];
    update_sol(inst, tour, tot_cost);
}

void nearest_neighbor(instance* inst) {
	nearest_neighbor_from(inst, 0);
}

// inserts into `arr` currently of length `len` at position `pos` the value `value` in O(len)
void array_insert(int* arr, int len, int pos, int value) {
	memmove(arr + pos + 1, arr + pos, sizeof(int) * (len - pos));
	arr[pos] = value;
}

double find_max_segment(const instance* inst, int* max_start, int* max_end) {
	double max_dist = 0;
	*max_start = *max_end = -1;
	for (int i = 0; i < inst->num_nodes; i++) {
		for (int j = i + 1; j < inst->num_nodes; j++) {
			double dist = inst->costs[i][j];
			if (dist > max_dist) {
				max_dist = dist;
				*max_start = i;
				*max_end = j;
			}
		}
	}
	return max_dist;
}

/*
    Not extra point if we implement this!
    Extra-mialage heuristic: choose 2 point (consider far away) and connect them twice (to create a cycle) and then connect the other points to the cycle 
    rewritting the cycle (greedy policy avoiding the points already in the cycle) and then return the solution
    the time complexity is O(n^3) because you have to check every node for every node

*/
void extra_milage(instance* inst) {
	int max_start, max_end;
	double max_dist = find_max_segment(inst, &max_start, &max_end);
	int tour[inst->num_nodes + 1];
	tour[0] = max_start;
	tour[1] = max_end;
	tour[2] = max_start;
	double tot_cost = 2 * max_dist;
	int count = 2; // number of nodes in the tour, excluding the starting point at the end
	bool used[inst->num_nodes];
	memset(used, false, sizeof(bool) * inst->num_nodes);
	used[max_start] = used[max_end] = true;
	while (count < inst->num_nodes) {
		// for node to insert in the tour
		double min_cost = INFINITY;
		int min_h = -1, min_pos = -1;
		for (int h = 0; h < inst->num_nodes; h++) {
			if (used[h]) continue;
			// for each segment to be extended
			for (int pos = 0; pos < count; pos++) {
				int i = tour[pos], j = tour[pos + 1];
				assert(i != j && i != h && j != h);
				double cij = inst->costs[i][j];
				double cihj = inst->costs[i][h] + inst->costs[h][j];
				double extra = cihj - cij;
				assert(extra >= -1e-9); // triangle inequality
				if (extra < min_cost) {
					min_cost = extra;
					min_h = h;
					min_pos = pos;
				}
			}
		}
		debug(80, "extending tour with %d -> %d -> %d [+%f]\n", tour[min_pos], min_h, tour[min_pos + 1], min_cost);
		array_insert(tour, count + 1, min_pos + 1, min_h);
		tot_cost += min_cost;
		count++;
		used[min_h] = true;
		if (inst->verbose >= 90)
			plot_partial_sol(inst, tour, count);
	}
	update_sol(inst, tour, tot_cost);
}


/*
    TODO: 2-opt
    BEST CARD OF THE PROF: Refinement method (2-opt) -> take a solution and try to improve it
    2-opt: take 2 edges and swap them it's called 2-opt move, if the cost is better (computing the delta cost. if delta positive -> cost is increased by delta so 
    no better sol, if negative totalcost + delta < totalcost -> better sol) then keep the solution 
    the method become good only with 2-opt, error of 2% of the optimal solution (we can sold the code ahah)
*/
