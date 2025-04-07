#include "tsp.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <string.h>

double nearest_neighbor_from(const Instance* inst, int* tour, int start) {
    for(int i = 0; i < inst->num_nodes; i++) {
        tour[i] = i;
    }
    swap(tour, 0, start);
    double tot_cost = 0;
    for (int i = 0; i < inst->num_nodes - 1; i++) {
        int min_index = -1;
        double min_cost = INF_COST;
        for (int j = i + 1; j < inst->num_nodes; j++) {
            double cost = get_cost(inst, tour[i], tour[j]);
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
    tot_cost += get_cost(inst, last_node, start);
	return tot_cost;
}

void nearest_neighbor(Instance* inst) {
	int start = rand() % inst->num_nodes;
	int tour[inst->num_nodes];
	double c = nearest_neighbor_from(inst, tour, start);
	update_sol(inst, tour, c);
}

void nearest_neighbor_all_starts(Instance* inst) {
	int tour[inst->num_nodes];
	int s0 = rand() % inst->num_nodes;
	for (int delta = 0; delta < inst->num_nodes && !is_out_of_time(inst); delta++) {
		int start = (s0 + delta) % inst->num_nodes;
		double c = nearest_neighbor_from(inst, tour, start);
		if (inst->two_opt)
			two_opt_from(inst, tour, &c);
		debug(50, "NN from %d: cost = %f\n", start, c);
		update_sol(inst, tour, c);
	}
}

// inserts into `arr` currently of length `len` at position `pos` the value `value` in O(len)
void array_insert(int* arr, int len, int pos, int value) {
	memmove(arr + pos + 1, arr + pos, sizeof(int) * (len - pos));
	arr[pos] = value;
}

double find_max_segment(const Instance* inst, int* max_start, int* max_end) {
	double max_dist = 0;
	*max_start = *max_end = -1;
	for (int i = 0; i < inst->num_nodes; i++) {
		for (int j = i + 1; j < inst->num_nodes; j++) {
			double dist = get_cost(inst, i, j);
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
void extra_milage(Instance* inst) {
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
	while (count < inst->num_nodes && !is_out_of_time(inst)) {
		// for node to insert in the tour
		double min_cost = INFINITY;
		int min_h = -1, min_pos = -1;
		for (int h = 0; h < inst->num_nodes; h++) {
			if (used[h]) continue;
			// for each segment to be extended
			for (int pos = 0; pos < count; pos++) {
				int i = tour[pos], j = tour[pos + 1];
				assert(i != j && i != h && j != h);
				double cij = get_cost(inst, i, j);
				double cihj = get_cost(inst, i, h) + get_cost(inst, h, j);
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
	if (is_out_of_time(inst) && count < inst->num_nodes) {
		debug(30, "Extra-milage: out of time before completing the tour, completing with a random solution\n");
		for (int h = 0; h < inst->num_nodes; h++) {
			if (used[h]) continue;
			tour[count++] = h;
			used[h] = true;
		}
		tot_cost = compute_tour_cost(inst, tour);
		assert(count == inst->num_nodes);
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
