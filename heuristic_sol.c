#include "tsp.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>

/*  idea: Lesson 4/03/2025

	nearest neighbor

    you have a set of cities, start from a given city
    find the *nearest* city to the current city
    add the city to the tour
    repeat until all cities are visited
    return to the starting city
    the tour is the solution
    essentially is a for loop that goes through all the cities and finds the nearest one
    use an array to store every node, then put on first the starting node switching the nodes, then check the distance between the 
    first and the other nodes on the other side of the array (unvisited nodes) then put in second place the second nearest node, and so on
    updating the total cost of the solution every time you add a node
    the time complexity is O(n^2) because you have to check every node for every node
    call the function with all the nodes as starting point and take the best solution (policy called multi-start)
    Cost of an edge that goes from node i to node j cost(i,j,inst) = sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2) but it's time consuming
    so i precomputing the cost of every edge and store it in a matrix
    so after reading the file i can compute the cost of every edge and store it in a matrix calling compute_all_costs(inst) adding a new field in the structure
    double** costs;
    then i can call the function cost(i,j,inst) = inst->costs[i][j]
    if i want to put all in array i can use a single array of size n*n and access the element costs[i*n+j] where n = inst->num_nodes (this is faster)
    time limit: we can have t_start = world clock seconds, not clock time do not check the time limit on greedy, put on specific points, after the greedy
    if the time limit is exceeded, return the best solution found so far.
    time: 60s or 60*10s
*/
void nearest_neighbor(instance* inst) {
	int tour[inst->num_nodes];
	double tot_cost = 0;
	int count = 1, prev = 0;
	tour[0] = prev;
	bool used[inst->num_nodes];
	memset(used, false, sizeof(bool) * inst->num_nodes);
	while (count < inst->num_nodes) {
		double min_cost = INFINITY;
		int min_index = -1;
		for (int i = 0; i < inst->num_nodes; i++) {
			if (used[i]) continue;
			double cost = inst->costs[prev][i];
			if (cost < min_cost) {
				min_cost = cost;
				min_index = i;
			}
		}
		debug(80, "extending %d -> %d [%f]\n", prev, min_index, min_cost);
		used[min_index] = true;
		tour[count++] = min_index;
		tot_cost += min_cost;
		prev = min_index;
		if (inst->verbose >= 90)
			plot_partial_sol(inst, tour, count);
	}
	tot_cost += inst->costs[prev][0];
	update_sol(inst, tour, tot_cost);
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
