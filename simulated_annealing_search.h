#include "edge.h"
#include "sparse_graph.h"

#include<unordered_set>
#include<vector>

#ifndef SYM__SIMULATED_ANNEALING_SEARCH_H
#define SYM__SIMULATED_ANNEALING_SEARCH_H

// Given a graph g, does a search to find good candidates for noise.
//
// Lists the top k candidates with their associated scores.
std::vector<std::pair<std::unordered_set<Edge, EdgeHash>, long double>>>
                            simulated_annealing_search(const SparseGraph& g,
                                                       size_t num_iterations,
                                                       size_t k);

#endif
