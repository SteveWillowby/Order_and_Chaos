#include "edge.h"
#include "sparse_graph.h"

#include<unordered_set>

#ifndef SYM__SIMULATED_ANNEALING_SEARCH_H
#define SYM__SIMULATED_ANNEALING_SEARCH_H

// Given a graph g, does a search to find
std::unordered_set<Edge, EdgeHash>
                    simulated_annealing_search(const SparseGraph& g,
                                               size_t num_iterations);

#endif
