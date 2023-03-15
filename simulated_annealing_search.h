#include "edge.h"
#include "graph.h"

#include<unordered_set>
#include<vector>

#ifndef SYM__SIMULATED_ANNEALING_SEARCH_H
#define SYM__SIMULATED_ANNEALING_SEARCH_H

// Given a graph `g`, does a search to find good candidates for noise.
//
// Lists the top `k` candidates with their associated scores.
//
// `expected_additions` is the number of edges we expect would be
//  added to the HYPOTHESIS graph (i.e. removed from `g`).
// `expected_removals` is the same but in reverse.
std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                         simulated_annealing_search(const Graph& g,
                                                    size_t num_iterations,
                                                    size_t k);
                                                    // size_t expected_additions,
                                                    // size_t expected_removals);
#endif
