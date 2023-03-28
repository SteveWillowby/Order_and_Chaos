#include "edge.h"
#include "graph.h"

#include<unordered_set>
#include<vector>

#ifndef SYM__GENETIC_ALG_SEARCH_H
#define SYM__GENETIC_ALG_SEARCH_H

// Given a graph `g`, does a search to find good candidates for noise.
//
// Lists the top `k` candidates with their associated scores.
//
// `nt` is the number of threads to be used by the code. Pass 0 to use
//      a default value of std::threads::hardware_concurrency();
std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                                 genetic_alg_search(const Graph& g,
                                                    size_t num_iterations,
                                                    size_t k,
                                                    size_t nt);

#endif
