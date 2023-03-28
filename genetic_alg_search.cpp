#include "edge.h"
#include "genetic_alg_search.h"
#include "graph.h"

#include<rand>
#include<thread>
#include<unordered_set>
#include<vector>

void graph_score_worker() {
    while (true) {
        
    }
}

// Given a graph `g`, does a search to find good candidates for noise.
//
// Lists the top `k` candidates with their associated scores.
//
// `nt` is the number of threads to be used by the code.
std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                                 genetic_alg_search(const Graph& g,
                                                    size_t num_iterations,
                                                    size_t k,
                                                    size_t nt) {
    if (nt == 0) {
        nt = std::thread::hardware_concurrency();
    }

    std::vector<std::thread> thread_pool;
    for (size_t i = 0; i < nt; i++) {
        thread_pool.push_back(std::thread(thread_pool_worker));
    }
}
