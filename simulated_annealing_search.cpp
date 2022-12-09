#include "coloring.h"
#include "edge.h"
#include "nauty_traces.h"
#include "nt_sparse_graph.h"
#include "sparse_graph.h"

#include<unordered_set>
#include<vector>

std::vector<std::pair<std::unordered_set<Edge, EdgeHash>, long double>>>
                            simulated_annealing_search(const SparseGraph& g,
                                                       size_t num_iterations,
                                                       size_t k) {
    NTSparseGraph g_main = g;

    NautyTracesOptions o;
    o.get_node_orbits = true;
    o.get_edge_orbits = true;
    o.get_canonical_node_order = false;

    // Perform preliminary calculations to speed up the code.
    NautyTracesResult orbits_info = traces(g_main, o);
    
    Coloring<int> node_orbits& = orbits_info.node_orbits;
    Coloring<Edge, EdgeHash>& edge_orbits = orbits_info.edge_orbits;

    // Create a copy of the edge orbits.
    Coloring<Edge, EdgeHash> editable_edge_orbits = Coloring<Edge, EdgeHash>();
    for (auto c_itr = edge_orbits.colors().begin();
              c_itr != edge_orbits.colors().end(); c_itr++) {
        for (auto cell_itr = edge_orbits.cell(*c_itr).begin();
                  cell_itr != edge_orbits.cell(*c_itr).end(); cell_itr++) {
            editable_edge_orbits.set(*cell_itr, *c_itr);
        }
    }

    std::unordered_set<Edge, EdgeHash> candidate_additions =
                                    std::unordered_set<Edge, EdgeHash>();
    std::unordered_set<Edge, EdgeHash> candidate_removals =
                                    std::unordered_set<Edge, EdgeHash>();

    // Get a score for the initial candidate noise.
    long double prev_score = score(g_main, node_orbits, edge_orbits,
                                   editable_edge_orbits,
                                   candidate_additions, candidate_removals);
    long double curr_score;

    std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                          long double>>> top_k_results =
        std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                              long double>>>();
    for (int i = 0; i < k; i++) {
        // initialize top_k_results
    }

    for (size_t itr = 0; itr < num_iterations; itr++) {
        // Main loop.
    }
}
