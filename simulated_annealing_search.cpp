#include "coloring.h"
#include "edge.h"
#include "nauty_traces.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"
#include "sparse_graph.h"

#include<cmath>
#include<random>
#include<unordered_set>
#include<vector>

std::vector<std::pair<std::unordered_set<Edge, EdgeHash>, long double>>>
                            simulated_annealing_search(const SparseGraph& g,
                                                       size_t num_iterations,
                                                       size_t k) {
    NTSparseGraph g_main = g;

    bool use_self_loops = g_main.num_loops() > 0;

    size_t max_possible_edges = (size_t(use_self_loops) * g.num_nodes()) +
            (g.num_nodes() * (g.num_nodes() - 1)) / (size_t(!g.directed) + 1);


    NautyTracesOptions o;
    o.get_node_orbits = true;
    o.get_edge_orbits = true;
    o.get_canonical_node_order = false;

    // Perform preliminary calculations to speed up the code.
    NautyTracesResult orbits_info = traces(g_main, o);
    
    Coloring<int> node_orbits& = orbits_info.node_orbits;
    Coloring<Edge, EdgeHash>& edge_orbits = orbits_info.edge_orbits;

    // Get a score for the initial candidate noise (i.e. no noise).
    long double prev_score = score(g_main, node_orbits, edge_orbits,
                                   editable_edge_orbits,
                                   candidate_additions, candidate_removals);
    long double curr_score;

    std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                          long double>>> top_k_results =
        std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                              long double>>>();
    for (int i = 0; i < k; i++) {
        // initialize top_k_results with empty edge-flip sets and their score
        top_k_results.push_back(
            std::pair<std::unordered_set<Edge, EdgeHash>, long double>(
                std::unordered_set<Edge, EdgeHash>(), prev_score));
    }

    if (g.num_edges() == max_possible_edges || g.num_edges() == 0) {
        // It's a clique or an empty graph. There's nothing to look for.
        return top_k_results;
    }

    // Create an editable copy of the edge orbits.
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

    std::random_device rd;  // Will provide a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // chance of adding a new (non-)edge to the candidate
    float prob_new = 0.5;
    // chance of choosing an edge vs. a non-edge in when adding or removing
    //  something from the candidate set.
    float prob_edge = float(g.num_edges()) / float(max_possible_edges);

    float rand;

    size_t num_non_edges = max_possible_edges - g_main.num_edges();
    size_t max_flip_or_edge =
            4 * (num_non_edges < g.num_edges() ? num_non_edges : g.num_edges());
    // O(max_possible_edges)   ---- O(well)
    __combinatoric_utility.set_max_access(max_possible_edges, max_flip_or_edge);

    bool add, is_edge;
    Edge flipped;

    size_t num_skipped_iterations = 0;

    float temperature, transition_prob;

    std::unordered_set<Edge, EdgeHash> flips;

    for (size_t itr = 0; itr < num_iterations + num_skipped_iterations; itr++) {
        // Main loop.
        add = dist(gen) < prob_new;
        is_edge = dist(gen) < prob_edge;
        if (add) {
            // Add an edge or a non-edge to the candidate set.
            if (is_edge) {
                // Add an edge to the candidate set.

                // TODO: Add the sampler code.
                flipped = sampler.sample_edge();
                g_main.add_edge(flipped->first, flipped->second);

                if (g_main.num_edges() == 0) {
                    // No edges to remove from the graph.
                    num_skipped_iterations++;
                    continue;
                }
            } else {
                // Add a non-edge to the candidate set.
                if (g_main.num_edges() == max_possible_edges) {
                    // No edges to add to the graph.
                    num_skipped_iterations++;
                    continue;
                }
            }
        } else {
            // Remove an edge or a non-edge from the candidate set.
            if (is_edge) {
                // Remove an edge from the candidate set.
                if (candidate_removals.size() == 0) {
                    // There is no edge in the candidate set.
                    num_skipped_iterations++;
                    continue;
                }
            } else {
                // Remove a non-edge from the candidate set.
                if (candidate_additions.size() == 0) {
                    // There is no non-edge in the candidate set.
                    num_skipped_iterations++;
                    continue;
                }
            }
        }

        // Currently the temperature decreases linearly.
        temperature =
            1.0 - (double(itr + 1 - num_skipped_iterations) / num_iterations);

        curr_score = score(g_main, node_orbits, edge_orbits,
                           editable_edge_orbits,
                           candidate_additions, candidate_removals);
        if (curr_score < prev_score) {
            // Consider rejecting the change.
            transition_prob = std::exp2((curr_score - prev_score) /temperature);
            if (dist(gen) >= transition_prob) {
                // Reject the change.
                // TODO: Undo the change.
                continue;
            }
        }

        // See if this goes in the top-k score list.
        if (curr_score > top_k_results[0][1]) {
            size_t i;
            for (i = 0; i < k; i++) {
                if (curr_score < top_k_results[i][1]) {
                    break;
                }
            }
            i--;
            for (size_t j = 0; j < i; j++) {
                top_k_results[j] = top_k_results[j + 1];
            }

            // Combine additions and removals.
            flips = std::unordered_set<Edge, EdgeHash>(candidate_additions);
            for (auto e_itr = candidate_removals.begin();
                      e_itr != candidate_removals.end(); e_itr++) {
                flips.insert(*e_itr);
            }
            top_k_results[i] = std::pair<std::unordered_set<Edge, EdgeHash>,
                                         long double>(flips, curr_score);
        }
        prev_score = curr_score;
    }

    return top_k_results;
}
