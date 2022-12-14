#include "coloring.h"
#include "edge.h"
#include "edge_sampler.h"
#include "nauty_traces.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"
#include "graph.h"

#include<cmath>
#include<random>
#include<tuple>
#include<unordered_set>
#include<vector>

size_t hash_edge_set(const std::unordered_set<Edge, EdgeHash>& edges,
                     const EdgeHash& edge_hasher);
bool edge_set_eq(const std::unordered_set<Edge, EdgeHash>& A,
                 const std::unordered_set<Edge, EdgeHash>& B);

// TODO: Move the graph and edge-coloring maintenance out of scoring_function
//  and into the simulated annealing search.
// Better yet, make two different scoring functions. The current one may be
//  useful for things like genetic algorithms, etc.

std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                            simulated_annealing_search(const Graph& g,
                                                       size_t num_iterations,
                                                       size_t k) {
    NTSparseGraph g_main = g;
    EdgeHash edge_hasher;

    bool use_self_loops = g.num_loops() > 0;

    size_t max_possible_edges = (size_t(use_self_loops) * g.num_nodes()) +
            (g.num_nodes() * (g.num_nodes() - 1)) / (size_t(!g.directed) + 1);

    size_t num_non_edges = max_possible_edges - g_main.num_edges();
    size_t max_flip_or_edge =
            4 * (num_non_edges < g.num_edges() ? num_non_edges : g.num_edges());
    // O(max_possible_edges)   ---- O(well)
    CombinatoricUtility comb_util(max_possible_edges, max_flip_or_edge);


    NautyTracesOptions o;
    o.get_node_orbits = true;
    o.get_edge_orbits = true;
    o.get_canonical_node_order = false;

    // Perform preliminary calculations to speed up the code.
    NautyTracesResults orbits_info = traces(g_main, o);

    // Create an editable copy of the edge orbits.
    Coloring<Edge, EdgeHash> editable_edge_orbits = Coloring<Edge, EdgeHash>();
    for (auto c_itr = orbits_info.edge_orbits.colors().begin();
              c_itr != orbits_info.edge_orbits.colors().end(); c_itr++) {
        for (auto cell_itr = orbits_info.edge_orbits.cell(*c_itr).begin();
                  cell_itr != orbits_info.edge_orbits.cell(*c_itr).end();
                  cell_itr++) {
            editable_edge_orbits.set(*cell_itr, *c_itr);
        }
    }

    std::unordered_set<Edge, EdgeHash> candidate_additions =
                                    std::unordered_set<Edge, EdgeHash>();
    std::unordered_set<Edge, EdgeHash> candidate_removals =
                                    std::unordered_set<Edge, EdgeHash>();

    // Get a score for the initial candidate noise (i.e. no noise).
    long double prev_score = score(g_main, comb_util, orbits_info.node_orbits,
                                   orbits_info.edge_orbits,
                                   editable_edge_orbits,
                                   candidate_additions, candidate_removals);
    long double curr_score;

    std::vector<std::tuple<std::unordered_set<Edge, EdgeHash>,
                           long double, size_t>> top_k_results =
        std::vector<std::tuple<std::unordered_set<Edge, EdgeHash>,
                               long double, size_t>>();
    std::vector<std::pair<std::unordered_set<Edge, EdgeHash>, long double>>
        return_value =
            std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                                  long double>>();

    for (size_t i = 0; i < k; i++) {
        // initialize top_k_results with empty edge-flip sets and their score
        top_k_results.push_back(
            std::make_tuple<std::unordered_set<Edge, EdgeHash>,
                            long double, size_t>(
                                    std::unordered_set<Edge, EdgeHash>(),
                                    (long double)(prev_score), 0));
    }

    if (g.num_edges() == max_possible_edges || g.num_edges() == 0) {
        // It's a clique or an empty graph. There's nothing to look for.
        for (size_t i = 0; i < k; i++) {
            return_value.push_back(std::pair<std::unordered_set<Edge, EdgeHash>,
                                             long double>(
                                                std::get<0>(top_k_results[i]),
                                                std::get<1>(top_k_results[i])));
        }
        return return_value;
    }

    std::random_device rd;  // Will provide a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    // chance of adding a new (non-)edge to the candidate
    float prob_new = 0.5;
    // chance of choosing an edge vs. a non-edge in when adding or removing
    //  something from the candidate set.
    float prob_edge = float(g.num_edges()) / float(max_possible_edges);

    bool add, is_edge;
    Edge flipped;

    size_t num_skipped_iterations = 0;

    float temperature, transition_prob;

    std::unordered_set<Edge, EdgeHash> flips;

    EdgeSampler sampler(g, gen);

    for (size_t itr = 0; itr < num_iterations + num_skipped_iterations; itr++) {
        // Main loop.
        add = dist(gen) < prob_new;
        is_edge = dist(gen) < prob_edge;
        if (add) {
            // Add an edge or a non-edge to the candidate set.
            if (is_edge) {
                // Add an edge to the candidate removals set.
                //  I.e. Remove the edge from the graph.
                if (g.num_edges() == candidate_removals.size()) {
                    // No edges to remove from the graph.
                    num_skipped_iterations++;
                    continue;
                }

                flipped = sampler.sample_edge();
                // g_main.delete_edge(flipped->first, flipped->second);
                candidate_removals.insert(flipped);
            } else {
                // Add a non-edge to the candidate additions set.
                //  I.e. Add a (non-)edge to the graph as an edge.
                if (max_possible_edges - g.num_edges() ==
                                    candidate_additions.size()) {
                    // No edges to add to the graph.
                    num_skipped_iterations++;
                    continue;
                }

                flipped = sampler.sample_non_edge();
                // g_main.add_edge(flipped->first, flipped->second);
                candidate_additions.insert(flipped);
            }
        } else {
            // Remove an edge or a non-edge from the candidate set.
            if (is_edge) {
                // Remove an edge from the candidate removals set.
                //  I.e. Put the edge back in the graph.
                if (candidate_removals.size() == 0) {
                    // There is no edge in the candidate set.
                    num_skipped_iterations++;
                    continue;
                }

                flipped = sampler.un_sample_edge();
                // g_main.add_edge(flipped->first, flipped->second);
                candidate_removals.erase(flipped);
            } else {
                // Remove a non-edge from the candidate additions set.
                //  I.e. Make it a non-edge again.
                if (candidate_additions.size() == 0) {
                    // There is no non-edge in the candidate set.
                    num_skipped_iterations++;
                    continue;
                }

                flipped = sampler.un_sample_non_edge();
                // g_main.delete_edge(flipped->first, flipped->second);
                candidate_additions.erase(flipped);
            }
        }

        // Currently the temperature decreases linearly.
        temperature =
            1.0 - (double(itr + 1 - num_skipped_iterations) / num_iterations);

        curr_score = score(g_main, comb_util, orbits_info.node_orbits,
                           orbits_info.edge_orbits,
                           editable_edge_orbits,
                           candidate_additions, candidate_removals);
        if (curr_score < prev_score) {
            // Consider rejecting the change.
            transition_prob = std::exp2((curr_score - prev_score) /temperature);
            if (dist(gen) >= transition_prob) {
                // Reject the change.
                if (add) {
                    if (is_edge) {
                        candidate_removals.erase(flipped);
                    } else {
                        candidate_additions.erase(flipped);
                    }
                } else {
                    if (is_edge) {
                        candidate_removals.insert(flipped);
                    } else {
                        candidate_additions.insert(flipped);
                    }
                }
                sampler.undo();
                continue;
            }
        }

        // See if this goes in the top-k score list.
        if (curr_score > std::get<1>(top_k_results[0])) {

            // Combine additions and removals.
            flips = std::unordered_set<Edge, EdgeHash>(candidate_additions);
            for (auto e_itr = candidate_removals.begin();
                      e_itr != candidate_removals.end(); e_itr++) {
                flips.insert(*e_itr);
            }
            size_t hash_value = hash_edge_set(flips, edge_hasher);

            // Verify that this edge set is not already in top_k_results.
            bool is_new = true;
            for (size_t i = 1; i < k; i++) {
                if (std::get<2>(top_k_results[i]) == hash_value &&
                        edge_set_eq(std::get<0>(top_k_results[i]), flips)) {
                    is_new = false;
                    break;
                }
            }
            if (!is_new) {
                // We still move to the `flips` edge set, but we do not update
                //  top_k_results.
                prev_score = curr_score;
                continue;
            }

            size_t i;
            for (i = 0; i < k; i++) {
                if (curr_score < std::get<1>(top_k_results[i])) {
                    break;
                }
            }
            i--;
            for (size_t j = 0; j < i; j++) {
                top_k_results[j] = top_k_results[j + 1];
            }
            top_k_results[i] =
               std::make_tuple<std::unordered_set<Edge, EdgeHash>,
                               long double, size_t>(
                                 std::unordered_set<Edge, EdgeHash>(flips),
                                 (long double)(curr_score), size_t(hash_value));
        }
        prev_score = curr_score;
    }

    for (size_t i = 0; i < k; i++) {
        return_value.push_back(std::pair<std::unordered_set<Edge, EdgeHash>,
                                         long double>(
                                                std::get<0>(top_k_results[i]),
                                                std::get<1>(top_k_results[i])));
    }
    return return_value;
}


size_t hash_edge_set(const std::unordered_set<Edge, EdgeHash>& edges,
                     const EdgeHash& edge_hasher) {
    size_t h = 0;
    // X-OR all the edge hashes.
    for (auto e_itr = edges.begin(); e_itr != edges.end(); e_itr++) {
        h ^= edge_hasher(*e_itr);
    }
    return h;
}

bool edge_set_eq(const std::unordered_set<Edge, EdgeHash>& A,
                 const std::unordered_set<Edge, EdgeHash>& B) {
    if (A.size() != B.size()) {
        return false;
    }
    for (auto e_itr = A.begin(); e_itr != A.end(); e_itr++) {
        if (B.find(*e_itr) == B.end()) {
            return false;
        }
    }
    return true;
}
