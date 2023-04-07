#include<cmath>
#include<iostream>
#include<random>
#include<string>

#include "basic_edge_set.h"
#include "edge.h"
#include "edge_sampler.h"
#include "file_utils.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"
#include "scoring_function.h"
#include "thread_pool_scorer.h"

int main(void) {
    const bool directed = true;
    std::cout<<"## Directed?   "<<(directed ? "Yes" : "No")<<std::endl;
    const size_t num_nodes = 5;
    const size_t num_edges = 5;
    const size_t num_additions = 1;
    const size_t num_deletions = 1;
    const size_t num_noise_sets = 8;
    const size_t num_rounds = 5;
    const size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed));
    const size_t max_flip_or_edge = num_edges * 2;

    NTSparseGraph g(directed, num_nodes);

    std::random_device rd;  // Will provide a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    std::cout<<"The graph has "<<num_nodes<<" nodes and the following edges:"
             <<std::endl;
    EdgeSampler creation_sampler(g, gen);
    for (size_t i = 0; i < num_edges; i++) {
        Edge e = creation_sampler.sample_non_edge();
        std::cout<<" ("<<e.first<<", "<<e.second<<"),";
        g.add_edge(e.first, e.second);
    }
    std::cout<<std::endl;

    EdgeSampler noise_sampler(g, gen);

    CombinatoricUtility comb_util(max_possible_edges, max_flip_or_edge);

    
    NautyTracesOptions nt_options;
    nt_options.get_node_orbits = true;
    nt_options.get_edge_orbits = true;
    nt_options.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g, nt_options);

    // Calculate the noise probability at which the full graph is equally
    //  likely to be the noise as it is to be the structure.
    long double alpha_plus =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(g.num_nodes()))) /
                                (long double)(g.num_edges()));
    long double log2_p_plus = std::log2l(alpha_plus) -
                              std::log2l(1.0 + alpha_plus);
    long double log2_1_minus_p_plus = -std::log2l(1.0 + alpha_plus);

    long double alpha_minus =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(g.num_nodes()))) /
                        (long double)(max_possible_edges - g.num_edges()));
    long double log2_p_minus = std::log2l(alpha_minus) -
                               std::log2l(1.0 + alpha_minus);
    long double log2_1_minus_p_minus = -std::log2l(1.0 + alpha_minus);

    std::cout<<"Creating Thread Pool Scorer..."<<std::endl;

    ThreadPoolScorer TPS(2, g, comb_util,
                         nt_result.node_orbits, nt_result.edge_orbits,
                         log2_p_plus, log2_p_minus,
                         log2_1_minus_p_plus, log2_1_minus_p_minus);

    std::cout<<" ...Finished Creating Thread Pool Scorer"<<std::endl;
    
    std::vector<std::pair<std::unique_ptr<EdgeSet>,
                          std::unique_ptr<EdgeSet>>> tasks =
        std::vector<std::pair<std::unique_ptr<EdgeSet>,
                              std::unique_ptr<EdgeSet>>>();

    std::unordered_set<Edge, EdgeHash> edge_additions, edge_deletions;
    for (size_t r = 0; r < num_rounds; r++) {
        tasks.clear();
        for (size_t i = 0; i < num_noise_sets; i++) {
            edge_additions = std::unordered_set<Edge, EdgeHash>();
            edge_deletions = std::unordered_set<Edge, EdgeHash>();
            for (size_t j = 0; j < num_additions; j++) {
                edge_additions.insert(noise_sampler.sample_non_edge());
            }
            for (size_t j = 0; j < num_additions; j++) {
                noise_sampler.un_sample_non_edge();
            }
            for (size_t j = 0; j < num_deletions; j++) {
                edge_deletions.insert(noise_sampler.sample_edge());
            }
            for (size_t j = 0; j < num_deletions; j++) {
                noise_sampler.un_sample_edge();
            }
            tasks.push_back(std::pair<std::unique_ptr<EdgeSet>,
                                      std::unique_ptr<EdgeSet>>(
                   std::unique_ptr<EdgeSet>(new BasicEdgeSet(edge_additions)),
                   std::unique_ptr<EdgeSet>(new BasicEdgeSet(edge_deletions))));
        }

        const std::vector<long double>& scores = TPS.get_scores(&tasks);

        for (size_t i = 0; i < num_noise_sets; i++) {
            std::cout<<"- ";
            for (auto e = tasks[i].first->edges().begin();
                      e != tasks[i].first->edges().end(); e++) {
                std::cout<<"("<<e->first<<", "<<e->second<<"), ";
            }
            std::cout<<" \t + ";
            for (auto e = tasks[i].second->edges().begin();
                      e != tasks[i].second->edges().end(); e++) {
                std::cout<<"("<<e->first<<", "<<e->second<<"), ";
            }
            std::cout<<" \t "<<scores[i]<<std::endl;
        }
    }

    return 0;
}
