#include<cmath>
#include<iostream>
#include<random>
#include<string>

#include "nt_wrappers/nauty_traces.h"
#include "scheno/scheno.h"

int main(void) {

    const size_t num_threads = 1;

    const bool directed = false;
    std::cout<<"## Directed?   "<<(directed ? "Yes" : "No")<<std::endl;
    const size_t num_nodes = 5;
    const size_t num_edges = 5;
    const size_t num_additions = 1;
    const size_t num_deletions = 1;
    const size_t num_noise_sets = 1;
    const size_t num_rounds = 18;
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

    IntEdgeConverterAndSampler iecas(g);
    const std::vector<long double>& heuristics = iecas.get_heuristic_scores();
    for (size_t i = 0; i < heuristics.size(); i++) {
        if (heuristics[i] == 0.0) {
            continue;
        }
        Edge e = iecas.edge(i);
        std::cout<<"("<<e.first<<", "<<e.second<<"): "<<heuristics[i]
                 <<std::endl;
    }

    std::cout<<"Creating Thread Pool Scorer..."<<std::endl;

    ThreadPoolScorer TPS(num_threads, g, comb_util,
                         nt_result.node_orbits, nt_result.edge_orbits,
                         log2_p_plus, log2_p_minus,
                         log2_1_minus_p_plus, log2_1_minus_p_minus,
                         (size_t) -1, true, true);

    std::cout<<" ...Finished Creating Thread Pool Scorer"<<std::endl;
    
    std::vector<std::unique_ptr<EdgeSetPair>> tasks =
        std::vector<std::unique_ptr<EdgeSetPair>>();

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
            tasks.push_back(std::unique_ptr<EdgeSetPair>(
                        new BasicEdgeSetPair(edge_deletions, edge_additions)));
        }

        const std::vector<std::pair<long double, long double>>& scores =
                    TPS.get_scores(&tasks);

        for (size_t i = 0; i < num_noise_sets; i++) {
            std::cout<<"- ";
            std::pair<std::unordered_set<Edge, EdgeHash>,
                      std::unordered_set<Edge, EdgeHash>> r =
                            tasks[i]->edges_and_non_edges();
            for (auto e = r.first.begin();
                      e != r.first.end(); e++) {
                std::cout<<"("<<e->first<<", "<<e->second<<"), ";
            }
            std::cout<<" \t + ";
            for (auto e = r.second.begin();
                      e != r.second.end(); e++) {
                std::cout<<"("<<e->first<<", "<<e->second<<"), ";
            }
            std::cout<<" \t "<<scores[i].first<<", "
                             <<scores[i].second<<std::endl;
        }
    }

    return 0;
}
