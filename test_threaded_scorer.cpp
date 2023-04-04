#include<cmath>
#include<iostream>
#include<random>
#include<string>

#include "edge.h"
#include "edge_sampler.h"
#include "file_utils.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"
#include "scoring_function.h"
#include "thread_pool_scorer.h"

int main(void) {
    const bool directed = false;
    const size_t num_nodes = 7;
    const size_t num_edges = 15;
    const size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed));
    const size_t max_flip_or_edge = num_edges * 2;

    NTSparseGraph g(directed, num_nodes);

    std::random_device rd;  // Will provide a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    EdgeSampler creation_sampler(g, gen);
    for (size_t i = 0; i < num_edges; i++) {
        Edge e = creation_sampler.sample_non_edge();
        g.add_edge(e.first, e.second);
    }

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

    // TODO: Passing 0 for num_threads does not work!
    ThreadPoolScorer TPS(2, g, comb_util,
                         nt_result.node_orbits, nt_result.edge_orbits,
                         log2_p_plus, log2_p_minus,
                         log2_1_minus_p_plus, log2_1_minus_p_minus);

    std::cout<<"Finished Creating Thread Pool Scorer..."<<std::endl;
    TPS.terminate();
    std::cout<<"Terminated."<<std::endl;

    return 0;
}
