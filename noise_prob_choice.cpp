#include<cmath>

#include "noise_prob_choice.h"

#include "nauty_traces.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"

std::vector<long double> log2_noise_probs_empty_g(NTSparseGraph& g,
                                      const CombinatoricUtility& comb_util) {
    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g, o);

    bool directed = g.directed;
    size_t num_nodes = g.num_nodes();
    size_t num_edges = g.num_edges();

    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed)) +
            (num_nodes * size_t(g.num_loops() > 0));

    long double alpha;
    if (num_edges < (max_possible_edges / 2)) {
        alpha = std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(num_edges));
    } else {
        alpha = std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(max_possible_edges - num_edges));
    }
    long double log2_p_plus = std::log2l(alpha) - std::log2l(1 + alpha);
    long double log2_1_minus_p_plus = -std::log2l(1 + alpha);
    long double log2_p_minus = log2_p_plus;
    long double log2_1_minus_p_minus = log2_1_minus_p_plus;

    return {log2_p_plus, log2_1_minus_p_plus,
            log2_p_minus, log2_1_minus_p_minus};
}

std::vector<long double> log2_noise_probs_empty_g_full(NTSparseGraph& g,
                                      const CombinatoricUtility& comb_util) {
    NautyTracesOptions o;

    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g, o);

    bool directed = g.directed;
    size_t num_nodes = g.num_nodes();
    size_t num_edges = g.num_edges();

    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed));

    long double k1 =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(num_edges));

    long double k2 =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                        (long double)(max_possible_edges - num_edges));

    // Here's hoping the decimal places work OK!
    //  TODO: Double-check the floating point math for safety
    long double log2_p_plus = std::log2l(k1 - (k1 * k2)) -
                              std::log2l(1.0 - (k1 * k2));
    long double log2_1_minus_p_plus = std::log2l(1.0 - k1) -
                                      std::log2l(1.0 - (k1 * k2));

    long double log2_p_minus = std::log2l(k2 - (k1 * k2)) -
                               std::log2l(1.0 - (k1 * k2));
    long double log2_1_minus_p_minus = std::log2l(1.0 - k2) -
                                       std::log2l(1.0 - (k1 * k2));

    return {log2_p_plus, log2_1_minus_p_plus,
            log2_p_minus, log2_1_minus_p_minus};
}
