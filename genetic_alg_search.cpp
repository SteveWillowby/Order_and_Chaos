#include "edge.h"
#include "genetic_alg_search.h"
#include "graph.h"
#include "thread_pool_scorer.h"

#include<rand>
#include<unordered_set>
#include<vector>

class GeneCollection;

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
    // Initialize Basics

    NTSparseGraph g_nt(g);
    bool directed = g_nt.directed();
    size_t num_nodes = g_nt.num_nodes();
    size_t num_edges = g_nt.num_edges();

    // Genetic Alg Variables

    const size_t DEPTH = 2;  // How many layers of gene heirarchy are there?
    size_t POP_SIZE = num_nodes * 10;
    // The top TOP_POP_SIZE performers are guaranteed to be kept.
    //  The remaining (POP_SIZE - TOP_POP_SIZE) are selected randomly from
    //      the broad population.
    size_t TOP_POP_SIZE = POP_SIZE / 2;

    size_t MATE_SIZE_A = (POP_SIZE * POP_SIZE) / 1000;
    size_t MATE_SIZE_B = POP_SIZE * 10;
    size_t MATE_SIZE = (MATE_SIZE_A < MATE_SIZE_B ? MATE_SIZE_B : MATE_SIZE_A);
    size_t MUTATE_SIZE = POP_SIZE / 10;

    // Initialization of Scoring Thread Pool
    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed));
    size_t max_flip_or_edge = num_edges * 5;

    CombinatoricUtility comb_util(max_possible_edges, max_flip_or_edge);

    NautyTracesOptions nt_options;
    nt_options.get_node_orbits = true;
    nt_options.get_edge_orbits = true;
    nt_options.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g_nt, nt_options);

    // Calculate the noise probability at which the full graph is equally
    //  likely to be the noise as it is to be the structure.
    long double alpha_plus =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(num_edges));
    long double log2_p_plus = std::log2l(alpha_plus) -
                              std::log2l(1.0 + alpha_plus);
    long double log2_1_minus_p_plus = -std::log2l(1.0 + alpha_plus);

    long double alpha_minus =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                        (long double)(max_possible_edges - num_edges));
    long double log2_p_minus = std::log2l(alpha_minus) -
                               std::log2l(1.0 + alpha_minus);
    long double log2_1_minus_p_minus = -std::log2l(1.0 + alpha_minus);

    ThreadPoolScorer TPS(nt, g_nt, comb_util,
                         nt_result.node_orbits, nt_result.edge_orbits,
                         log2_p_plus, log2_p_minus,
                         log2_1_minus_p_plus, log2_1_minus_p_minus);


    // Initialization of Gene Population

    auto top_results =
       std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>();

    // TODO: At lower levels, make edges represented as ints, and represent
    //  genes as ints as well. That way we only need to convert to actual edges
    //  at the top level.

    // TODO: Keep track of when a gene or meta-gene first appears so that we can
    //  give them unique numeric IDs. That way, we can convert the population to
    //  an actual edge set only at the time of evaluation. (space savings)

    // TODO: Make a recursive merge function that handles conflicts.
}

class GeneCollection {
public:
    GeneCollection(bool directed, size_t num_nodes, size_t depth);

protected:

};
