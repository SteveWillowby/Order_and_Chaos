#include "edge.h"
#include "graph.h"

#include<random>
#include<unordered_set>
#include<vector>

// Converts, Weights, and Samples edges as ints

#ifndef SYM__INT_EDGE_SAMPLER_H
#define SYM__INT_EDGE_SAMPLER_H

class IntEdgeConverterAndSampler {
public:
    IntEdgeConverterAndSampler(const Graph& g);

    // dist should be std::uniform_real_distribution<long double>(0, 1.0)
    SYM__edge_int_type weighted_sample(std::mt19937& gen,
                std::uniform_real_distribution<long double>& dist) const;

    // dist should be
    //  std::uniform_int_distribution<SYM__edge_int_type>(0, n * n)
    SYM__edge_int_type unweighted_sample(std::mt19937& gen,
                std::uniform_int_distribution<SYM__edge_int_type>& dist) const;

    bool is_edge(SYM__edge_int_type e) const;
    Edge edge(SYM__edge_int_type e) const;

    // Maps an edge int to its weighted probability of being sampled
    const std::vector<long double>& get_heuristic_scores() const;

protected:
    const bool directed;
    const size_t n;
    const bool self_loops;
    std::unordered_set<SYM__edge_int_type> edges;

    // Maps an edge int to its weighted probability of being sampled
    std::vector<long double> heuristic_scores;

    std::vector<long double> cumulative_sums;

    // Used to account for the fact that due to floating point issues,
    //  the cumulative sum might be < 1.0
    long double total_sum;
};

#endif
