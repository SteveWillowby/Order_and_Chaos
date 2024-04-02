#include "nt_wrappers/nauty_traces.h"

#include<random>
#include<unordered_set>
#include<vector>

// Converts, Weights, and Samples edges as ints

#ifndef SCHENO__INT_EDGE_SAMPLER_H
#define SCHENO__INT_EDGE_SAMPLER_H

class IntEdgeConverterAndSampler {
public:
    IntEdgeConverterAndSampler(const Graph& g, bool weight_samples);

    // When this constructor is used, then only edges from `legal_edges` are
    //  randomly sampled by sample()
    //
    // IMPORTANT: simple_sample() will not respect the `legal_edges` constraint
    IntEdgeConverterAndSampler(const Graph& g,
                               const Graph& legal_edges, bool weight_samples);

    // dist should be std::uniform_real_distribution<long double>(0, 1.0)
    SCHENO__edge_int_type sample(std::mt19937& gen,
                std::uniform_real_distribution<long double>& dist) const;

    // dist should be
    //  std::uniform_int_distribution<SCHENO__edge_int_type>(0, (n * n) - 1)
    //
    // Quicker. Does not use weights. Does not respect `legal_edges`.
    SCHENO__edge_int_type simple_sample(std::mt19937& gen,
                std::uniform_int_distribution<SCHENO__edge_int_type>& dist) const;

    bool is_edge(SCHENO__edge_int_type e) const;
    Edge edge(SCHENO__edge_int_type e) const;
    SCHENO__edge_int_type edge(const Edge& e) const;

    // Maps an edge int to its weighted probability of being sampled
    const std::vector<long double>& get_heuristic_scores() const;

protected:
    const bool directed;
    const size_t n;
    const bool self_loops;
    std::unordered_set<SCHENO__edge_int_type> edges;

    // Maps an edge int to its weighted probability of being sampled
    std::vector<long double> heuristic_scores;

    std::vector<long double> cumulative_sums;

    // Used to account for the fact that due to floating point issues,
    //  the cumulative sum might be < 1.0
    long double total_sum;
};

#endif
