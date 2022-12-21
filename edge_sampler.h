#include "edge.h"
#include "sparse_graph.h"

#include<cstdint>
#include<random>

#ifndef SYM__EDGE_SAMPLER_H
#define SYM__EDGE_SAMPLER_H

// Note: if 2^32 < n^2, increase this size and adjust max edge label
typedef uint32_t SYM__edge_int_type;
#define SYM__MAX_EDGE_LABEL 0xFFFFFFFF

class EdgeSampler {
public:
    // O(n^2) time and space
    EdgeSampler(const Graph& g, std::mt19937& gen);

    // Randomly selects an edge which has not yet been sampled.
    // O(1)
    Edge sample_edge();
    // Randomly selects a non-edge which has not yet been sampled.
    // O(1)
    Edge sample_non_edge();
    // Puts an edge into the list of edges that may be sampled.
    // O(1)
    Edge un_sample_edge();
    // Puts a non-edge into the list of non-edges that may be sampled.
    // O(1)
    Edge un_sample_non_edge();

    // Replaces a sampled edge with a different sampled edge.
    //  returns <Newly Sampled, Newly Un-Sampled>
    // O(1)
    std::pair<Edge, Edge> swap_edge_samples();

    // O(1)
    std::pair<Edge, Edge> swap_non_edge_samples();

    // Un-does the last (and only the last) (un-)sample operation.
    // O(1)
    void undo();

protected:
    const bool directed;
    const size_t n;
    const bool has_self_loops;
    std::mt19937& generator;
    std::uniform_real_distribution<double> dist;

    // 0 -- none
    // 1 -- sampled edge
    // 2 -- sampled non-edge
    // 3 -- un-sampled edge
    // 4 -- un-sampled non-edge
    // 5 -- swap edge
    // 6 -- swap non-edge
    uint8_t last_op;

    std::vector<SYM__edge_int_type> edges_un_sampled;
    std::vector<SYM__edge_int_type> edges_sampled;
    std::vector<SYM__edge_int_type> non_edges_un_sampled;
    std::vector<SYM__edge_int_type> non_edges_sampled;

    Edge int_to_edge(SYM__edge_int_type e) const;
    // SYM__edge_int_type edge_to_int(const Edge& e) const;
};

#endif
