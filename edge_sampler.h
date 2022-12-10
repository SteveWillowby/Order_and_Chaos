#include "edge.h"
#include "sparse_graph.h"

#include<cstdint>
#include<random>

#ifndef SYM__EDGE_SAMPLER_H
#define SYM__EDGE_SAMPLER_H

// Note: if 2^32 < n^2, increase this size.
typedef edge_int_type uint32_t;

class EdgeSampler {
public:
    EdgeSampler(const SparseGraph& g, std::mt19937& generator);

    Edge sample_edge();
    Edge sample_non_edge();
    void un_sample_edge(const Edge& e);
    void un_sample_non_edge(const Edge& e);

protected:
    const bool directed;
    const size_t n;
    const bool has_self_loops;

    std::vector<edge_int_type> edges_un_sampled;
    std::vector<edge_int_type> edges_sampled;
    std::vector<edge_int_type> non_edges_un_sampled;
    std::vector<edge_int_type> non_edges_sampled;

    Edge size_t_to_edge(edge_int_type e) const;
    edge_int_type edge_to_size_t(const Edge& e) const;
};

#endif
