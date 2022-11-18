#include "edge.h"
#include "nt_sparse_graph.h"

#include<unordered_map>
#include<vector>

#ifndef SYM__TRACES_H
#define SYM__TRACES_H

struct SYMTracesOptions {
    // Set to true to collect the orbits of the nodes and the number of distinct
    //  (node) orbits.
    bool get_node_orbits;
    // Set to true to collect the orbits of the edges and the number of distinct
    //  (edge) orbits.
    bool get_edge_orbits;
    // Set to true to get the canonical node order.
    bool get_canonical_node_order;
};

struct SYMTracesResults {
    // Will be non-zero if Traces encountered an error
    int error_status;

    // Number of automorphism orbits of the nodes.
    size_t num_node_orbits;
    // Number of automorphism orbits of the edges.
    size_t num_edge_orbits;

    // The number of automorphisms is roughly
    //  num_aut_base * (10 ^ num_aut_exponent)
    double num_aut_base;
    int num_aut_exponent;

    // Lists nodes in a canonical order. Used iff get_canonical_node_order is
    //  true.
    std::vector<int> canonical_node_order;

    // Maps nodes to their canonical automorphism orbit ID. Used iff
    //  get_node_orbits is true.
    //
    // NOTE: The orbit ids are not forced to be 0 through num-orbits-minus-1.
    //  Rather, they can be anything. However, they will not overlap with edge
    //  orbit ids.
    std::vector<int> node_orbits;

    // Maps edges to their canonical automorphism orbit ID. Used iff
    //  get_edge_orbits is true.
    //
    // NOTE: The orbit ids are not forced to be 0 through num-orbits-minus-1.
    //  Rather, they can be anything. However, they will not overlap with node
    //  orbit ids.
    std::unordered_map<Edge, int, EdgeHash> edge_orbits;
};

// Even though g is not passed as a const, it is left un-modified.
SYMTracesResults traces(NTSparseGraph& g, const SYMTracesOptions& o);

#endif
