#include "coloring.h"
#include "edge.h"
// file_utils.h is not needed for this code, but is included to that all
//  references to this library can be included in a single header file
#include "file_utils.h"
// same with debugging.h
#include "debugging.h"
#include "nt_sparse_graph.h"

#include<unordered_map>
#include<vector>

#ifndef SCHENO__TRACES_H
#define SCHENO__TRACES_H

struct NautyTracesOptions {
    // Set to true to collect the orbits of the nodes and the number of distinct
    //  (node) orbits.
    bool get_node_orbits;

    // Set to true to collect the orbits of the edges and the number of distinct
    //  (edge) orbits.
    bool get_edge_orbits;

    // Set to true to get the canonical node order.
    bool get_canonical_node_order;
};

struct NautyTracesResults {
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

    // Maps nodes to their automorphism orbit ID. Used iff
    //  get_node_orbits is true.
    //
    // NOTE: The orbit ids are not canonical. To make them canonical, you must
    //  use the canonical node order to relabel these orbits.
    //
    // NOTE: The orbit ids are not forced to be 0 through num-orbits-minus-1.
    //  Rather, they can be anything. However, they will not overlap with edge
    //  orbit ids.
    Coloring<int> node_orbits;

    // Maps edges to their automorphism orbit ID. Used iff
    //  get_edge_orbits is true.
    //
    // NOTE: The orbit ids are not canonical. To make them canonical, you must
    //  use the canonical node order to relabel these orbits.
    //
    // NOTE: The orbit ids are not forced to be 0 through num-orbits-minus-1.
    //  Rather, they can be anything. However, they will not overlap with node
    //  orbit ids.
    Coloring<Edge, EdgeHash> edge_orbits;
};

// Even though g is not passed as a const, it is left un-modified.
NautyTracesResults nauty(NTSparseGraph& g, const NautyTracesOptions& o);

// Even though g is not passed as a const, it is left un-modified.
//  However, p might be modified.
NautyTracesResults nauty(NTSparseGraph& g, const NautyTracesOptions& o,
                         NTPartition& p);

// Even though g is not passed as a const, it is left un-modified.
NautyTracesResults traces(NTSparseGraph& g, const NautyTracesOptions& o);

// Even though g is not passed as a const, it is left un-modified.
//  However, p might be modified.
NautyTracesResults traces(NTSparseGraph& g, const NautyTracesOptions& o,
                          NTPartition& p);

// fake_iso() is a Weisfeiler-Lehman-based algorithm that is not guaranteed to
// produce correct results. However, it is faster and will work for most graphs.
//
// It is not part of the original nauty/traces code. Rather it was added by
//  Justus Hibshman.

// Even though g is not passed as a const, it is left un-modified.
NautyTracesResults fake_iso(NTSparseGraph& g, const NautyTracesOptions& o);

// Even though g is not passed as a const, it is left un-modified.
//  However, p might be modified.
NautyTracesResults fake_iso(NTSparseGraph& g, const NautyTracesOptions& o,
                            NTPartition& p);


////////////////// These Variables are for Internal Use Only ///////////////////

// This class is used to prevent unnecessary de- and re-allocation of space.
//
// It manages the space for where Nauty and Traces write out the canonical
//  sparsegraph and the node orbits.
class __NTRunSpace {
public:
    __NTRunSpace();
    ~__NTRunSpace();

    void set_size(size_t n, size_t nde);
    sparsegraph g;
    int *orbits;

protected:
    size_t actual_n;
    size_t actual_elen;
};

static thread_local __NTRunSpace __nt_run_space = __NTRunSpace();

#endif
