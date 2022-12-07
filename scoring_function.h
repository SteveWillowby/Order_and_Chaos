#include "edge.h"
#include "nt_sparse_graph.h"

#include<unordered_set>
#include<vector>

#ifndef SYM__SCORING_FUNCTION_H
#define SYM__SCORING_FUNCTION_H

// WARNING: THE IMPLEMENTATION CURRENTLY USES O(n^2) SPACE.
//  It does so for the sake of memoization.


// `g` will be modified but will be returned to its original state,
//      externally speaking (i.e. the internal details might change, but the
//      graph it represents will be the same).
//
// Requires that `log2_g_aut` = log2(|Aut(G)|)
//
// `node_orbit_coloring` should describe the automorphism orbits of `g`'s nodes.
// `edge_orbit_coloring` should describe the automorphism orbits of `g`'s edges.
// `editable_edge_orbit_coloring` should be identical to `edge_orbit_coloring`.
//      note that though `editable_edge_orbit_coloring` will be temporarily
//      modified, it will represent the same coloring when the function is done.
double score(double log2_g_aut,
             NTSparseGraph& g,
             const Coloring<int>& node_orbit_coloring,
             const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
             Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
             const std::unordered_set<Edge,EdgeHash>& edge_additions,
             const std::unordered_set<Edge,EdgeHash>& edge_removals);


/* TODO: Implement the NTSparseGraph features needed for this faster version.

// `g` will be modified but will be returned to its original state,
//      externally speaking (i.e. the internal details might change, but the
//      graph it represents will be the same).
// `g` should have the structure_coloring turned on.
//
// `g_edge_tracker` should be identitcal to `g` but should have the change
//      highlights coloring turned on.
//  For an extra speedup, `g_edge_tracker` can have its reference coloring set
//      to the automorphism orbits; this is not required however.
//
// Requires that `log2_g_aut` = log2(|Aut(G)|)
double score(double g_num_aut_base, int g_num_aut_exponent,
             NTSparseGraph& g,
             NTSparseGraph& g_edge_tracker,
             const std::unordered_set<Edge,EdgeHash>& edge_additions,
             const std::unordered_set<Edge,EdgeHash>& edge_removals);
*/


// This class is used to memoize log and factorial operations.
class __CombinatoricUtility {
public:
    __CombinatoricUtility();

    double log2(size_t x);
    double log2_factorial(size_t x);

    // Calling with b > a will cause an error.
    double log2_a_choose_b(size_t a, size_t b);

protected:
    std::vector<double> log2_factorials;
    std::vector<double> log2_s;
};

static thread_local __CombinatoricUtility
            __combinatoric_utility = __CombinatoricUtility();

#endif
