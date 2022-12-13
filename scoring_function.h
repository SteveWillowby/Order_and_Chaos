#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"

#include<unordered_set>
#include<vector>

#ifndef SYM__SCORING_FUNCTION_H
#define SYM__SCORING_FUNCTION_H

// WARNING: THE IMPLEMENTATION CURRENTLY USES A LOT OF SPACE
//  It does so for the sake of memoizing log and factorial computations.
//  (See below for more details)
//
// *IN ORDER TO USE THIS CODE* on n-node graphs, each thread must have called
//  __combinatoric_utility.set_max_access(max_e, max_f) where:
//    * max_e >= the max number of edges your graph COULD ever have
//                  (i.e. n choose 2 if undirected, n * (n - 1) if directed,
//                      [+n if you allow self-loops])
//    * max_f >= min(the max number of edges your graph WILL have PLUS
//                       the number of edges you will delete from it,
//                   (the max number of edges you COULD have MINUS
//                        the number of edges your graph WILL have) PLUS
//                            the number of edges you will delete from it)
//
//       It's probably a safe bet to make:
//                        max_f >= min(3 * edges, max_e - (3 * (max_e - edges)))
//
//  NOTE: This call will take O(max_e) time, and the data structure will use
//   O(max_f) space.


// `g` will be modified but will be returned to its original state,
//      externally speaking (i.e. the internal details might change, but the
//      graph it represents will be the same).
//
// `node_orbit_coloring` should describe the automorphism orbits of `g`'s nodes.
// `edge_orbit_coloring` should describe the automorphism orbits of `g`'s edges.
// `editable_edge_orbit_coloring` should be identical to `edge_orbit_coloring`.
//      note that though `editable_edge_orbit_coloring` will be temporarily
//      modified, it will represent the same coloring when the function is done.
//
// `edge_additions` are the edges added to `g`. Thus they are the edges deleted
//      from the hypothesis graph.
// `edge_removals` are the edges deleted from `g`. Thus they are the edges
//      added to the hypothesis graph.
long double score(NTSparseGraph& g,
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
long double score(NTSparseGraph& g,
                  NTSparseGraph& g_edge_tracker,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals);
*/


// This class is used to memoize log and factorial operations.
class __CombinatoricUtility {
public:
    __CombinatoricUtility();

    long double log2(size_t x);
    long double log2_factorial(size_t x);

    // Documentation is at the top of this file.
    void set_max_access(size_t max_e, size_t max_f);

    // Calling with b > a will cause an error.
    long double log2_a_choose_b(size_t a, size_t b);

protected:
    // edge_flip_start_1 = 0
    size_t edge_flip_end_1;  // exclusive index
    size_t edge_flip_start_2;
    size_t edge_flip_end_2;  // exclusive index
    std::vector<long double>[2] log2_factorials;
    std::vector<long double>[2] log2_s;
};

static thread_local __CombinatoricUtility
            __combinatoric_utility = __CombinatoricUtility();

#endif
