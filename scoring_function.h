#include "nt_wrappers/nauty_traces.h"

#include<array>
#include<unordered_set>
#include<vector>

#ifndef SCHENO__SCORING_FUNCTION_H
#define SCHENO__SCORING_FUNCTION_H

// WARNING: THE IMPLEMENTATION CURRENTLY USES A LOT OF SPACE
//  It does so for the sake of memoizing log and factorial computations.
//  (See below for more details)
//
// *IN ORDER TO USE THIS CODE* on n-node graphs, you must have initialized
//  a CombinatoricUtility object with max_e and max_f set such that:
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
//
//  NOTE: The CombinatoricUtility methods are thread-safe EXCEPT for
//      update_max_access()

class CombinatoricUtility;

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
//
// Against what might be expected, `log2_p_plus` is NOT the log-probability that
//  an edge is added to `g`. RATHER, `log2_p_plus` is the log of the probability
//  that any given noise edge is added to the hypothesis graph (i.e. removed
//  from `g`).
// The same (in reverse) is true for `log2_p_minus`.
//
// `log2_1_minus_p_plus` = log2(1 - p_plus)
long double score(NTSparseGraph& g, const CombinatoricUtility& comb_util,
                  const Coloring<int>& node_orbit_coloring,
                  const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
                  Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals,
                  const long double log2_p_plus, const long double log2_p_minus,
                  const long double log2_1_minus_p_plus,
                  const long double log2_1_minus_p_minus,
                  const size_t max_change,
                  bool full_iso);

// This function is similar to the above, except it gives a breakdown
//  of the score into components. Here 'H' means hypothesis graph, 'N' means
//  the noise set, and 'H U N' means the edgeset union of the two.
//
// a[0] -- the full score
// a[1] -- log(symmetry of connected components)             (positive number)
//              i.e. log(symmetry) - a[2]
// a[2] -- log(symmetry of singletons)                       (positive number)
//              i.e. log((# singletons in H)!)
// a[3] -- log(AO of noise ignoring singleton swapping)      (positive number)
//              i.e. log(AO of noise) - a[4]
// a[4] -- log(AO of noise due to connecting to singletons)  (positive number)
//              i.e. log((# sing. in H)! / (# sing. in H U N)!)
// a[5] -- log(probability of noise size)                    (negative number)
//
// a[0] should equal the sum of a[1] through a[5]
//
// Further, this version does all its edits to the graph BEFORE turning the
//  graph into an NTSparseGraph
//
// If `full_iso` is set to true, the program will use traces like normal.
// If it is set to false, then the program will use the WL-based "fake_iso".
std::array<long double, 6>
  score_breakdown(SparseGraph& g, const CombinatoricUtility& comb_util,
                  const Coloring<int>& node_orbit_coloring,
                  const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
                  Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals,
                  const long double log2_p_plus, const long double log2_p_minus,
                  const long double log2_1_minus_p_plus,
                  const long double log2_1_minus_p_minus,
                  const size_t max_change,
                  bool full_iso);
                  
std::pair<long double, long double>
            score(NTSparseGraph& g, const CombinatoricUtility& comb_util,
                  const Coloring<int>& node_orbit_coloring,
                  const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
                  Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals,
                  const long double log2_p_plus, const long double log2_p_minus,
                  const long double log2_1_minus_p_plus,
                  const long double log2_1_minus_p_minus,
                  const size_t max_change,
                  const double* precomputed_wl_diff,
                  double* cost_matrix,
                  ptrdiff_t* col_for_row, ptrdiff_t* row_for_col,
                  double* u, double* v,
                  void* workspace, double* pw_scores_1, double* pw_scores_2,
                  size_t* start_indices,
                  bool full_iso);

// Outdated
long double score(NTSparseGraph& g, const CombinatoricUtility& comb_util,
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
// Against what might be expected, `log2_p_plus` is NOT the log-probability that
//  an edge is added to `g`. RATHER, `log2_p_plus` is the log of the probability
//  that an edge in the hypothesis graph is removed.
// The same (in reverse) is true for `log2_p_minus`.
//
// `log2_1_minus_p_plus` = log2(1 - p_plus)
long double score(NTSparseGraph& g, const CombinatoricUtility& comb_util,
                  NTSparseGraph& g_edge_tracker,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals,
                  const double log2_p_plus, const double log2_p_minus,
                  const double log2_p_plus, const double log2_p_minus,
                  const double log2_1_minus_p_plus,
                  const double log2_1_minus_p_minus);
*/


// This class is used to memoize log and factorial operations.
class CombinatoricUtility {
public:
    // Documentation is at the top of this file.
    CombinatoricUtility(size_t max_e, size_t max_f);

    long double log2(size_t x) const;
    long double log2_factorial(size_t x) const;

    // Documentation is at the top of this file.
    //
    // NOT thread-safe at present.
    void update_max_access(size_t max_e, size_t max_f);

    // Calling with b > a will cause an error.
    long double log2_a_choose_b(size_t a, size_t b) const;

protected:
    // edge_flip_start_1 = 0
    size_t edge_flip_end_1;  // exclusive index
    size_t edge_flip_start_2;
    size_t edge_flip_end_2;  // exclusive index
    std::vector<long double> log2_factorials[2];
    std::vector<long double> log2_s[2];
};

#endif
