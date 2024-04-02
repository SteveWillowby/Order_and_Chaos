#include "nt_wrappers/nauty_traces.h"

#ifndef SCHENO__WL_MEASURES_H
#define SCHENO__WL_MEASURES_H

// Returns either pw_scores_1 or pw_scores_2, whichever contain the result
//
// `node` is a node for which node-centric similarity measures will be
//   calculated -- if -1 is passed, no node-centrality is used
//
// `cost_matrix` should be an array of length n*n
// `col_for_row` and `row_for_col` should be ptrdiff_t arrays of length n
// `u` and `v` should be arrays of length n
// `workspace` should be initialized as malloc(assign2DCBufferSize(n, n))
// `pw_scores_1` and `pw_scores_2` should be arrays of length (n choose 2) + n
// `start_indices` should be an array of n entries where si_0 = 0 and
//      si_i = (si_{i-1} + (n - i) + 1)
//   e.g. if n = 10, si_0 = 0, si_1 = 10, si_2 = 19, si_3 = 27, si_4 = 34, ...
//
// Total difference between nodes a and b will be stored at
//  result[start_indices[min(a, b)] + abs(a - b)]
//
// `node_coloring` and `edge_coloring` can be NULL, in which case they are not
//      used.
//
//  If `node_coloring` is not NULL, then any node with a different color from
//      another automatically has a difference of 1.
//
//  If `edge_coloring` is not NULL, then a pair of edges with different colors
//      from one another cannot line up in the bipartite matching.
//    For example:
//      node A connects to node C with edge color red
//      node B connects to node D with edge color blue
//          When computing the difference between A and B, the connection cost
//            between C and D will be 1, regardless of C's difference from D.
double* wl_node_differences(const Graph& g,
                            const Coloring<int>* node_coloring,
                            const Coloring<Edge, EdgeHash>* edge_coloring,
                            double* cost_matrix,
                            ptrdiff_t* col_for_row,
                            ptrdiff_t* row_for_col,
                            double* u, double* v,
                            void* workspace,
                            double* pw_scores_1,
                            double* pw_scores_2,
                            const size_t* start_indices);

// Same arguments as above*. This function returns a measure loosly representing
//  the log of the number of automorphisms of the graph, but since it is a fuzzy
//  measure, it will be larger than the actual log of automorphisms. That's part
//  of the point.
//
// *The one extra argument is `precomputed_basics`.
//      `precomputed_basics` can be NULL, or it can be a matrix like pw_scores_1
//          If it's a matrix, then it is used as the start difference costs.
//          This can be useful when the edge coloring or node coloring is a
//            refinement for a previous run of wl_node_differences()
long double wl_symmetry_measure(const Graph& g,
                                const Coloring<int>* node_coloring,
                                const Coloring<Edge, EdgeHash>* edge_coloring,
                                const double* precomputed_basics,
                                double* cost_matrix,
                                ptrdiff_t* col_for_row,
                                ptrdiff_t* row_for_col,
                                double* u, double* v,
                                void* workspace,
                                double* pw_scores_1,
                                double* pw_scores_2,
                                const size_t* start_indices);

// Given the result of wl_node_differences() as `pw_scores`,
//  tells you loosely speaking what `node`'s fuzzy orbit size is.
// This number will always be >= `node`'s actual orbit size.
long double wl_orbit_size(const size_t node, const size_t num_nodes,
                          const double* pw_scores,
                          const size_t* start_indices);

#endif
