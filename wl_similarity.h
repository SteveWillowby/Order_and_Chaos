#include "graph.h"

#ifndef SYM__WL_SIMILARITY_H
#define SYM__WL_SIMILARITY_H

// Returns either pw_scores_1 or pw_scores_2, whichever contain the result
//
// `node` is a node for which node-centric similarity measures will be
//   calculated -- if -1 is passed, no node-centrality is used
//
// `cost_matrix` should be an array of length n*n
// `col_for_row` and `row_for_col` should be a ptrdiff_t array of length n
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
// If sum_result is true, then the negation of all the entries in the matrix
//  will be summed together and put in neg_sum.
double* wl_similarity_measure(bool sum_result, long double* neg_sum,
                              const Graph& g, const size_t node,
                              double* cost_matrix,
                              ptrdiff_t* col_for_row,
                              ptrdiff_t* row_for_col,
                              double* u, double* v,
                              void* workspace,
                              double* pw_scores_1,
                              double* pw_scores_2,
                              size_t* start_indices);

#endif
