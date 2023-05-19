#include "nt_sparse_graph.h"
#include "scoring_function.h"

#ifndef SYM__NOISE_PROB_CHOICE_H
#define SYM__NOISE_PROB_CHOICE_H

// Returns, in order:
//  log2_p_plus -- log2(prob an edge was added to original)
//  log2_1_minus_p_plus
//  log2_p_minus -- log2(prob an edge was removed from original)
//  log2_1_minus_p_minus
#define default_log2_noise_probs log2_noise_probs_empty_g

// Returns, in order:
//  log2_p_plus -- log2(prob an edge was added to original)
//  log2_1_minus_p_plus
//  log2_p_minus -- log2(prob an edge was removed from original)
//  log2_1_minus_p_minus
//
// This particular function chooses p_plus = p_minus such that the graph
//  is equally likely to be the data as it is to be noise added to the empty
//  graph.
//
// Even though g is not passed as const, it is left un-modified.
std::vector<long double> log2_noise_probs_empty_g(NTSparseGraph& g,
                                          const CombinatoricUtility& comb_util);

// Returns, in order:
//  log2_p_plus -- log2(prob an edge was added to original)
//  log2_1_minus_p_plus
//  log2_p_minus -- log2(prob an edge was removed from original)
//  log2_1_minus_p_minus
//
// This particular function chooses p_plus and p_minus such that the graph
//  is equally likely to be the data as it is to be noise added to the empty
//  graph AND is equally likely to be noise taken from the full graph.
//
// Even though g is not passed as const, it is left un-modified.
std::vector<long double> log2_noise_probs_empty_g_full(NTSparseGraph& g,
                                          const CombinatoricUtility& comb_util);

#endif
