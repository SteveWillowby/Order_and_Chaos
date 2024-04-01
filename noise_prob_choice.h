#include "scoring_function.h"

#include "nt_wrappers/nauty_traces.h"

#ifndef SCHENO__NOISE_PROB_CHOICE_H
#define SCHENO__NOISE_PROB_CHOICE_H
#define default_log2_noise_probs log2_noise_probs_fancy_equality


// Returns, in order:
//  log2_p_plus -- log2(prob an edge was added to original)
//  log2_1_minus_p_plus
//  log2_p_minus -- log2(prob an edge was removed from original)
//  log2_1_minus_p_minus
//
// This function chooses a p_plus = p_minus such that the following holds:
//
//  P(g | ER w/ p)        P(Data = 0 | S) * P(Noise = g | 0, ER w/ p)
// -----------------   = ---------------------------------------------
//  P(g | S)              P(Data = g | S) * P(Noise = 0 | g, ER w/ p)
//
// Where "ER w/ p" means Erdos-Renyi noise generation with probability p, and
//  S means the structure distribution (the inverse of ER w/ .5).
//
// IMPORTANT: Assumes that the combinatoric utility can return log2(n!)
//      (i.e. that max_flip_or_edge is at least n + 1)
std::vector<long double> log2_noise_probs_fancy_equality(NTSparseGraph& g,
                                      const CombinatoricUtility& comb_util);
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
