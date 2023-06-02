#include "graph.h"

#ifndef SYM__SCORING_HEURISTIC_H
#define SYM__SCORING_HEURISTIC_H

long double wl_symmetry_measure(const Graph& g, double* cost_matrix,
                                ptrdiff_t* col_for_row, ptrdiff_t* row_for_col,
                                double* u, double* v, void* workspace,
                                double* pw_scores_1, double* pw_scores_2,
                                size_t* start_indices);

#endif
