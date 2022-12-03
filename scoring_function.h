#include "nt_sparse_graph.h"

// `g` will be modified but will be returned to its original state,
//      externally speaking (i.e. the internal details might change, but the
//      graph it represents will be the same).
// `g` should have the structure_coloring turned on.
//
// `g_edge_tracker` should be identitcal to `g` but should have the change
//      highlights coloring turned on.
//
// `changes` is the set of edges and/or non-edges that are to be flipped to get
//      the candidate graph.
double score(NTSparseGraph& g,
             NTSparseGraph& g_edge_tracker,
             const std::unordered_set<Edge,EdgeHash>& changes);
