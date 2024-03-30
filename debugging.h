#include<iostream>
#include<string>
#include<vector>

#include "nt_sparse_graph.h"

std::string vec_as_string(const std::vector<int> &v);


#ifdef SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
bool consistency_check(const NTSparseGraph &g);

std::vector<int> cleaned_out_N_vec(const NTSparseGraph &g);

std::string nt_graph_as_string(const NTSparseGraph &g);

void print_graph(const NTSparseGraph &g);

void rand_test(float add_node_prob, float delete_node_prob,
               float add_edge_prob, float delete_edge_prob,
               const bool directed, size_t iterations,
               size_t reconstruct_frequency);
#endif
