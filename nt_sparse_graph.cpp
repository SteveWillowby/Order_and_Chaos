#include "edge.h"
#include "sparse_graph.h"
#include "nt_sparse_graph.h"

#include<stdexcept>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<utility>
#include<vector>


NTSparseGraph::NTSparseGraph(const bool directed) : SparseGraph(directed) {

}

NTSparseGraph::NTSparseGraph(const bool directed, size_t n)
                    : SparseGraph(directed, n) {
}

NTSparseGraph::NTSparseGraph(const Graph &g) : SparseGraph(g) {

}


int NTSparseGraph::add_node() {
    return SparseGraph::add_node();
}

int NTSparseGraph::delete_node(const int a) {
    return SparseGraph::delete_node(a);
}


bool NTSparseGraph::add_edge(const int a, const int b) {
    return SparseGraph::add_edge(a, b);
}

bool NTSparseGraph::delete_edge(const int a, const int b) {
    return SparseGraph::delete_edge(a, b);
}

void NTSparseGraph::flip_edge(const int a, const int b) {
    SparseGraph::flip_edge(a, b);
}
