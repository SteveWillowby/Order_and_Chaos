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

    internal_n = n;
    undirected_m = 0;

    out_degrees = std::vector<int>(n, 0);

    // Begin with ORIG_EDGE_SPACE_PER_NODE edge slots per node.
    node_to_startpoint = std::vector<int>();
    node_to_endpoint = std::vector<int>();
    startpoint_to_node = std::unordered_map<int, int>();
    for (size_t i = 0; i < n; i++) {
        node_to_startpoint.push_back(ORIG_EDGE_SPACE_PER_NODE * i);
        node_to_endpoint.push_back(ORIG_EDGE_SPACE_PER_NODE * (i + 1));
        startpoint_to_node.insert(
            std::pair<int, int>(ORIG_EDGE_SPACE_PER_NODE * i, i));
    }

    out_neighbors_vec = std::vector<int>(n * ORIG_EDGE_SPACE_PER_NODE, 0);

    chunks = std::vector<std::pair<size_t, size_t>>();
}

NTSparseGraph::NTSparseGraph(const Graph &g) : SparseGraph(g) {

}


int NTSparseGraph::add_node() {
    int new_node = SparseGraph::add_node();

    return new_node;
}

int NTSparseGraph::delete_node(const int a) {
    int replacement_node = SparseGraph::delete_node(a);

    return replacement_node;
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
