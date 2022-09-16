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
    nodes_with_extra_space = std::vector<int>();

    // Begin with MIN_EDGE_SPACE_PER_NODE edge slots per node.
    node_to_startpoint = std::vector<int>();
    node_to_endpoint = std::vector<int>();
    startpoint_to_node = std::unordered_map<int, int>();
    for (size_t i = 0; i < n; i++) {
        node_to_startpoint.push_back(MIN_EDGE_SPACE_PER_NODE * i);
        node_to_endpoint.push_back(MIN_EDGE_SPACE_PER_NODE * (i + 1));
        startpoint_to_node.insert(
            std::pair<int, int>(MIN_EDGE_SPACE_PER_NODE * i, i));
    }

    out_neighbors_vec = std::vector<int>(n * MIN_EDGE_SPACE_PER_NODE, 0);
}

NTSparseGraph::NTSparseGraph(const Graph &g) : SparseGraph(g) {

}


int NTSparseGraph::add_node() {
    int new_node = SparseGraph::add_node();

    if (undirected_m == 0) {
        out_degrees.push_back(0);

        if (nodes_with_extra_space.size() > 0) {
            int node_with_space =
                    nodes_with_extra_space[nodes_with_extra_space.size() - 1];
            int new_start_point = (node_to_startpoint[node_with_space] + 
                                   out_degrees[node_with_space]) / 2;
        } else {
            
        }

        return new_node;
    }

    // Add the new real node in slot n.
    // Move the first edge node to the back of the list.
    out_degrees.push_back(out_degrees[n]);
    out_degrees[n] = 0;

    return new_node;
}

int NTSparseGraph::delete_node(const int a) {
    int replacement_node = SparseGraph::delete_node(a);

    return replacement_node;
}


bool NTSparseGraph::add_edge(const int a, const int b) {
    bool added = SparseGraph::add_edge(a, b);

    return added;
}

bool NTSparseGraph::delete_edge(const int a, const int b) {
    bool deleted = SparseGraph::delete_edge(a, b);

    return deleted;
}

void NTSparseGraph::flip_edge(const int a, const int b) {
    // The implementation makes a call to add_edge or to delete_edge.
    SparseGraph::flip_edge(a, b);
}
