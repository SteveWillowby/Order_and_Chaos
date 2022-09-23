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
    num_edge_nodes = 0;

    out_degrees = std::vector<int>(n, 0);
    nodes_with_extra_space = std::vector<int>();

    // Begin with MIN_EDGE_SPACE_PER_NODE edge slots per node.
    node_to_startpoint = std::vector<int>();
    node_to_endpoint = std::vector<int>();
    startpoint_to_node = std::unordered_map<int, int>();
    for (size_t i = 0; i < n; i++) {
        node_to_startpoint.push_back(MIN_EDGE_SPACE_PER_NODE * i);
        node_to_endpoint.push_back(MIN_EDGE_SPACE_PER_NODE * (i + 1));
        startpoint_to_node[MIN_EDGE_SPACE_PER_NODE * i] = i;
    }

    out_neighbors_vec = std::vector<int>(n * MIN_EDGE_SPACE_PER_NODE, 0);
}

NTSparseGraph::NTSparseGraph(const Graph &g) : SparseGraph(g) {

}


// O(log n) amortized
int NTSparseGraph::add_node() {
    int new_node = SparseGraph::add_node();

    if (num_edge_nodes == 0) {
        // No edges in the graph.
        out_degrees.push_back(0);
    } else {
        int new_edge_node = allocate_edge_node();
        relabel_edge_node(n - 1, new_edge_node);
        out_degrees[n - 1] = 0;
    }

    auto space = extra_space_and_node.lower_bound(-MIN_EDGE_SPACE_PER_NODE);
    size_t edge_node_start = out_neighbors_vec.size() - (num_edge_nodes * 2);
    if (space == extra_space_and_node.end()) {
        // Create space for the new node's neighbors.
        for (size_t i = 0; i < MIN_EDGE_SPACE_PER_NODE; i += 2) {
            if (num_edge_nodes == 0) {
                out_neighbors_vec.push_back(0);
                out_neighbors_vec.push_back(0);
            } else {
                // Slide edge node info to the back.

                edge_node_start += 2; // Now the endpoint of the first edge node
                int edge_node = endpoint_to_node.find(edge_node_start)->second;
                endpoint_to_node[out_neighbors_vec.size()] = edge_node;

                out_neighbors_vec.push_back(
                    out_neighbors_vec[node_to_startpoint[edge_node]]);
                out_neighbors_vec.push_back(
                    out_neighbors_vec[node_to_startpoint[edge_node] + 1]);
                if (directed) {
                    // One of this edge node's neighbors is itself an edge node.
                    //  Update *that* edge node's `places` info.
                    neighbor_edge_node =
                        out_neighbors_vec[out_neighbors_vec.size() - 1];
                    const std::pair<int, int> &nen_locs =
                        edge_node_to_places(neighbor_edge_node);
                    edge_node_to_places[neighbors_edge_node] =
                          std::pair<int,int>(nen_locs.first,
                                             node_to_startpoint[edge_node] + 1);
                }

                node_to_endpoint[edge_node] = out_neighbors_vec.size();
                node_to_startpoint[edge_node] = out_neighbors_vec.size() - 2;
            }
        }
    } else {
        // Collect info.
        int capacity_for_new_node = space->first;
        int old_node = space->second;
        // Split existing node space to accommodate new node.
        extra_space_and_node.erase(space);

        int old_node_new_capacity = 0;
        if (old_node == -1) {
            // Update start- and end-points.
            node_to_endpoint[new_node] = capacity_for_new_node;
            node_to_startpoint[new_node] = 0;

        } else {
            old_node_new_capacity = (node_to_endpoint[old_node] -
                                      node_to_startpoint[old_node]) -
                                        capacity_for_new_node;
            // Update start- and end-points.
            node_to_endpoint[new_node] = node_to_endpoint[old_node];
            node_to_startpoint[new_node] = node_to_endpoint[old_node] - 
                                            capacity_for_new_node;
            node_to_endpoint[old_node] = node_to_startpoint[new_node];
        }


        // Update the extra capacity information.
        int capacity;

        // New Node
        //  The new node currently has zero edges.
        if (capacity_for_new_node >= 2 * MIN_EDGE_SPACE_PER_NODE) {
            capacity = capacity_for_new_node / 2;
            extra_space_and_node.insert(
                std::pair<int, int>(capacity, new_node));
            has_extra_space.push_back(true);
        } else {
            has_extra_space.push_back(false);
        }

        // Old Node
        if (old_node != -1) {
            if (old_node_new_capacity >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                    old_node_new_capacity >= 4 * out_degrees[old_node]) {
                capacity = old_node_new_capacity / 2;
                extra_space_and_node.insert(
                    std::pair<int, int>(capacity, old_node));
                // has_extra_space[old_node] is already true
            } else {
                has_extra_space[old_node] = false;
            }
        }
        // Done updating extra capacity information.
    }

    return new_node;
}

// O()
int NTSparseGraph::delete_node(const int a) {
    int replacement_node = SparseGraph::delete_node(a);

    return replacement_node;
}

// O(1) amortized
bool NTSparseGraph::add_edge(const int a, const int b) {
    bool added = SparseGraph::add_edge(a, b);

    return added;
}

// O(1) amortized
bool NTSparseGraph::delete_edge(const int a, const int b) {
    bool deleted = SparseGraph::delete_edge(a, b);

    return deleted;
}

void NTSparseGraph::flip_edge(const int a, const int b) {
    // The implementation makes a call to add_edge or to delete_edge.
    // Since these are virtual functions, they should call the NTSparseGraph
    //  version.
    SparseGraph::flip_edge(a, b);
}


// Returns the new internal label for the new edge node.
int NTSparseGraph::allocate_edge_node() {

    int new_label = out_degrees.size();

    // An edge node always has degree 2.
    out_degrees.push_back(2);
    out_neighbors_vec.push_back(0);
    out_neighbors_vec.push_back(0);
    node_to_startpoint.push_back(out_neighbors.size() - 2);
    node_to_endpoint.push_back(out_neighbors.size());
    edge_node_to_places.push_back(std::pair<int,int>(0, 0));

    return new_label;
}

// Given space allocated for an edge node at slot b, move a's edge info to that
//  slot, relabeling things as necessary.
void NTSparseGraph::relabel_edge_node(const int a, const int b) {

    const std::pair<int, int> locations = edge_node_to_places.find(a)->second;

    // Update out_neighbors_vec
    out_neighbors_vec[locations.first] = b;
    out_neighbors_vec[locations.second] = b;

    out_neighbors_vec[node_to_startpoint[b]] =
        out_neighbors_vec[node_to_startpoint[a]];
    out_neighbors_vec[node_to_startpoint[b] + 1] =
        out_neighbors_vec[node_to_startpoint[a] + 1];

    // Update edge_node_to_places
    edge_node_to_places.erase(a);
    edge_node_to_places[b] = locations;

    // Update edge_to_edge_node and edge_node_to_edge:
    auto e_itr = edge_node_to_edge.find(a);
    //  Overwrite the old label.
    edge_to_edge_node[e_itr->second] = b;
    edge_node_to_edge[b] = e_itr->second;
    edge_node_to_edge.erase(a);

    // Do not need to update out_degrees since all edge nodes have degree 2.
}
