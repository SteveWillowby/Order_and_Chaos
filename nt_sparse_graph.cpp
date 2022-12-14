#include "augmented_multimap.h"
#include "coloring.h"
#include "edge.h"
#include "sparse_graph.h"
#include "nt_partition.h"
#include "nt_sparse_graph.h"

#include<stdexcept>
#include<unordered_map>
#include<unordered_set>
#include<utility>
#include<vector>

NTSparseGraph::NTSparseGraph(const bool directed) : NTSparseGraph(directed,1) {}

NTSparseGraph::NTSparseGraph(const bool directed, size_t n)
                    : SparseGraph(directed, n) {

    internal_n = n;
    num_edge_nodes = 0;

    out_degrees = std::vector<int>(n, 0);
    has_self_loop = std::vector<bool>(n, false);

    // Begin with MIN_EDGE_SPACE_PER_NODE edge slots per node.
    node_to_startpoint = std::vector<size_t>();
    node_to_endpoint = std::vector<size_t>();
    endpoint_to_node = std::unordered_map<size_t, int>();
    for (size_t i = 0; i < n; i++) {
        node_to_startpoint.push_back(MIN_EDGE_SPACE_PER_NODE * i);
        node_to_endpoint.push_back(MIN_EDGE_SPACE_PER_NODE * (i + 1));
        endpoint_to_node[MIN_EDGE_SPACE_PER_NODE * (i + 1)] = i;
    }

    out_neighbors_vec = std::vector<int>(n * MIN_EDGE_SPACE_PER_NODE, 0);

    extra_space_and_node = AugmentedMultimap<size_t, int>();

    edge_to_edge_node = std::unordered_map<Edge, int, EdgeHash>();
    edge_node_to_edge = std::unordered_map<int, Edge>();

    edge_node_to_places = std::unordered_map<int, std::pair<size_t,size_t>>();
}

// O(n + m)
//
// The 1-node effects of the call to SparseGraph(d, 1) are overwritten by the
//  call to copy_assignment(g) which in turn calls SparseGraph::operator=(g).
NTSparseGraph::NTSparseGraph(const Graph &g) : SparseGraph(g.directed, 1) {
    copy_assignment(g);
}

const sparsegraph NTSparseGraph::as_nauty_traces_graph() {
    sparsegraph g;
    SG_INIT(g);
    g.nv = internal_n;
    if (directed) {
        // For every 2 edge nodes there are 3 undirected edges,
        //  which means that for every 2 edge nodes there are 6 directed edges.
        g.nde = num_edge_nodes * 3;
    } else {
        // For every edge node there are 4 directed edges.
        g.nde = num_edge_nodes * 4;
    }

    g.v = &(node_to_startpoint[0]);
    g.vlen = internal_n;

    g.d = &(out_degrees[0]);
    g.dlen = internal_n;

    g.e = &(out_neighbors_vec[0]);
    g.elen = out_neighbors_vec.size();

    g.w = NULL;
    g.wlen = 0;

    return g;
}


NTSparseGraph& NTSparseGraph::operator=(const NTSparseGraph& g) {
    return copy_assignment(g);
}
NTSparseGraph& NTSparseGraph::operator=(const SparseGraph& g) {
    return copy_assignment(g);
}
NTSparseGraph& NTSparseGraph::operator=(const Graph& g) {
    return copy_assignment(g);
}

NTSparseGraph& NTSparseGraph::copy_assignment(const Graph& g) {
    if (g.directed != directed) {
        throw std::logic_error(
                std::string("Error! Cannot perform copy assignment when one") +
                " graph is directed and the other is not.");
    }
    
    if (this == &g) {
        return *this;
    }
    
    SparseGraph::operator=(g);
    
    // SparseGraph's = operator does not make any calls to functions like
    //  add_node or add_edge. Thus we can call it and then fill out the NT info
    //  separately.
    if (directed) {
        num_edge_nodes = 0;
        for (int i = 0; i < int(n); i++) {
            num_edge_nodes += _neighbors[i].size() -
                size_t(_neighbors[i].find(i) != _neighbors[i].end());
        }
        internal_n = n + num_edge_nodes;

        if (internal_n > NAUTY_TRACES_MAXN) {
            throw std::logic_error(
               std::string("Error! Too many nodes for Nauty/Traces. Directed") +
               std::string(" graphs must have N + 2 * (undirected) M <= ") +
               std::to_string(NAUTY_TRACES_MAXN));
        }

        // Add basic node info.
        out_degrees = std::vector<int>(internal_n, 2);
        node_to_startpoint = std::vector<size_t>(internal_n, 0);
        node_to_endpoint = std::vector<size_t>(internal_n, 0);
        endpoint_to_node = std::unordered_map<size_t, int>();
        has_self_loop = std::vector<bool>(n, false);

        size_t total_space_needed = 0;

        for (int i = 0; i < int(n); i++) {
            out_degrees[i] = _neighbors[i].size();
            if (_neighbors[i].find(i) != _neighbors[i].end()) {
                out_degrees[i]--;
                has_self_loop[i] = true;
            }
            node_to_startpoint[i] = total_space_needed;
            total_space_needed +=
                (size_t(out_degrees[i]) < MIN_EDGE_SPACE_PER_NODE ?
                    MIN_EDGE_SPACE_PER_NODE : out_degrees[i]);
            node_to_endpoint[i] = total_space_needed;
            endpoint_to_node[total_space_needed] = i;
        }

        size_t edge_node_start = total_space_needed;
        total_space_needed += num_edge_nodes * 2;
        out_neighbors_vec = std::vector<int>(total_space_needed, 0);

        edge_to_edge_node = std::unordered_map<Edge, int, EdgeHash>();
        edge_node_to_edge = std::unordered_map<int, Edge>();
        edge_node_to_places =
                std::unordered_map<int, std::pair<size_t,size_t>>();
        
        extra_space_and_node = AugmentedMultimap<size_t, int>();

        size_t loc = edge_node_start;
        for (size_t en = n; en < internal_n; en++) {
            node_to_startpoint[en] = loc;
            loc += 2;
            node_to_endpoint[en] = loc;
            endpoint_to_node[loc] = en;
        }

        // Add edge and edge node info.

        std::vector<size_t> slots_used = std::vector<size_t>(n, 0);
        int next_edge_node = n;
        size_t loc_1, loc_2, en_loc_1, en_loc_2;
        int en_1, en_2;

        for (int i = 0; i < int(n); i++) {
            for (auto nbr = _neighbors[i].begin();
                            nbr != _neighbors[i].end(); nbr++) {
                if (*nbr <= i) {
                    // Only count an edge once, and ignore self-loops.
                    continue;
                }

                en_1 = next_edge_node++;
                en_2 = next_edge_node++;
                loc_1 = node_to_startpoint[i] + slots_used[i];
                out_neighbors_vec[loc_1] = en_1;
                slots_used[i]++;
                loc_2 = node_to_startpoint[*nbr] + slots_used[*nbr];
                out_neighbors_vec[loc_2] = en_2;
                slots_used[*nbr]++;
                
                edge_to_edge_node[EDGE(i, *nbr, true)] = en_1;
                edge_to_edge_node[EDGE(*nbr, i, true)] = en_2;
                edge_node_to_edge[en_1] = EDGE(i, *nbr, true);
                edge_node_to_edge[en_2] = EDGE(*nbr, i, true);
                
                en_loc_1 = edge_node_start + (en_1 - n) * 2;
                en_loc_2 = en_loc_1 + 2;
                
                out_neighbors_vec[en_loc_1] = i;
                out_neighbors_vec[en_loc_1 + 1] = en_2;
                out_neighbors_vec[en_loc_2] = *nbr;
                out_neighbors_vec[en_loc_2 + 1] = en_1;
                
                edge_node_to_places[en_1] =
                    std::pair<size_t, size_t>(loc_1, en_loc_2 + 1);
                edge_node_to_places[en_2] =
                    std::pair<size_t, size_t>(loc_2, en_loc_1 + 1);
            }
        }
    } else {
        // Undirected
        num_edge_nodes = 0;
        for (int i = 0; i < int(n); i++) {
            num_edge_nodes += _neighbors[i].size() -
                size_t(_neighbors[i].find(i) != _neighbors[i].end());
        }
        num_edge_nodes /= 2;
        internal_n = n + num_edge_nodes;

        if (internal_n > NAUTY_TRACES_MAXN) {
            throw std::logic_error(
               std::string("Error! Too many nodes for Nauty/Traces. ") +
               std::string("Undirected graphs must have N + M <= ") +
               std::to_string(NAUTY_TRACES_MAXN));
        }

        // Add basic node info.
        out_degrees = std::vector<int>(internal_n, 2);
        node_to_startpoint = std::vector<size_t>(internal_n, 0);
        node_to_endpoint = std::vector<size_t>(internal_n, 0);
        endpoint_to_node = std::unordered_map<size_t, int>();
        has_self_loop = std::vector<bool>(n, false);

        size_t total_space_needed = 0;

        for (int i = 0; i < int(n); i++) {
            out_degrees[i] = _neighbors[i].size();
            if (_neighbors[i].find(i) != _neighbors[i].end()) {
                out_degrees[i]--;
                has_self_loop[i] = true;
            }
            node_to_startpoint[i] = total_space_needed;
            total_space_needed +=
                (size_t(out_degrees[i]) < MIN_EDGE_SPACE_PER_NODE ?
                    MIN_EDGE_SPACE_PER_NODE : out_degrees[i]);
            node_to_endpoint[i] = total_space_needed;
            endpoint_to_node[total_space_needed] = i;
        }

        size_t edge_node_start = total_space_needed;
        total_space_needed += num_edge_nodes * 2;
        out_neighbors_vec = std::vector<int>(total_space_needed, 0);

        edge_to_edge_node = std::unordered_map<Edge, int, EdgeHash>();
        edge_node_to_edge = std::unordered_map<int, Edge>();
        edge_node_to_places =
                std::unordered_map<int, std::pair<size_t,size_t>>();
        
        extra_space_and_node = AugmentedMultimap<size_t, int>();

        size_t loc = edge_node_start;
        for (size_t en = n; en < internal_n; en++) {
            node_to_startpoint[en] = loc;
            loc += 2;
            node_to_endpoint[en] = loc;
            endpoint_to_node[loc] = en;
        }

        // Add edge and edge node info.

        std::vector<size_t> slots_used = std::vector<size_t>(n, 0);
        int next_edge_node = n;
        size_t loc_1, loc_2, en_loc;
        int en; // edge node

        for (int i = 0; i < int(n); i++) {
            for (auto nbr = _neighbors[i].begin();
                            nbr != _neighbors[i].end(); nbr++) {
                if (*nbr <= i) {
                    // Only count an edge once, and ignore self-loops.
                    continue;
                }

                en = next_edge_node++;
                loc_1 = node_to_startpoint[i] + slots_used[i];
                out_neighbors_vec[loc_1] = en;
                slots_used[i]++;
                loc_2 = node_to_startpoint[*nbr] + slots_used[*nbr];
                out_neighbors_vec[loc_2] = en;
                slots_used[*nbr]++;
                
                edge_to_edge_node[EDGE(i, *nbr, false)] = en;
                edge_node_to_edge[en] = EDGE(i, *nbr, false);
                
                en_loc = edge_node_start + (en - n) * 2;
                
                // Remember that i > *nbr.
                out_neighbors_vec[en_loc] = i;
                out_neighbors_vec[en_loc + 1] = *nbr;
                
                edge_node_to_places[en] =
                    std::pair<size_t, size_t>(loc_1, loc_2);
            }
        }
    }
    
    return *this;
}


// O(log n) amortized
int NTSparseGraph::add_node() {
    if (internal_n >= NAUTY_TRACES_MAXN) {
        if (directed) {
            throw std::logic_error(
               std::string("Error! Too many nodes for Nauty/Traces. Directed") +
               std::string(" graphs must have N + 2 * (undirected) M <= ") +
               std::to_string(NAUTY_TRACES_MAXN));
        } else {
            throw std::logic_error(
               std::string("Error! Too many nodes for Nauty/Traces. Undirec") +
               std::string("ted graphs must have N + M <= ") +
               std::to_string(NAUTY_TRACES_MAXN));
        }
    }

    int new_node = SparseGraph::add_node();

    has_self_loop.push_back(false);

    if (num_edge_nodes > 0) {
        // Take the edge node with label n-1 and turn it into node internal_n-1
        //  so that the regular nodes can be labeled 0 through n-1.

        out_degrees.push_back(out_degrees[n-1]);
        node_to_startpoint.push_back(node_to_startpoint[n-1]);
        node_to_endpoint.push_back(node_to_endpoint[n-1]);

        // Update the following for the edge node:
        //  * endpoint_to_node
        //  * edge_to_edge_node
        //  * edge_node_to_edge
        //  * edge_node_to_places
        //  * out_neighbors_vec
        relabel_edge_node(n - 1, out_degrees.size() - 1);

        out_degrees[n - 1] = 0;
        node_to_startpoint[n - 1] = -1;  // placeholder values - uninitialized
        node_to_endpoint[n - 1] = -1;  // placeholder values - uninitialized
    } else {
        out_degrees.push_back(0);
        node_to_startpoint.push_back(-1);  // placeholder values - uninitialized
        node_to_endpoint.push_back(-1);  // placeholder values - uninitialized
    }

    move_node_to_more_space(new_node);

    internal_n++;

    return new_node;
}

// O(N(a) + N(n - 1))
//
//  The runtime stems from the fact that a's edges have to be deleted and
//   node n-1's edges have to be relabeled to now refer to "a".
int NTSparseGraph::delete_node(const int a) {

    // This info is about to be erased. Save it for now.
    std::vector<int> neighbors = std::vector<int>(_neighbors[a].begin(),
                                                  _neighbors[a].end());

    int replacement_node = SparseGraph::delete_node(a);

    // Delete all edge nodes for a's edges.
    for (auto itr = neighbors.begin(); itr < neighbors.end(); itr++) {
        delete_edge_node_or_nodes(a, *itr);
    }

    // Give away node a's extra space.
    size_t a_startpoint = node_to_startpoint[a];
    size_t a_endpoint = node_to_endpoint[a];
    size_t slot_size = a_endpoint - a_startpoint;
    // All of a's neighbors got deleted in the loop above. Thus we know that
    //  slot_size >= 4 * out_degree, so we only need the MIN_... condition.
    if (slot_size >= 2 * MIN_EDGE_SPACE_PER_NODE) {
        // Node a was listed as having extra capacity. Delete that.
        extra_space_and_node.erase(slot_size / 2, a);
    }

    auto left_node_itr = endpoint_to_node.find(a_startpoint);
    if (left_node_itr == endpoint_to_node.end()) {
        // Left node is '-1'
        if (a_startpoint > 0) {
            extra_space_and_node.erase(a_startpoint, -1);
        }
        extra_space_and_node.insert(a_endpoint, -1);
        endpoint_to_node.erase(a_endpoint);
    } else {
        int left_node = left_node_itr->second;

        // Hand over the space.
        size_t ln_start = node_to_startpoint[left_node];
        size_t ln_end_old = node_to_endpoint[left_node];
        node_to_endpoint[left_node] = a_endpoint;
        endpoint_to_node[a_endpoint] = left_node;
        endpoint_to_node.erase(a_startpoint);

        // Perform extra capacity updates.
        size_t ln_end = a_endpoint;
        size_t old_size = ln_end_old - ln_start;
        size_t new_size = ln_end - ln_start;
        size_t degree = out_degrees[left_node];

        if (old_size >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                old_size >= 4 * degree) {
            extra_space_and_node.erase(old_size / 2, left_node);
        }
        if (new_size >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                new_size >= 4 * degree) {
            extra_space_and_node.insert(new_size / 2, left_node);
        }
    }

    if (replacement_node != a) {
        // We need to relabel replacement_node to be called 'a'.

        // Copy basic info.
        has_self_loop[a] = has_self_loop[replacement_node];
        out_degrees[a] = out_degrees[replacement_node];
        node_to_startpoint[a] = node_to_startpoint[replacement_node];
        node_to_endpoint[a] = node_to_endpoint[replacement_node];
        endpoint_to_node[node_to_endpoint[a]] = a;

        // If replacement node had extra space, relabel that as belonging to a.
        size_t capacity = node_to_endpoint[a] - node_to_startpoint[a];
        if (capacity >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                capacity >= 4 * size_t(out_degrees[a])) {
            extra_space_and_node.insert(capacity / 2, a);
            extra_space_and_node.erase(capacity / 2, replacement_node);
        }

        // Handle the edges.
        if (directed) {
            for (auto itr = _neighbors[a].begin();
                      itr != _neighbors[a].end(); itr++) {

                if (*itr == a) {
                    // A self-loop. Has no edge nodes.
                    //  NOTE: This self-loop was originally a self-loop from
                    //  replacement_node to itself, but that got relabeled by
                    //  SparseGraph::delete_node(a).
                    continue;
                }

                auto en_A_itr =
                     edge_to_edge_node.find(EDGE(replacement_node, *itr, true));
                auto en_B_itr =
                     edge_to_edge_node.find(EDGE(*itr, replacement_node, true));
                int en_A = en_A_itr->second;
                int en_B = en_B_itr->second;
                edge_to_edge_node.erase(en_A_itr);
                edge_to_edge_node.erase(en_B_itr);
                edge_to_edge_node[EDGE(a, *itr, true)] = en_A;
                edge_to_edge_node[EDGE(*itr, a, true)] = en_B;

                edge_node_to_edge[en_A] = EDGE(a, *itr, true);
                edge_node_to_edge[en_B] = EDGE(*itr, a, true);

                out_neighbors_vec[node_to_startpoint[en_A]] = a;
            }
        } else {
            // In an undirected graph, the relabeling of a node might mean that
            //  for some edge_node_to_places pairs, the order of the two places
            //  must be swapped if the order of the places' node ids have been
            //  swapped.
            for (auto itr = _neighbors[a].begin();
                      itr != _neighbors[a].end(); itr++) {

                if (*itr == a) {
                    // A self-loop. Has no edge nodes.
                    //  NOTE: This self-loop was originally a self-loop from
                    //  replacement_node to itself, but that got relabeled by
                    //  SparseGraph::delete_node(a).
                    continue;
                }

                int min_node = (a < *itr ? a : *itr);
                int max_node = (a < *itr ? *itr : a);

                auto en_itr =
                     edge_to_edge_node.find(EDGE(replacement_node, *itr,false));
                int edge_node = en_itr->second;
                edge_to_edge_node.erase(en_itr);
                edge_to_edge_node[EDGE(a, *itr, false)] = edge_node;

                edge_node_to_edge[edge_node] = EDGE(a, *itr, false);

                out_neighbors_vec[node_to_startpoint[edge_node]] = min_node;
                out_neighbors_vec[node_to_startpoint[edge_node] + 1] = max_node;

                if ((a < *itr) != (replacement_node < *itr)) {
                    // The order of the node IDs swapped. Swap them in
                    //  edge_node_to_places.
                    const auto &locs = edge_node_to_places[edge_node];
                    edge_node_to_places[edge_node] =
                        std::pair<size_t, size_t>(locs.second, locs.first);
                }
            }
        }
    }

    // Lastly, now that we have relabeled replacement_node to 'a', we should
    //  relabel the largest edge node (if one exists) to 'replacement_node'.
        
    internal_n--;

    //  Relabel the largest edge node to replacement_node.
    if (num_edge_nodes > 0) {
        node_to_startpoint[replacement_node] = node_to_startpoint[internal_n];
        node_to_endpoint[replacement_node] = node_to_endpoint[internal_n];
        out_degrees[replacement_node] = 2;
        relabel_edge_node(internal_n, replacement_node);
    }

    // Delete surplus vector info.
    out_degrees.pop_back();
    node_to_startpoint.pop_back();
    node_to_endpoint.pop_back();
    has_self_loop.pop_back();

    return replacement_node;
}

// O(1) amortized
bool NTSparseGraph::add_edge(const int a, const int b) {
    bool added = SparseGraph::add_edge(a, b);
    if (!added) {
        return false;
    }

    if (a == b) {
        has_self_loop[a] = true;
        return true;
    }

    if (directed && _out_neighbors[b].find(a) != _out_neighbors[b].end()) {
        // The internal nodes already exist.
        return true;
    }

    if (directed && internal_n >= NAUTY_TRACES_MAXN - 1) {
        throw std::logic_error(
            std::string("Error! Too many nodes for Nauty/Traces. Directed") +
            std::string(" graphs must have N + 2 * (undirected) M <= ") +
            std::to_string(NAUTY_TRACES_MAXN));
    } else if (!directed && internal_n >= NAUTY_TRACES_MAXN) {
        throw std::logic_error(
            std::string("Error! Too many nodes for Nauty/Traces. Undirected") +
            std::string(" graphs must have N + M <= ") +
            std::to_string(NAUTY_TRACES_MAXN));
    }

    if (size_t(out_degrees[a]) ==
            (node_to_endpoint[a] - node_to_startpoint[a])) {
        // Move node's edge info to a new place.
        move_node_to_more_space(a);
    }
    if (size_t(out_degrees[b]) ==
            (node_to_endpoint[b] - node_to_startpoint[b])) {
        // Move node's edge info to a new place.
        move_node_to_more_space(b);
    }

    if (directed) {
        internal_n += 2;
        num_edge_nodes += 2;

        int edge_node_A = allocate_edge_node();
        int edge_node_B = allocate_edge_node();

        out_neighbors_vec[node_to_startpoint[a] +
                          size_t(out_degrees[a])] = edge_node_A;
        out_neighbors_vec[node_to_startpoint[b] +
                          size_t(out_degrees[b])] = edge_node_B;
        out_neighbors_vec[node_to_startpoint[edge_node_A]] = a;
        out_neighbors_vec[node_to_startpoint[edge_node_A] + 1] = edge_node_B;
        out_neighbors_vec[node_to_startpoint[edge_node_B]] = b;
        out_neighbors_vec[node_to_startpoint[edge_node_B] + 1] = edge_node_A;

        edge_node_to_places[edge_node_A] =
            std::pair<size_t, size_t>(node_to_startpoint[a] +
                                        size_t(out_degrees[a]),
                                      node_to_startpoint[edge_node_B] + 1);
        edge_node_to_places[edge_node_B] =
            std::pair<size_t, size_t>(node_to_startpoint[b] +
                                        size_t(out_degrees[b]),
                                      node_to_startpoint[edge_node_A] + 1);

        edge_to_edge_node[EDGE(a, b, directed)] = edge_node_A;
        edge_to_edge_node[EDGE(b, a, directed)] = edge_node_B;
        edge_node_to_edge[edge_node_A] = EDGE(a, b, directed);
        edge_node_to_edge[edge_node_B] = EDGE(b, a, directed);

        endpoint_to_node[node_to_endpoint[edge_node_A]] = edge_node_A;
        endpoint_to_node[node_to_endpoint[edge_node_B]] = edge_node_B;

    } else {
        internal_n++;
        num_edge_nodes++;

        int edge_node = allocate_edge_node();
        out_neighbors_vec[node_to_startpoint[a] +
                          size_t(out_degrees[a])] = edge_node;
        out_neighbors_vec[node_to_startpoint[b] +
                          size_t(out_degrees[b])] = edge_node;

        // Put the smaller node ID first.
        out_neighbors_vec[node_to_startpoint[edge_node]] = (a < b ? a : b);
        out_neighbors_vec[node_to_startpoint[edge_node] + 1] = (a < b ? b : a);

        edge_to_edge_node[EDGE(a, b, directed)] = edge_node;
        edge_node_to_edge[edge_node] = EDGE(a, b, directed);

        int min_node = (a < b ? a : b);
        int max_node = (a < b ? b : a);
        edge_node_to_places[edge_node] = std::pair<size_t, size_t>(
                  node_to_startpoint[min_node] + size_t(out_degrees[min_node]),
                  node_to_startpoint[max_node] + size_t(out_degrees[max_node]));

        endpoint_to_node[node_to_endpoint[edge_node]] = edge_node;
    }

    out_degrees[a]++;
    out_degrees[b]++;

    // Update extra_space_and_node
    size_t capacity_A = node_to_endpoint[a] - node_to_startpoint[a];
    size_t capacity_B = node_to_endpoint[b] - node_to_startpoint[b];
    size_t extra_capacity;
    if (size_t(out_degrees[a]) * 4 > capacity_A &&
            size_t(out_degrees[a] - 1) * 4 <= capacity_A &&
                size_t(capacity_A) >= MIN_EDGE_SPACE_PER_NODE * 2) {
        // Used to have space but no longer does.
        extra_capacity = capacity_A / 2;
        extra_space_and_node.erase(extra_capacity, a);
    }
    if (size_t(out_degrees[b]) * 4 > capacity_B &&
            size_t(out_degrees[b] - 1) * 4 <= capacity_B &&
                size_t(capacity_B) >= MIN_EDGE_SPACE_PER_NODE * 2) {
        // Used to have space but no longer does.
        extra_capacity = capacity_B / 2;
        extra_space_and_node.erase(extra_capacity, b);
    }

    return true;
}

// O(1) amortized
bool NTSparseGraph::delete_edge(const int a, const int b) {
    bool deleted = SparseGraph::delete_edge(a, b);

    if (!deleted) {
        return false;
    }

    if (directed &&
            _out_neighbors[b].find(a) != _out_neighbors[b].end()) {
        // The edge was deleted _but_ we need the same edge nodes for the
        //  reverse edge.
        return true;
    }

    delete_edge_node_or_nodes(a, b);

    return true;
}

void NTSparseGraph::delete_edge_node_or_nodes(const int a, const int b) {
    if (a == b) {
        has_self_loop[a] = false;
        return;
    }

    if (directed) {
        // Get edge node labels from edge_to_edge_node, then delete this edge
        //  from edge_to_edge_node.
        auto edge_node_itr = edge_to_edge_node.find(EDGE(a, b, directed));
        int edge_node_A = edge_node_itr->second;
        edge_to_edge_node.erase(edge_node_itr);
        edge_node_itr = edge_to_edge_node.find(EDGE(b, a, directed));
        int edge_node_B = edge_node_itr->second;
        edge_to_edge_node.erase(edge_node_itr);

        edge_node_to_edge.erase(edge_node_A);
        edge_node_to_edge.erase(edge_node_B);

        // Update out_neighbors_vec
        //  First do the regular nodes area
        const std::pair<size_t, size_t> &locations_A =
                            edge_node_to_places[edge_node_A];
        const std::pair<size_t, size_t> &locations_B =
                            edge_node_to_places[edge_node_B];
        //      Removes the refs to the edge nodes from the regular nodes' part
        //      of out_neighbors_vec.
        // NOTE: In the process, decrements a's and b's out_degrees
        remove_edge_node_ref(a, locations_A.first);
        remove_edge_node_ref(b, locations_B.first);
        edge_node_to_places.erase(edge_node_A);
        edge_node_to_places.erase(edge_node_B);

        //  Second, update the edge nodes area.
        //    Move the last edge node(s)' data into the slots of the deleted
        //      edge node(s).
        //
        //  Make sure to move the larger edge node first so as to avoid the case
        //  where by moving one edge node you relabel the other one (bec. the
        //  other was the largest).
        if (edge_node_A > edge_node_B) {
            slide_back_edge_node_to_slot(edge_node_A);
            slide_back_edge_node_to_slot(edge_node_B);
        } else {
            slide_back_edge_node_to_slot(edge_node_B);
            slide_back_edge_node_to_slot(edge_node_A);
        }

        internal_n -= 2;
        num_edge_nodes -= 2;
    } else {
        // Get edge node label from edge_to_edge_node, then delete this edge
        //  from edge_to_edge_node.
        auto edge_node_itr = edge_to_edge_node.find(EDGE(a, b, directed));
        int edge_node = edge_node_itr->second;
        edge_to_edge_node.erase(edge_node_itr);

        edge_node_to_edge.erase(edge_node);

        // Update out_neighbors_vec
        //  First do the regular nodes area
        const auto &locations = edge_node_to_places.find(edge_node);
        //      Removes the refs to the edge node from the regular nodes' part
        //      of out_neighbors_vec.
        // NOTE: In the process, decrements a's and b's out_degrees
        remove_edge_node_ref((a < b ? a : b), locations->second.first);
        remove_edge_node_ref((a < b ? b : a), locations->second.second);
        edge_node_to_places.erase(locations);

        //  Second, update the edge nodes area.
        //    Move the last edge node(s)' data into the slots of the deleted
        //      edge node(s).
        slide_back_edge_node_to_slot(edge_node);

        internal_n--;
        num_edge_nodes--;
    }


    // Update extra_space_and_node
    size_t capacity_A = node_to_endpoint[a] - node_to_startpoint[a];
    size_t capacity_B = node_to_endpoint[b] - node_to_startpoint[b];
    size_t extra_capacity;
    if (size_t(out_degrees[a]) * 4 <= capacity_A &&
            size_t(out_degrees[a] + 1) * 4 > capacity_A &&
                capacity_A >= MIN_EDGE_SPACE_PER_NODE * 2) {
        // Just acquired enough extra space.
        extra_capacity = capacity_A / 2;
        extra_space_and_node.insert(extra_capacity, a);
    }
    if (size_t(out_degrees[b]) * 4 <= capacity_B &&
            size_t(out_degrees[b] + 1) * 4 > capacity_B &&
                capacity_B >= MIN_EDGE_SPACE_PER_NODE * 2) {
        // Just acquired enough extra space.
        extra_capacity = capacity_B / 2;
        extra_space_and_node.insert(extra_capacity, b);
    }
}

void NTSparseGraph::flip_edge(const int a, const int b) {
    // The implementation makes a call to add_edge or to delete_edge.
    // Since these are virtual functions, they should call the NTSparseGraph
    //  version.
    SparseGraph::flip_edge(a, b);
}

int NTSparseGraph::edge_node(const int a, const int b) const {
    if (a == b) {
        throw std::logic_error(std::string("Error in call to edge_node()! ")
                               + "Self-loops do not have edge nodes.");
    }
    const auto& itr = edge_to_edge_node.find(EDGE(a, b, directed));
    if (itr == edge_to_edge_node.end()) {
        throw std::range_error(std::string("Error in call to edge_node()! ")
                               + "Edge (" + std::to_string(a) + ", "
                               + std::to_string(b) + ") does not exist.");
    }
    return itr->second;
}

// Returns the new internal label for the new edge node.
int NTSparseGraph::allocate_edge_node() {

    int new_label = out_degrees.size();

    // An edge node always has degree 2.
    out_degrees.push_back(2);
    out_neighbors_vec.push_back(0);
    out_neighbors_vec.push_back(0);
    node_to_startpoint.push_back(out_neighbors_vec.size() - 2);
    node_to_endpoint.push_back(out_neighbors_vec.size());
    edge_node_to_places[out_degrees.size() - 1] = std::pair<size_t,size_t>(0,0);

    return new_label;
}

// Change the edge node's label. REQUIRES that node_to_startpoint and
//  node_to_endpoint have already been updated for node b to point to node a.
void NTSparseGraph::relabel_edge_node(const int a, const int b) {

    const auto &locations_itr = edge_node_to_places.find(a);
    if (locations_itr == edge_node_to_places.end()) {
        // We are also in the process of deleting node b, so there's not much
        //  we should do in relabeling node a. 
        endpoint_to_node[node_to_endpoint[b]] = b;
        return;
    }
    const std::pair<size_t, size_t> locations = locations_itr->second;

    // Update out_neighbors_vec
    out_neighbors_vec[locations.first] = b;
    out_neighbors_vec[locations.second] = b;

    // Update edge_node_to_places
    edge_node_to_places.erase(a);
    edge_node_to_places[b] = locations;

    // Update edge_to_edge_node and edge_node_to_edge:
    auto e_itr = edge_node_to_edge.find(a);
    //  Overwrite the old label.
    edge_to_edge_node[e_itr->second] = b;
    edge_node_to_edge[b] = e_itr->second;
    edge_node_to_edge.erase(a);

    endpoint_to_node[node_to_endpoint[b]] = b;
}

// Used when moving a regular node's list of edge nodes.
//  Does not change the allocation or ID of the edge node.
void NTSparseGraph::move_edge_node_reference(const size_t init_loc,
                                             const size_t target_loc) {
    int edge_node = out_neighbors_vec[init_loc];
    out_neighbors_vec[target_loc] = edge_node;

    const std::pair<size_t, size_t> locations = edge_node_to_places[edge_node];

    if (init_loc == locations.first) {
        edge_node_to_places[edge_node] =
            std::pair<size_t, size_t>(target_loc, locations.second);
    } else {
        // init_loc == locations.second
        edge_node_to_places[edge_node] =
            std::pair<size_t, size_t>(locations.first, target_loc);
    }
}

void NTSparseGraph::slide_first_edge_node_to_back() {
    size_t edge_node_start = out_neighbors_vec.size() - (num_edge_nodes * 2);

    edge_node_start += 2; // Now the endpoint of the first edge node
    int edge_node = endpoint_to_node[edge_node_start];
    endpoint_to_node.erase(edge_node_start);

    out_neighbors_vec.push_back(
                    out_neighbors_vec[node_to_startpoint[edge_node]]);
    out_neighbors_vec.push_back(
                    out_neighbors_vec[node_to_startpoint[edge_node] + 1]);

    node_to_endpoint[edge_node] = out_neighbors_vec.size();
    node_to_startpoint[edge_node] = out_neighbors_vec.size() - 2;
    endpoint_to_node[out_neighbors_vec.size()] = edge_node;

    if (directed) {
        // One of this edge node's neighbors is itself an edge node.
        //  Update *that* edge node's `places` info.
        int neighbor_edge_node =
                        out_neighbors_vec[out_neighbors_vec.size() - 1];
        const std::pair<size_t, size_t> &nen_locs =
                        edge_node_to_places[neighbor_edge_node];
        edge_node_to_places[neighbor_edge_node] =
                    std::pair<size_t,size_t>(nen_locs.first,
                                             node_to_startpoint[edge_node] + 1);
    }
}

void NTSparseGraph::slide_back_edge_node_to_slot(int edge_node_of_slot) {
    size_t new_startpoint = node_to_startpoint[edge_node_of_slot];

    size_t old_endpoint = out_neighbors_vec.size();
    size_t old_startpoint = old_endpoint - 2;

    int moving_node = endpoint_to_node[old_endpoint];

    int largest_node = out_degrees.size() - 1;
    size_t largest_node_startpoint = node_to_startpoint[largest_node];
    size_t largest_node_endpoint = largest_node_startpoint + 2;


    out_neighbors_vec[new_startpoint] = out_neighbors_vec[old_startpoint];
    out_neighbors_vec[new_startpoint +1] = out_neighbors_vec[old_startpoint +1];

    endpoint_to_node[new_startpoint + 2] = moving_node;
    node_to_endpoint[moving_node] = new_startpoint + 2;
    node_to_startpoint[moving_node] = new_startpoint;

    out_neighbors_vec.pop_back();
    out_neighbors_vec.pop_back();
    out_degrees.pop_back();
    node_to_startpoint.pop_back();
    node_to_endpoint.pop_back();
    endpoint_to_node.erase(old_endpoint);

    int neighbor_edge_node = out_neighbors_vec[new_startpoint + 1];
    if (directed && neighbor_edge_node != edge_node_of_slot &&
                    moving_node != edge_node_of_slot) {
        // One of this edge node's neighbors is itself an edge node.
        // Further, we are not removing that edge node.
        //  Update *that* edge node's `places` info.
        const std::pair<size_t, size_t> &nen_locs =
                        edge_node_to_places[neighbor_edge_node];
        edge_node_to_places[neighbor_edge_node] =
                    std::pair<size_t,size_t>(nen_locs.first,
                                             new_startpoint + 1);
    }

    // Only relabel if edge_node_of_slot is less than the new largest node.
    if (edge_node_of_slot < largest_node) {
        if (moving_node != largest_node) {
            node_to_startpoint[edge_node_of_slot] = largest_node_startpoint;
            node_to_endpoint[edge_node_of_slot] = largest_node_endpoint;
        }

        relabel_edge_node(largest_node, edge_node_of_slot);
    }

    edge_node_to_places.erase(largest_node);

    // edge_node_to_edge and edge_to_edge_node were updated by
    //  relabel_edge_node()
}

void NTSparseGraph::move_node_to_more_space(const int a) {
    // Determine required amount of space
    size_t required_capacity = (node_to_endpoint[a] - node_to_startpoint[a]) *2;
    if (required_capacity < MIN_EDGE_SPACE_PER_NODE) {
        required_capacity = MIN_EDGE_SPACE_PER_NODE;
    }

    size_t old_startpoint = node_to_startpoint[a];
    size_t old_endpoint = node_to_endpoint[a];

    // Check if such space is available.
    const auto &space = extra_space_and_node.lower_bound(required_capacity);
    size_t edge_node_start = out_neighbors_vec.size() - (num_edge_nodes * 2);


    if (space.is_none()) {
        // Create space for the new node's neighbors.
        for (size_t i = 0; i < required_capacity; i += 2) {
            if (num_edge_nodes == 0) {
                out_neighbors_vec.push_back(0);
                out_neighbors_vec.push_back(0);
            } else {
                slide_first_edge_node_to_back();
            }
        }

        // Update startpoint and endpoint
        node_to_startpoint[a] = edge_node_start;
        node_to_endpoint[a] = edge_node_start + required_capacity;

    } else {
        // Space is available. Use it.

        size_t capacity_for_new_node = space.key();
        int old_node = space.value();

        // Split existing node space to accommodate new node.
        extra_space_and_node.erase(capacity_for_new_node, old_node);

        size_t old_node_new_capacity = 0;
        if (old_node == -1) {
            // Update start- and end-points.
            node_to_endpoint[a] = capacity_for_new_node;
            node_to_startpoint[a] = 0;

        } else {
            old_node_new_capacity = (node_to_endpoint[old_node] -
                                      node_to_startpoint[old_node]) -
                                        capacity_for_new_node;
            // Update start- and end-points.
            node_to_endpoint[a] = node_to_endpoint[old_node];
            node_to_startpoint[a] = node_to_endpoint[old_node] - 
                                            capacity_for_new_node;
            node_to_endpoint[old_node] = node_to_startpoint[a];

            endpoint_to_node[node_to_endpoint[old_node]] = old_node;
        }


        // Update the extra capacity information.
        size_t capacity;

        // Moving Node
        if (capacity_for_new_node >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                capacity_for_new_node >= 4 * size_t(out_degrees[a])) {
            capacity = capacity_for_new_node / 2;
            extra_space_and_node.insert(capacity, a);
        }

        // Old Node
        if (old_node != -1) {
            if (old_node_new_capacity >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                   old_node_new_capacity >= 4 * size_t(out_degrees[old_node])) {
                capacity = old_node_new_capacity / 2;
                extra_space_and_node.insert(capacity, old_node);
            }
        }
        // Done updating extra capacity information.
    }

    endpoint_to_node[node_to_endpoint[a]] = a;

    // Now that space has been found, copy the edge info over.
    for (size_t i = 0; i < size_t(out_degrees[a]); i++) {
        move_edge_node_reference(old_startpoint + i, node_to_startpoint[a] + i);
    }

    // Give the old space to another node, if there was space.
    if (old_endpoint > old_startpoint) {

        // extra_space_and_node.erase(old_endpoint - old_startpoint, a);

        if (old_startpoint == 0) {
            // There was no left node. Add the old slot to extra capacity as
            //  node "-1".
            extra_space_and_node.insert(old_endpoint, -1);

            endpoint_to_node.erase(old_endpoint);  // No longer refers to a node
        } else if (extra_space_and_node.contains(old_startpoint, -1)) {
            // There was no left node AND there was empty space to the left.
            extra_space_and_node.erase(old_startpoint, -1);
            extra_space_and_node.insert(old_endpoint, -1);

            endpoint_to_node.erase(old_endpoint);  // No longer refers to a node
        } else {
            // Note that the node could have been moved to the slot immediately
            //  to its left, such that "left node" is the same as node a.
            int left_node = endpoint_to_node[old_startpoint];
            endpoint_to_node.erase(old_startpoint);
            endpoint_to_node[old_endpoint] = left_node;
            node_to_endpoint[left_node] = old_endpoint;

            // Patch in space.
            size_t ln_deg = out_degrees[left_node];
            size_t ln_start = node_to_startpoint[left_node];
            size_t ln_end = node_to_endpoint[left_node];
            size_t ln_old_end = old_startpoint; // Old startpoint of a.

            size_t ln_size = ln_end - ln_start;
            size_t ln_old_size = ln_old_end - ln_start;

            // Determine whether left_node USED to have extra capacity.
            if (ln_old_size >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                    ln_old_size >= 4 * ln_deg) {
                size_t ln_old_capacity = ln_old_size / 2;
                extra_space_and_node.erase(ln_old_capacity, left_node);
            }
            if (ln_size >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                    ln_size >= 4 * ln_deg) {
                size_t ln_capacity = ln_size / 2;
                extra_space_and_node.insert(ln_capacity, left_node);
            }
        }
    }
}


void NTSparseGraph::remove_edge_node_ref(const int main_node,
                                         const size_t ref_loc) {
    size_t last_local_edge_node_loc = node_to_startpoint[main_node] +
                                      size_t(out_degrees[main_node] - 1);
    if (ref_loc < last_local_edge_node_loc) {
        move_edge_node_reference(last_local_edge_node_loc, ref_loc);
    }
    out_degrees[main_node]--;
}

NTPartition NTSparseGraph::nauty_traces_coloring() const {
    NTPartition part = NTPartition(internal_n);

    // Has nodes 0 through internal_n in numeric order
    int* node_ids_ptr = part.get_node_ids();

    // All 1's except for the last entry
    int* partition_ints_ptr = part.get_partition_ints();

    // Separate the regular nodes from the edge nodes.
    partition_ints_ptr[n - 1] = 0;

    // Separate out the nodes with self-loops from the others.
    // Put the nodes without self-loops first.
    if (num_self_loops > 0) {
        size_t normal_idx = 0;
        size_t self_loop_idx = n - num_self_loops;
        partition_ints_ptr[self_loop_idx - 1] = 0;
        for (int i = 0; i < int(n); i++) {
            if (has_self_loop[i]) {
                node_ids_ptr[self_loop_idx] = i;
                self_loop_idx++;
            } else {
                node_ids_ptr[normal_idx] = i;
                normal_idx++;
            }
        }
    }

    if (!directed) {
        return part;
    }

    // Remember that when a connects to b, then under the hood there are two
    //  edge nodes in an undirected chain linking a -- x -- y -- b.
    // We need to distinguish between two kinds of edge nodes:
    // those that correspond to the head of a directed edge and those that don't
    size_t front_idx = n;
    size_t back_idx = internal_n - 1;
    for (int edge_node = n; edge_node < int(internal_n); edge_node++) {
        const Edge& e = edge_node_to_edge.find(edge_node)->second;
        if (_out_neighbors[e.first].find(e.second) ==
                                            _out_neighbors[e.first].end()) {
            // Edge (first, second) is not in the graph.
            node_ids_ptr[back_idx] = edge_node;
            back_idx--;
        } else {
            node_ids_ptr[front_idx] = edge_node;
            front_idx++;
        }
    }
    partition_ints_ptr[front_idx - 1] = 0;

    return part;
}

NTPartition NTSparseGraph::nauty_traces_coloring(
                                    const Coloring<int> &node_colors) const {
    if (node_colors.size() != n) {
        throw std::range_error(std::string("Error! `node_colors` has ")
                               + std::to_string(node_colors.size())
                               + " nodes whereas the graph has "
                               + std::to_string(n) + " nodes.");
    }

    NTPartition part = NTPartition(internal_n);

    // Has nodes 0 through internal_n in numeric order
    int* node_ids_ptr = part.get_node_ids();

    // All 1's except for the last entry
    int* partition_ints_ptr = part.get_partition_ints();

    // Separate the regular nodes from the edge nodes.
    partition_ints_ptr[n - 1] = 0;

    // Within the regular nodes, list them in order by color, and within a
    //  given color, list the nodes with self-loops first.
    size_t num_nodes_listed = 0;
    size_t cell_start;
    size_t next_self_loop_node_idx;
    for (auto color_itr = node_colors.colors().begin();
              color_itr != node_colors.colors().end(); color_itr++) {

        cell_start = num_nodes_listed;
        next_self_loop_node_idx = num_nodes_listed;

        for (auto node_itr = node_colors.cell(*color_itr).begin();
                  node_itr != node_colors.cell(*color_itr).end(); node_itr++) {
            if (*node_itr < 0 || *node_itr >= int(n)) {
                throw std::range_error(
                        std::string("Error! Element of `node_colors` "
                                    + std::to_string(*node_itr)
                                    + " is not one of the graph's nodes."));
            }
            if (has_self_loop[*node_itr]) {
                if (next_self_loop_node_idx < num_nodes_listed) {
                    node_ids_ptr[num_nodes_listed] =
                                    node_ids_ptr[next_self_loop_node_idx];
                }
                node_ids_ptr[next_self_loop_node_idx] = *node_itr;
                next_self_loop_node_idx++;
            } else {
                node_ids_ptr[num_nodes_listed] = *node_itr;
            }
            num_nodes_listed++;
        }
        partition_ints_ptr[num_nodes_listed - 1] = 0;
        if (next_self_loop_node_idx > cell_start) {
            partition_ints_ptr[next_self_loop_node_idx - 1] = 0;
        }
    }

    if (!directed) {
        return part;
    }

    // Remember that when a connects to b, then under the hood there are two
    //  edge nodes in an undirected chain linking a -- x -- y -- b.
    // We need to distinguish between two kinds of edge nodes:
    // those that correspond to the head of a directed edge and those that don't
    size_t front_idx = n;
    size_t back_idx = internal_n - 1;
    for (int edge_node = n; edge_node < int(internal_n); edge_node++) {
        const Edge& e = edge_node_to_edge.find(edge_node)->second;
        if (_out_neighbors[e.first].find(e.second) ==
                                            _out_neighbors[e.first].end()) {
            // Edge (first, second) is not in the graph.
            node_ids_ptr[back_idx] = edge_node;
            back_idx--;
        } else {
            node_ids_ptr[front_idx] = edge_node;
            front_idx++;
        }
    }
    partition_ints_ptr[front_idx - 1] = 0;

    return part;
}

NTPartition NTSparseGraph::nauty_traces_coloring(
                            const Coloring<Edge, EdgeHash> &edge_colors) const {
    if (edge_colors.size() != m) {
        throw std::range_error(std::string("Error! `edge_colors` has ")
                               + std::to_string(edge_colors.size())
                               + " edges whereas the graph has "
                               + std::to_string(m) + " edges.");
    }

    NTPartition part = NTPartition(internal_n);

    // Has nodes 0 through internal_n in numeric order
    int* node_ids_ptr = part.get_node_ids();

    // All 1's except for the last entry
    int* partition_ints_ptr = part.get_partition_ints();

    // Separate the regular nodes from the edge nodes.
    partition_ints_ptr[n - 1] = 0;

    // Separate out the nodes with self-loops from the others.
    // Put the nodes without self-loops first.
    //
    // We will place the nodes with self-loops during the edges part of the
    //  computation.
    if (num_self_loops > 0) {
        size_t normal_idx = 0;
        size_t self_loop_idx = n - num_self_loops;
        partition_ints_ptr[self_loop_idx - 1] = 0;
        for (int i = 0; i < int(n); i++) {
            if (has_self_loop[i]) {
                continue;
            }
            node_ids_ptr[normal_idx] = i;
            normal_idx++;
        }
    }

    if (!directed) {
        size_t next_idx = n;
        size_t self_loop_idx = n - num_self_loops;
        int edge_node;
        for (auto color_itr = edge_colors.colors().begin();
                    color_itr != edge_colors.colors().end(); color_itr++) {
            for (auto edge_itr = edge_colors.cell(*color_itr).begin();
                   edge_itr != edge_colors.cell(*color_itr).end(); edge_itr++) {
                if (edge_itr->first == edge_itr->second) {
                    // A self-loop
                    if (!has_self_loop[edge_itr->first]) {
                        throw std::range_error(std::string("Error! Node ")
                                               + std::to_string(edge_itr->first)
                                               + " does not have a self-loop.");
                    }
                    node_ids_ptr[self_loop_idx] = edge_itr->first;
                    self_loop_idx++;
                    continue;
                }
                // A regular edge
                auto en_itr = edge_to_edge_node.find(*edge_itr);
                if (en_itr == edge_to_edge_node.end()) {
                    throw std::range_error(
                            std::string("Error! Edge (")
                            + std::to_string(edge_itr->first) + ", "
                            + std::to_string(edge_itr->second)
                            + ") is in the coloring but not the graph.");
                }
                edge_node = en_itr->second;
                node_ids_ptr[next_idx] = edge_node;
                next_idx++;
            }
            if (next_idx > 0) {
                partition_ints_ptr[next_idx - 1] = 0;
            }
            if (self_loop_idx > 0) {
                partition_ints_ptr[self_loop_idx - 1] = 0;
            }
        }
        return part;
    }

    // directed
    size_t next_edge_idx = n;
    size_t next_non_edge_idx =
                internal_n - (num_edge_nodes - (m - num_self_loops));
    size_t self_loop_idx = n - num_self_loops;

    int edge_node;
    for (auto color_itr = edge_colors.colors().begin();
             color_itr != edge_colors.colors().end(); color_itr++) {
        for (auto edge_itr = edge_colors.cell(*color_itr).begin();
                 edge_itr != edge_colors.cell(*color_itr).end(); edge_itr++) {
            if (edge_itr->first == edge_itr->second) {
                // A self-loop
                if (!has_self_loop[edge_itr->first]) {
                    throw std::range_error(std::string("Error! Node ")
                                           + std::to_string(edge_itr->first)
                                           + " does not have a self-loop.");
                }
                node_ids_ptr[self_loop_idx] = edge_itr->first;
                self_loop_idx++;
                continue;
            }
            // A regular edge
            auto en_itr = edge_to_edge_node.find(*edge_itr);
            if (en_itr == edge_to_edge_node.end() ||
                    (_out_neighbors[edge_itr->first].find(edge_itr->second) ==
                     _out_neighbors[edge_itr->first].end())) {
                throw std::range_error(
                        std::string("Error! Edge (")
                        + std::to_string(edge_itr->first) + ", "
                        + std::to_string(edge_itr->second)
                        + ") is in the coloring but not the graph.");
            }
            edge_node = en_itr->second;
            node_ids_ptr[next_edge_idx] = edge_node;
            next_edge_idx++;

            if (_out_neighbors[edge_itr->second].find(edge_itr->first) ==
                    _out_neighbors[edge_itr->second].end()) {
                // The reverse edge is not in the graph. Put the
                //  corresponding edge node at the end.
                edge_node = out_neighbors_vec[node_to_startpoint[edge_node]+1];
                node_ids_ptr[next_non_edge_idx] = edge_node;
                next_non_edge_idx++;
            }
        }
        if (next_edge_idx > 0) {
            partition_ints_ptr[next_edge_idx - 1] = 0;
        }
        if (self_loop_idx > 0) {
            partition_ints_ptr[self_loop_idx - 1] = 0;
        }
    }

    return part;
}

NTPartition NTSparseGraph::nauty_traces_coloring(
                                    const Coloring<int> &node_colors,
                                    const Coloring<Edge, EdgeHash>
                                                        &edge_colors) const {
    if (node_colors.size() != n) {
        throw std::range_error(std::string("Error! `node_colors` has ")
                               + std::to_string(node_colors.size())
                               + " nodes whereas the graph has "
                               + std::to_string(n) + " nodes.");
    }
    if (edge_colors.size() != m) {
        throw std::range_error(std::string("Error! `edge_colors` has ")
                               + std::to_string(edge_colors.size())
                               + " edges whereas the graph has "
                               + std::to_string(m) + " edges.");
    }

    NTPartition part = NTPartition(internal_n);

    // Has nodes 0 through internal_n in numeric order
    int* node_ids_ptr = part.get_node_ids();

    // All 1's except for the last entry
    int* partition_ints_ptr = part.get_partition_ints();

    // Separate the regular nodes from the edge nodes.
    partition_ints_ptr[n - 1] = 0;

    // Keep track of where we put nodes with self-loops so that they can be
    //  further refined when we get to the self-loop edge colors.
    //
    // Maps a node with a self-loop to the start of the partition containing
    //  that node.
    std::unordered_map<int, size_t> self_loop_node_to_slot_start =
                                        std::unordered_map<int, size_t>();

    // Within the regular nodes, list them in order by color, and within a
    //  given color, list the nodes with self-loops first.
    size_t num_nodes_listed = 0;
    size_t cell_start;
    size_t next_self_loop_node_idx;
    size_t self_loop_slot_start;
    for (auto color_itr = node_colors.colors().begin();
              color_itr != node_colors.colors().end(); color_itr++) {

        cell_start = num_nodes_listed;
        self_loop_slot_start = num_nodes_listed;
        next_self_loop_node_idx = self_loop_slot_start;

        for (auto node_itr = node_colors.cell(*color_itr).begin();
                  node_itr != node_colors.cell(*color_itr).end(); node_itr++) {
            if (*node_itr < 0 || *node_itr >= int(n)) {
                throw std::range_error(
                        std::string("Error! Element of `node_colors` "
                                    + std::to_string(*node_itr)
                                    + " is not one of the graph's nodes."));
            }
            if (has_self_loop[*node_itr]) {
                if (next_self_loop_node_idx < num_nodes_listed) {
                    node_ids_ptr[num_nodes_listed] =
                                    node_ids_ptr[next_self_loop_node_idx];
                }
                node_ids_ptr[next_self_loop_node_idx] = *node_itr;
                self_loop_node_to_slot_start[*node_itr] = self_loop_slot_start;
                next_self_loop_node_idx++;
            } else {
                node_ids_ptr[num_nodes_listed] = *node_itr;
            }
            num_nodes_listed++;
        }
        partition_ints_ptr[num_nodes_listed - 1] = 0;
        if (next_self_loop_node_idx > cell_start) {
            partition_ints_ptr[next_self_loop_node_idx - 1] = 0;
        }
    }

    ////////////// Now Edge Colors //////////////

    // For a given collection of nodes that have the same node color and all
    //  have self-loops, we have recorded their "self loop start idx," which is
    //  the idx of the first node in that collection.
    //
    // Now, as we go through edge colors on those 
    std::unordered_map<size_t, size_t> sl_start_to_curr_idx =
                                           std::unordered_map<size_t, size_t>();
    std::unordered_set<size_t> referenced_sl_start_locs;
    if (!directed) {
        size_t next_idx = n;
        size_t self_loop_cell_curr_idx;
        size_t self_loop_start_idx;
        int edge_node;
        int a, b;
        for (auto color_itr = edge_colors.colors().begin();
                    color_itr != edge_colors.colors().end(); color_itr++) {
            referenced_sl_start_locs = std::unordered_set<size_t>();
            for (auto edge_itr = edge_colors.cell(*color_itr).begin();
                   edge_itr != edge_colors.cell(*color_itr).end(); edge_itr++) {
                a = edge_itr->first;
                b = edge_itr->second;
                if (a == b) {
                    // A self-loop
                    if (!has_self_loop[a]) {
                        throw std::range_error(std::string("Error! Node ")
                                               + std::to_string(a)
                                               + " does not have a self-loop.");
                    }
                    self_loop_start_idx = self_loop_node_to_slot_start[a];
                    self_loop_cell_curr_idx =
                        sl_start_to_curr_idx[self_loop_start_idx];
                    node_ids_ptr[self_loop_cell_curr_idx] = a;
                    sl_start_to_curr_idx[self_loop_start_idx]++;
                    referenced_sl_start_locs.insert(self_loop_start_idx);
                    continue;
                }
                // A regular edge
                auto en_itr = edge_to_edge_node.find(*edge_itr);
                if (en_itr == edge_to_edge_node.end()) {
                    throw std::range_error(
                            std::string("Error! Edge (")
                            + std::to_string(edge_itr->first) + ", "
                            + std::to_string(edge_itr->second)
                            + ") is in the coloring but not the graph.");
                }
                edge_node = en_itr->second;
                node_ids_ptr[next_idx] = edge_node;
                next_idx++;
            }
            if (next_idx > 0) {
                partition_ints_ptr[next_idx - 1] = 0;
            }
            for (auto idx_itr = referenced_sl_start_locs.begin();
                      idx_itr != referenced_sl_start_locs.end(); idx_itr++) {
                self_loop_cell_curr_idx = sl_start_to_curr_idx[*idx_itr];
                // We know self_loop_cell_curr_idx > 0 bec. it was incremented
                partition_ints_ptr[self_loop_cell_curr_idx - 1] = 0;
            }
        }
        return part;
    }

    // directed
    size_t next_edge_idx = n;
    size_t next_non_edge_idx =
                internal_n - (num_edge_nodes - (m - num_self_loops));
    size_t self_loop_cell_curr_idx;
    size_t self_loop_start_idx;

    int edge_node;
    int a, b;
    for (auto color_itr = edge_colors.colors().begin();
             color_itr != edge_colors.colors().end(); color_itr++) {
        referenced_sl_start_locs = std::unordered_set<size_t>();
        for (auto edge_itr = edge_colors.cell(*color_itr).begin();
                 edge_itr != edge_colors.cell(*color_itr).end(); edge_itr++) {
            a = edge_itr->first;
            b = edge_itr->second;
            if (a == b) {
                // A self-loop
                if (!has_self_loop[a]) {
                    throw std::range_error(std::string("Error! Node ")
                                           + std::to_string(a)
                                           + " does not have a self-loop.");
                }
                self_loop_start_idx = self_loop_node_to_slot_start[a];
                self_loop_cell_curr_idx =
                    sl_start_to_curr_idx[self_loop_start_idx];
                node_ids_ptr[self_loop_cell_curr_idx] = a;
                sl_start_to_curr_idx[self_loop_start_idx]++;
                referenced_sl_start_locs.insert(self_loop_start_idx);
                continue;
            }
            // A regular edge
            auto en_itr = edge_to_edge_node.find(*edge_itr);
            if (en_itr == edge_to_edge_node.end() ||
                    (_out_neighbors[edge_itr->first].find(edge_itr->second) ==
                     _out_neighbors[edge_itr->first].end())) {
                throw std::range_error(
                        std::string("Error! Edge (")
                        + std::to_string(edge_itr->first) + ", "
                        + std::to_string(edge_itr->second)
                        + ") is in the coloring but not the graph.");
            }
            edge_node = en_itr->second;
            node_ids_ptr[next_edge_idx] = edge_node;
            next_edge_idx++;

            if (_out_neighbors[edge_itr->second].find(edge_itr->first) ==
                    _out_neighbors[edge_itr->second].end()) {
                // The reverse edge is not in the graph. Put the
                //  corresponding edge node at the end.
                edge_node = out_neighbors_vec[node_to_startpoint[edge_node]+1];
                node_ids_ptr[next_non_edge_idx] = edge_node;
                next_non_edge_idx++;
            }
        }
        if (next_edge_idx > 0) {
            partition_ints_ptr[next_edge_idx - 1] = 0;
        }
        for (auto idx_itr = referenced_sl_start_locs.begin();
                  idx_itr != referenced_sl_start_locs.end(); idx_itr++) {
            self_loop_cell_curr_idx = sl_start_to_curr_idx[*idx_itr];
            // We know self_loop_cell_curr_idx > 0 bec. it was incremented
            partition_ints_ptr[self_loop_cell_curr_idx - 1] = 0;
        }
    }

    return part;
}
