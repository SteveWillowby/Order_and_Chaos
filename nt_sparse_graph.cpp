#include "augmented_multimap.h"
// #include "coloring.h"
#include "edge.h"
#include "sparse_graph.h"
#include "nt_sparse_graph.h"

#include<unordered_map>
#include<utility>
#include<vector>

// TODO: Remove
#include<iostream>

NTSparseGraph::NTSparseGraph(const bool directed) : NTSparseGraph(directed,1) {}

NTSparseGraph::NTSparseGraph(const bool directed, size_t n)
                    : SparseGraph(directed, n) {

    internal_n = n;
    num_edge_nodes = 0;

    out_degrees = std::vector<int>(n, 0);
    has_self_loop = std::vector<bool>(n, false);

    // Begin with MIN_EDGE_SPACE_PER_NODE edge slots per node.
    node_to_startpoint = std::vector<int>();
    node_to_endpoint = std::vector<int>();
    endpoint_to_node = std::unordered_map<int, int>();
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
NTSparseGraph::NTSparseGraph(const Graph &g) :
        NTSparseGraph(g.directed, g.num_nodes()) {

    for (int i = 0; i < int(n); i++) {
        for (auto nbr = g.out_neighbors(i).begin();
                                    nbr != g.out_neighbors(i).end(); nbr++) {
            add_edge(i, *nbr);
        }
    }

    /*
    // This here is partially finished code for making this constructor more
    //  efficient.
    //  
    // SparseGraph's constructor does not make any calls to functions like add_node
    //  or add_edge. Thus we can call it and then fill out the NT info separately.
    if (directed) {
        num_edge_nodes = 2 * m;
        internal_n = n + num_edge_nodes;

        if (internal_n > NAUTY_TRACES_MAXN) {
            throw std::logic_error(
               std::string("Error! Too many nodes for Nauty/Traces. Directed") +
               std::string(" graphs must have N + 2 * (undirected) M <= ") +
               std::to_string(NAUTY_TRACES_MAXN));
        }

        size_t total_space_needed = 0;

        // Add basic node info.
        out_degrees = std::vector<int>(internal_n, 0);
        node_to_startpoint = std::vector<int>(internal_n, 0);
        node_to_endpoint = std::vector<int>(internal_n, 0);
        endpoint_to_node = std::unordered_map<int, int>();

        for (int i = 0; i < n; i++) {
            out_degrees[i] = _neighbors[i].size();
            node_to_startpoint[i] = total_space_needed;
            total_space_needed +=
                (out_degrees[i] < MIN_EDGE_SPACE_PER_NODE ?
                    MIN_EDGE_SPACE_PER_NODE : out_degrees[i]);
            node_to_endpoint[i] = total_space_needed;
            endpoint_to_node[total_space_needed] = i;
        }

        size_t edge_node_start = total_space_needed;
        total_space_needed += num_edge_nodes * 2;
        out_neighbors_vec = std::vector<int>(total_space_needed, 0);

        // Add edge and edge node info.

        std::vector<int> slots_used = std::vector<int>(n, 0);
        int next_edge_node = n;

        for (int i = 0; i < n; i++) {
            for (auto nbr = _neighbors[i].begin();
                            nbr < _neighbors[i].end(); nbr++) {
                if (*nbr < i) {
                    continue;
                }

                out_neighbors_vec[node_to_startpoint[i] + slots_used[i]] =
                                        next_edge_node;
                slots_used[i]++;
                out_neighbors_vec[node_to_startpoint[i]
                
            }
        }
    } else {

    }
    */
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
    int replacement_node = SparseGraph::delete_node(a);

    bool actually_replaced = replacement_node != int(has_self_loop.size() - 1);
    if (actually_replaced) {
        has_self_loop[a] = has_self_loop[replacement_node];
    }
    has_self_loop.pop_back();

    // TODO: In an undirected graph, the relabeling of a node might mean that
    //  for some edge_node_to_places pairs, the order of the two places must be
    //  swapped if the order of the places' node ids have been swapped.

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

    if (out_degrees[a] == (node_to_endpoint[a] - node_to_startpoint[a])) {
        // Move node's edge info to a new place.
        move_node_to_more_space(a);
    }
    if (out_degrees[b] == (node_to_endpoint[b] - node_to_startpoint[b])) {
        // Move node's edge info to a new place.
        move_node_to_more_space(b);
    }

    if (directed) {
        internal_n += 2;
        num_edge_nodes += 2;

        int edge_node_A = allocate_edge_node();
        int edge_node_B = allocate_edge_node();

        out_neighbors_vec[node_to_startpoint[a] + out_degrees[a]] = edge_node_A;
        out_neighbors_vec[node_to_startpoint[b] + out_degrees[b]] = edge_node_B;
        out_neighbors_vec[node_to_startpoint[edge_node_A]] = a;
        out_neighbors_vec[node_to_startpoint[edge_node_A] + 1] = edge_node_B;
        out_neighbors_vec[node_to_startpoint[edge_node_B]] = b;
        out_neighbors_vec[node_to_startpoint[edge_node_B] + 1] = edge_node_A;

        edge_node_to_places[edge_node_A] =
            std::pair<size_t, size_t>(node_to_startpoint[a] + out_degrees[a],
                                      node_to_startpoint[edge_node_B] + 1);
        edge_node_to_places[edge_node_B] =
            std::pair<size_t, size_t>(node_to_startpoint[b] + out_degrees[b],
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
        out_neighbors_vec[node_to_startpoint[a] + out_degrees[a]] = edge_node;
        out_neighbors_vec[node_to_startpoint[b] + out_degrees[b]] = edge_node;

        // Put the smaller node ID first.
        out_neighbors_vec[node_to_startpoint[edge_node]] = (a < b ? a : b);
        out_neighbors_vec[node_to_startpoint[edge_node] + 1] = (a < b ? b : a);

        edge_to_edge_node[EDGE(a, b, directed)] = edge_node;
        edge_node_to_edge[edge_node] = EDGE(a, b, directed);

        int min_node = (a < b ? a : b);
        int max_node = (a < b ? b : a);
        edge_node_to_places[edge_node] = std::pair<size_t, size_t>(
                    node_to_startpoint[min_node] + out_degrees[min_node],
                    node_to_startpoint[max_node] + out_degrees[max_node]);

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

    if (a == b) {
        has_self_loop[a] = false;
        return true;
    }

    if (directed &&
            _out_neighbors[b].find(a) != _out_neighbors[b].end()) {
        // The edge was deleted _but_ we need the same edge nodes for the
        //  reverse edge.
        return true;
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
                size_t(capacity_A) >= MIN_EDGE_SPACE_PER_NODE * 2) {
        // Just acquired enough extra space.
        extra_capacity = capacity_A / 2;
        extra_space_and_node.insert(extra_capacity, a);
    }
    if (size_t(out_degrees[b]) * 4 <= capacity_B &&
            size_t(out_degrees[b] + 1) * 4 > capacity_B &&
                size_t(capacity_B) >= MIN_EDGE_SPACE_PER_NODE * 2) {
        // Just acquired enough extra space.
        extra_capacity = capacity_B / 2;
        extra_space_and_node.insert(extra_capacity, b);
    }

    return true;
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
    /*
    std::cout<<"The locations for node "<<a<<" are thought to be "<<locations.first<<" and ";
    std::cout<<locations.second<<". At those places we find the values ";
    std::cout<<out_neighbors_vec[locations.first]<<" and "<<out_neighbors_vec[locations.second]<<"."<<std::endl;
    std::cout<<"Node "<<b<<" will now be listed as having those locations: ";
    */

    // Update out_neighbors_vec
    out_neighbors_vec[locations.first] = b;
    out_neighbors_vec[locations.second] = b;

    // Update edge_node_to_places
    edge_node_to_places.erase(a);
    edge_node_to_places[b] = locations;

    // std::cout<<edge_node_to_places[b].first<<" "<<edge_node_to_places[b].second<<std::endl;

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
    int new_startpoint = node_to_startpoint[edge_node_of_slot];

    int old_endpoint = out_neighbors_vec.size();
    int old_startpoint = old_endpoint - 2;

    int moving_node = endpoint_to_node[old_endpoint];

    int largest_node = out_degrees.size() - 1;
    int largest_node_startpoint = node_to_startpoint[largest_node];
    int largest_node_endpoint = largest_node_startpoint + 2;

    // TODO: remove ALL print statements
    // std::cout<<"Old. vs. New startpoint: "<<old_startpoint<<" "<<new_startpoint<<std::endl;
    // std::cout<<"###"<<moving_node<<" moves to "<<edge_node_of_slot<<"'s space."<<std::endl;

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

        // std::cout<<"Entering the gloom. Relabeling "<<largest_node<<" as "<<edge_node_of_slot<<std::endl;
        relabel_edge_node(largest_node, edge_node_of_slot);
    } else if (edge_node_of_slot > largest_node + 1) {
        // TODO: Remove this check and line.
        std::cout<<"Logic Error! size_t(edge_node_of_slot) != out_degrees.size())"<<std::endl;
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

    int old_startpoint = node_to_startpoint[a];
    int old_endpoint = node_to_endpoint[a];

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
        int capacity;

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
    for (int i = 0; i < out_degrees[a]; i++) {
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
            int ln_start = node_to_startpoint[left_node];
            int ln_end = node_to_endpoint[left_node];
            int ln_old_end = old_startpoint; // Old startpoint of a.

            size_t ln_size = ln_end - ln_start;
            size_t ln_old_size = ln_old_end - ln_start;

            // Determine whether left_node USED to have extra capacity.
            if (ln_old_size >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                    ln_old_size >= 4 * ln_deg) {
                int ln_old_capacity = ln_old_size / 2;
                extra_space_and_node.erase(ln_old_capacity, left_node);
            }
            if (ln_size >= 2 * MIN_EDGE_SPACE_PER_NODE &&
                    ln_size >= 4 * ln_deg) {
                int ln_capacity = ln_size / 2;
                extra_space_and_node.insert(ln_capacity, left_node);
            }
        }
    }
}


void NTSparseGraph::remove_edge_node_ref(const int main_node,
                                         const size_t ref_loc) {
    size_t last_local_edge_node_loc = node_to_startpoint[main_node] +
                                      out_degrees[main_node] - 1;
    if (ref_loc < last_local_edge_node_loc) {
        move_edge_node_reference(last_local_edge_node_loc, ref_loc);
    }
    out_degrees[main_node]--;
}
