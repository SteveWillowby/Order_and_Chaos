/* Provides a C++ graph class that interfaces with the Nauty & Traces sparse
 *  graph classes.
 *
 * Extends the Graph class defined in graph.h
 *
 * This class prioritizes efficient asymptotic (i.e. big-O) runtime at the
 *  expense of data storage and constant factors.
 *
 *
 * Undirected graphs can have 1 <= n + m <= 2,000,000,000.
 * Directed graphs are limited to 1 <= n + 2*m_undir <= 2,000,000,000 where
 *  m_undir is the number of edges in the undirected version of the graph.
 *
 *  These limits exist for the following reasons.
 *  2,000,000,000 is the maximum number of nodes that the Nauty and Traces code
 *      allows.
 *  We add 1 node per edge on undirected graphs so that the edges can be
 *      colored.
 *  We add 2 nodes per (unidrected) edge on directed graphs so that even if the
 *      edges are not colored, the code can use colors behind the scenes to run
 *      Traces on the graph. Traces does not natively support a directed input.
 *      Additionally, these extra nodes do allow for colorings on the edges.
 *  The extra nodes are entirely behind the scenes. For example, they do not
 *      change the output of num_nodes(). They only appear in the output of
 *      functions with "nauty_traces" in their name.
 */

#include "nauty27r4/nauty.h"
#include "nauty27r4/nausparse.h"

#include "augmented_multimap.h"
// #include "coloring.h"
#include "edge.h"
#include "sparse_graph.h"

#include<unordered_map>
#include<utility>
#include<vector>

#ifndef SYM__NT_SPARSE_GRAPH_H
#define SYM__NT_SPARSE_GRAPH_H

const size_t NAUTY_TRACES_MAXN = 2000000000;

class NTSparseGraph : public SparseGraph {

public:
    NTSparseGraph(const bool directed);
    // n is the number of nodes
    NTSparseGraph(const bool directed, const size_t n);
    NTSparseGraph(const Graph &g);

    // returns a `sparsegraph` struct that can be passed into nauty or traces
    // virtual const sparsegraph as_nauty_traces_graph() const;

    // Returns a coloring for the Nauty/Traces sparsegraph
    /*
    const Coloring<int> &
        nauty_traces_coloring(const Coloring<int> &node_coloring) const;
    const Coloring<int> &
        nauty_traces_coloring(const Coloring<Edge> &edge_coloring) const;
    const Coloring<int> &
        nauty_traces_coloring(const Coloring<int> &node_coloring,
                              const Coloring<Edge> &edge_coloring) const;
    */

    // size_t num_nodes() const; -- defined in graph.cpp
    // size_t num_edges() const; -- defined in graph.cpp

    // const bool directed; -- defined in graph.h

    // Returns the id of the new node.
    // O(1)
    virtual int add_node();
    // Deletes node a AND if a < n-1, relabels node n-1 to have label a.
    //  Returns the old label of the node that is now labeled a.
    // O(number of node a's neighbors + number of node n-1's neighbors)
    virtual int delete_node(const int a);

    //////////////////////////////////////////////////////////////
    // All edge-changing functions have amortized O(1) runtime. //
    //////////////////////////////////////////////////////////////

    // Returns true iff the edge was absent (and thus now added)
    virtual bool add_edge(const int a, const int b);
    // Returns true iff the edge was present (and thus now deleted)
    virtual bool delete_edge(const int a, const int b);

    // Deletes the edge if it exists and adds it if it does not.
    virtual void flip_edge(const int a, const int b);

    //  defined in sparse_graph.h/cpp
    // virtual bool has_edge(const int a, const int b) const;

    //  O(1) -- defined in sparse_graph.cpp
    // virtual const std::unordered_set<int> &neighbors(const int a) const;
    //  neighbors that node a points to
    //  O(1) -- defined in sparse_graph.cpp
    // virtual const std::unordered_set<int> &out_neighbors(const int a) const;
    //  neighbors that point to node a
    //  O(1) -- defined in sparse_graph.cpp
    // virtual const std::unordered_set<int> &in_neighbors(const int a) const;

protected:
    // size_t n; -- defined in graph.h
    size_t internal_n; // TODO: Check if actually used.
    // size_t m; -- defined in graph.h
    size_t num_edge_nodes; // TODO: Check if actually used.

    //  defined in sparse_graph.h
    // std::vector<std::unordered_set<int>> _neighbors;
    //  defined in sparse_graph.h
    //  Unused if graph is undirected
    // std::vector<std::unordered_set<int>> _out_neighbors;
    //  defined in sparse_graph.h
    //  Unused if graph is undirected
    // std::vector<std::unordered_set<int>> _in_neighbors;

    // The following functions edit the Nauty/Traces graph representation.
    /*
    virtual void _add_node();
    virtual int _delete_node(const int a);
    virtual void _add_edge(const int a, const int b);
    virtual void _delete_edge(const int a, const int b);
    */

    // Stores the ID of a node corresponding to edge (a, b).
    //  In a directed graph, each (undirected) edge really has two nodes. To
    //  access the second node, access the node for edge (b, a).
    std::unordered_map<Edge, int, EdgeHash> edge_to_edge_node;
    std::unordered_map<int, Edge> edge_node_to_edge;

    // We require that the nodes of the actual graph get the first n labels.
    //  Then the edge nodes get the remaining labels.
    std::vector<int> out_degrees;
    std::vector<int> node_to_startpoint;
    std::vector<int> node_to_endpoint;

    // Gives the two places in out_neighbors_vec where the internal node is
    //  listed as a neighbor.
    // When the graph is directed, we follow the
    //  convention that the first listed node is the real node and the second
    //  is the other edge node.
    std::unordered_map<int, std::pair<size_t,size_t>> edge_node_to_places;

    // Used for regular nodes.
    std::unordered_map<int, int> endpoint_to_node;

    // Only real nodes have extra space. Edge nodes always have a fixed amount.

    // A node is defined to have extra space when the following are all met:
    //  * (endpoint - startpoint) >= 2 * MIN_EDGE_SPACE_PER_NODE
    //  * out_degree <= (endpoint - startpoint) / 4
    //
    // We always give FLOOR(half the space) to the node that needs space.
    //  Further, we require that the new space have at least twice as much space
    //  as the node receiving that space needs.

    // std::vector<bool> has_extra_space; -- currently unused
    // std::vector<int> extra_space; -- currently unused.

    // Call lower_bound(target) to get the node with the smallest capacity that
    //  is at least as large as your target storage space.

    AugmentedMultimap<size_t, int> extra_space_and_node;

    // When a node is an edge node and the graph is directed, we follow the
    //  convention that the first listed node is the real node and the second
    //  is the other edge node.
    std::vector<int> out_neighbors_vec;

    // Must be an even number.
    const size_t MIN_EDGE_SPACE_PER_NODE = 8;

    int allocate_edge_node();
    void relabel_edge_node(const int a, const int b);
    void slide_first_edge_node_to_back();

    // Can be used for new nodes as well that currently have no space and no
    //  edges. (i.e. endpoint = startpoint)
    // Only for regular nodes.
    void move_node_to_more_space(const int a);

    void move_edge_node(const size_t init_loc, const size_t target_loc);
};

#endif
