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

#include "nauty27r4_modified/nauty.h"
#include "nauty27r4_modified/nausparse.h"

#include "augmented_multimap.h"
#include "coloring.h"
#include "nt_partition.h"
#include "edge.h"
#include "sparse_graph.h"

#include<cstddef>
#include<unordered_map>
#include<utility>
#include<vector>

#ifndef SCHENO__NT_SPARSE_GRAPH_H
#define SCHENO__NT_SPARSE_GRAPH_H

// Only define full debug mode to fully test out the class in nt_graph_test
//  When this flag is defined, the protected members of the class become public.
//  Comment out for normal use.
// #define SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE

const size_t NAUTY_TRACES_MAXN = 2000000000;

class NTSparseGraph : public SparseGraph {

public:
    NTSparseGraph(const bool directed);
    // n is the number of nodes
    NTSparseGraph(const bool directed, const size_t n);
    NTSparseGraph(const Graph &g);

    // Copy assignment operators.
    NTSparseGraph& operator=(const Graph& g);
    NTSparseGraph& operator=(const SparseGraph& g);
    NTSparseGraph& operator=(const NTSparseGraph& g);

    // returns a `sparsegraph` struct that can be passed into nauty or traces
    //
    // NOTE: This function is not a const function because the sparsegraph
    //  struct type requires non-const pointers to the data. However, Nauty and
    //  Traces do not modify the graph structure, so as long as you also do not
    //  use this info to modify the graph, then this function is safe.
    //
    // NOTE: nauty_traces.cpp assumes details about this class's output for this
    //  function. If the details of this output change, traces.cpp must change
    //  too.
    //  In particular, nauty_traces.cpp assumes the following:
    //      * The first n nodes in v are the regular nodes.
    //      * There are at most 3 nt edges per edge in the graph when the graph
    //          is directed and at most 2 nt edges per edge when the undirected.
    const sparsegraph as_nauty_traces_graph();

    // Returns a partitioning (i.e. coloring) for the Nauty/Traces sparsegraph.
    //
    // NOTE: traces.cpp assumes the following about all partition outputs:
    //      * The real nodes are the first n nodes in the node_ids list.
    NTPartition nauty_traces_coloring() const;
    NTPartition nauty_traces_coloring(const Coloring<int> &node_colors) const;
    NTPartition nauty_traces_coloring(const Coloring<Edge, EdgeHash>
                                                            &edge_colors) const;
    NTPartition nauty_traces_coloring(const Coloring<int> &node_colors,
                                      const Coloring<Edge, EdgeHash>
                                                            &edge_colors) const;

    /*

    // Featured colorings:

    // When on, this class keeps info ready to more quickly make an NTPartition
    //  that separates out edge nodes from regular nodes.
    //
    // NOTE: Can only be turned on if there is not already a featured coloring.
    void turn_on_structure_coloring();

    // When on, this class keeps info ready to more quickly make an NTPartition
    //  that highlights any edges (or non-edges) which have been changed.
    //
    // NOTE: when a change highlights coloring is set, the graph cannot add or
    //  remove any nodes. Only edges may be changed.
    //
    // NOTE: Can only be turned on if there is not already a featured coloring.
    void turn_on_change_highlights_coloring();
    void turn_on_change_highlights_coloring(const Coloring<int> &node_colors);
    void turn_on_change_highlights_coloring(const Coloring<Edge, EdgeHash>
                                                                &edge_colors);
    void turn_on_change_highlights_coloring(const Coloring<int> &node_colors,
                                            const Coloring<Edge, EdgeHash>
                                                                &edge_colors);

    void turn_off_featured_coloring();

    // If structure coloring is set, returns the structure coloring. 
    // If a change highlights coloring is set, returns that.
    // Otherwise, throws an error.
    //
    // Note that this still requires a memcopy of O(n + m) values. However, it
    //  does not need to *calculate* the coloring; that has already been done.
    NTPartition featured_coloring() const;
    */

    // size_t num_nodes() const; -- defined in graph.cpp
    // size_t num_edges() const; -- defined in graph.cpp

    // returns the number of self-loops
    // size_t num_loops() const; -- defined in sparse_graph.cpp

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

    // Returns the internal label of the node used to represent edge (a, b).
    //
    // NOTE: This label is subject to change if the graph is modified.
    virtual int edge_node(const int a, const int b) const;

#ifndef SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
protected:
#endif

    NTSparseGraph& copy_assignment(const Graph& g);

    // size_t n; -- defined in graph.h
    size_t internal_n;
    // size_t m; -- defined in graph.h
    // size_t num_self_loops; -- defined in graph.h
    size_t num_edge_nodes;

    //  defined in sparse_graph.h
    // std::vector<std::unordered_set<int>> _neighbors;
    //  defined in sparse_graph.h
    //  Unused if graph is undirected
    // std::vector<std::unordered_set<int>> _out_neighbors;
    //  defined in sparse_graph.h
    //  Unused if graph is undirected
    // std::vector<std::unordered_set<int>> _in_neighbors;

    ////////////////////////////////////////////////////////////////////////////
    // The following functions and data structures manage the Nauty/Traces
    //  graph representation.
    ////////////////////////////////////////////////////////////////////////////

    // Called within delete_edge() and delete_node()
    //  Handles all the Nauty/Traces edge deletion
    void delete_edge_node_or_nodes(const int a, const int b);

    // Stores the ID of a node corresponding to edge (a, b).
    //  In a directed graph, each (undirected) edge really has two nodes. To
    //  access the second node, access the node for edge (b, a).
    //  In a directed graph, the edge node associated with edge (a, b) is the
    //  edge node that connects to a and to the other edge node.
    std::unordered_map<Edge, int, EdgeHash> edge_to_edge_node;
    std::unordered_map<int, Edge> edge_node_to_edge;

    // We require that the nodes of the actual graph get the first n labels.
    //  Then the edge nodes get the remaining labels.
    std::vector<int> out_degrees;
    std::vector<size_t> node_to_startpoint;
    std::vector<size_t> node_to_endpoint;

    std::vector<bool> has_self_loop;

    // Gives the two places in out_neighbors_vec where the internal node is
    //  listed as a neighbor.
    // When the graph is directed, we follow the
    //  convention that the first place listed is for the real node and the
    //  second place is for the other edge node.
    // When the graph is undirected, we follow the convention that the first
    //  place listed corresponds to the smaller of the two real node IDs.
    std::unordered_map<int, std::pair<size_t,size_t>> edge_node_to_places;

    // Used for regular and internal nodes.
    std::unordered_map<size_t, int> endpoint_to_node;


    // Only real nodes have extra space. Edge nodes always have a fixed amount.

    // A node is defined to have extra space when the following are all met:
    //  * (endpoint - startpoint) >= 2 * MIN_EDGE_SPACE_PER_NODE
    //  * (out_degree * 4) <= (endpoint - startpoint)
    //
    // NOTE: requiring that out_degree <= (endpoint - startpoint) / 4 IS NOT THE
    //  SAME AS requiring that (out_degree * 4) <= (endpoint - startpoint) due
    //  to integer math. We choose the latter for computational efficiency.
    //
    // We always give FLOOR(half the space) to the node that needs new space.
    //  Further, we require that the new space have at least twice as much space
    //  as the node receiving that space needs.

    // Call lower_bound(target) to get the node with the smallest capacity that
    //  is at least as large as your target storage space.

    AugmentedMultimap<size_t, int> extra_space_and_node;

    // When a node is an edge node and the graph is directed, we follow the
    //  convention that the first listed node is the real node and the second
    //  is the other edge node.
    // When a node is an edge node and the graph is undirected, we follow the
    //  convention that the first node listed has the smaller of the two node
    //  IDs.
    std::vector<int> out_neighbors_vec;

    // Must be an even number.
    const size_t MIN_EDGE_SPACE_PER_NODE = 4;

    int allocate_edge_node();
    void relabel_edge_node(const int a, const int b);

    // Moves an edge node's out_neighbors_vec info to the back of
    //  out_neighbors_vec, extending out_neighbors_vec to create the needed
    //  space. Does not relabel or delete any edge nodes.
    void slide_first_edge_node_to_back();

    // Can be used for new nodes as well that currently have no space and no
    //  edges. (i.e. endpoint = startpoint)
    // Only for regular nodes.
    void move_node_to_more_space(const int a);

    void move_edge_node_reference(const size_t init_loc,
                                  const size_t target_loc);

    // Updates out_degrees of main_node and the out_neighbors_vec
    void remove_edge_node_ref(const int main_node, const size_t ref_loc);

    // Removes reference to edge_node_of_slot and puts the last edge_node of
    //  out_neighbors_vec its place. If the last edge node of out_neighbors_vec
    //  IS edge_node_of_slot, then this removes edge_node_of_slot.
    //
    // ADDITIONALLY, this function must change the ID of the edge node with the
    //  highest ID to now have the ID of the destroyed node.
    void slide_back_edge_node_to_slot(int edge_node_of_slot);


    /////////////////////// Code/Vars for Colorings ////////////////////////

    /*
    bool structure_coloring;
    bool change_highlights_coloring;

    std::vector<int> node_ids;
    std::vector<int> partition_ints;
    std::vector<int> node_to_partition_loc;
    std::vector<int> node_to_cell;
    std::unordered_map<int, int> cell_starts;
    std::unordered_map<int, int> cell_ends;
    */
};

#endif
