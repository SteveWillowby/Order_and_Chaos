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

#include "nauty.h"
#include "nausparse.h"

#include "coloring.h"
#include "edge.h"
#include "graph.h"

#include<pair>
#include<unordered_map>
#include<unordered_set>
#include<vector>

#ifndef SYM__NT_SPARSE_GRAPH_H
#define SYM__NT_SPARSE_GRAPH_H

const size_t NAUTY_TRACES_MAXN = 2000000000;

class NTSparseGraph : public Graph {

public:
    NTSparseGraph(const bool directed);
    // n is the number of nodes
    NTSparseGraph(const bool directed, const size_t n);
    NTSparseGraph(const Graph &g);

    // returns a `sparsegraph` struct that can be passed into nauty or traces
    const sparsegraph as_nauty_traces_graph() const;

    // Returns a coloring for the Nauty/Traces sparsegraph
    const Coloring<int> &
        nauty_traces_coloring(const Coloring<int> &node_coloring) const;
    const Coloring<int> &
        nauty_traces_coloring(const Coloring<Edge> &edge_coloring) const;
    const Coloring<int> &
        nauty_traces_coloring(const Coloring<int> &node_coloring,
                              const Coloring<Edge> &edge_coloring) const;

    // size_t num_nodes() const; -- defined in graph.cpp
    // size_t num_edges() const; -- defined in graph.cpp

    // const bool directed; -- defined in graph.h

    // Returns the id of the new node.
    // O(1)
    int add_node();
    // Deletes node a AND if a < n-1, relabels node n-1 to have label a.
    //  Returns the old label of the node that is now labeled a.
    // O(number of node a's neighbors + number of node n-1's neighbors)
    int delete_node(const int a);

    //////////////////////////////////////////////////////////////
    // All edge-changing functions have amortized O(1) runtime. //
    //////////////////////////////////////////////////////////////

    void add_edge(const int a, const int b);
    void delete_edge(const int a, const int b);

    // Deletes the edge if it exists and adds it if it does not.
    void flip_edge(const int a, const int b);

    bool has_edge(const int a, const int b) const;

    // O(1)
    const std::unordered_set<int> &neighbors(const int a) const;
    // neighbors that node a points to
    // O(1)
    const std::unordered_set<int> &out_neighbors(const int a) const;
    // neighbors that point to node a
    // O(1)
    const std::unordered_set<int> &in_neighbors(const int a) const;

private:
    // size_t n; -- defined in graph.h
    size_t internal_n;
    // size_t m; -- defined in graph.h
    size_t m_undirected;

    std::vector<std::unordered_set<int>> neighbors;
    // Unused if graph is undirected
    std::vector<std::unordered_set<int>> out_neighbors;
    // Unused if graph is undirected
    std::vector<std::unordered_set<int>> in_neighbors;

    std::vector<int> out_degrees;
    std::vector<int> node_vec_placement_maps;
    std::vector<int> out_neighbors_vec;
    std::vector<std::pair<size_t, size_t>> available_chunks;
};

#endif
