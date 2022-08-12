/* Provides a C++ graph class that interfaces with the Nauty & Traces sparse
 *  graph classes.
 *
 * Extends the SparseGraph class defined in sparse_graph.h
 *
 * This class prioritizes efficient asymptotic (i.e. big-O) runtime at the
 *  expense of data storage and constant factors.
 *
 * Undirected graphs can have 1 <= n <= 2,000,000,000 and any m.
 * Directed graphs are limited to 1 <= n + 2m <= 2,000,000,000.
 */

#include "nauty.h"
#include "nausparse.h"

#include "coloring.h"
#include "edge.h"
#include "sparse_graph.h"

#include<pair>
#include<unordered_map>
#include<unordered_set>
#include<vector>

#ifndef SYM__NT_SPARSE_GRAPH_H
#define SYM__NT_SPARSE_GRAPH_H

const size_t NAUTY_TRACES_MAXN = 2000000000;

class NTSparseGraph : SparseGraph {

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
    // O(number of node a's neighbors + number of node n-1's neighbors)
    int delete_node(int a);

    //////////////////////////////////////////////////////////////
    // All edge-changing functions have amortized O(1) runtime. //
    //////////////////////////////////////////////////////////////

    // If the graph has edge types, sets the edge type to 0.
    void add_edge(const int a, const int b);
    // Only call on graphs with edge types.
    void add_edge(const int a, const int b, const int type);
    // Returns the edge type if there is one and zero otherwise.
    //  (note that 0 is an edge type, too)
    int del_edge(const int a, const int b);

    // Deletes the edge if it exists and adds it if it does not.
    //  Returns the type of the changed edge if it has one and zero otherwise.
    //  Defaults to adding edge types of zero if the graph has edge types.
    int flip_edge(const int a, const int b);
    int flip_edge(const int a, const int b, const int type);

    // Requires that the graph has edge types and that edge (a, b) is present.
    int set_edge_type(const int a, const int b, const int type);

    // O(1)
    const std::unordered_set<int> &neighbors(int a) const;
    // neighbors that node a points to
    // O(1)
    const std::unordered_set<int> &out_neighbors(int a) const;
    // neighbors that point to node a
    // O(1)
    const std::unordered_set<int> &in_neighbors(int a) const;

private:
    // size_t n; -- defined in graph.h
    size_t internal_n;
    // size_t m; -- defined in graph.h
    size_t m_undirected;

    std::vector<std::unordered_set<int>> neighbors_sets;
    // Unused if graph is undirected
    std::vector<std::unordered_set<int>> out_neighbors_sets;
    // Unused if graph is undirected
    std::vector<std::unordered_set<int>> in_neighbors_sets;
    // Unused if graph does not have edge types.
    //  Maps an edge (a, b) to ((n1, n2), type) where n1 and n2 are the nodes
    //      representing the edge and type is the edge type.
    std::unordered_map<Edge, std::pair<std::pair<int, int>, int>> edge_types;

    std::vector<int> out_degrees;
    std::vector<int> node_vec_placement_maps;
    std::vector<int> out_neighbors_vec;
    std::vector<std::pair<size_t, size_t>> available_chunks;
};

#endif
