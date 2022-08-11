/* Provides a C++ graph class that interfaces with the Nauty & Traces sparse
 *  graph classes.
 *
 * This class allows the use of edge types, as well as running Traces on
 *  directed graphs via node augmentation.
 *
 * This class prioritizes efficient asymptotic (i.e. big-O) runtime at the
 *  expense of data storage and constant factors.
 */

#include "nauty.h"
#include "nausparse.h"

#include "nt_mode.h"

#include<pair>
#include<unordered_map>
#include<unordered_set>
#include<vector>

#ifndef NT_GRAPH_H
#define NT_GRAPH_H

// TODO: Decide whether or not n must be constant.

const size_t NAUTY_TRACES_MAXN = 2000000000;

class SparseGraph {

public:
    // n is the number of nodes
    SparseGraph(const size_t n, const bool directed, const bool has_edge_types
                const NTMode &nt_mode);
    SparseGraph(const SparseGraph &g);
    SparseGraph(const SparseGraph &g, const NTMode &nt_mode);

    // returns a `sparsegraph` struct that can be passed into nauty or traces
    const sparsegraph as_nauty_traces_graph() const;

    size_t num_nodes() const;
    size_t num_edges() const;

    const bool directed;
    const bool has_edge_types;
    const NTMode &nt_mode;

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
    size_t n;
    size_t internal_n;
    size_t m;
    size_t m_undirected;

    std::unordered_set<int> 
    // Unused if graph is undirected
    std::unordered_set<int> out_neighbors_set;
    // Unused if graph is undirected
    std::unordered_set<int> in_neighbors_set;
    // Unused if graph does not have edge types.
    //  Maps an edge (a, b) to ((n1, n2), type) where n1 and n2 are the node(s)
    //      representing the edge and type is the edge type.
    std::unordered_map<std::pair<int, int>,
                       std::pair<std::pair<int, int>, int>> edge_types;

    std::vector<int> out_degrees;
    std::vector<int> node_vec_placement_maps;
    std::vector<int> out_neighbors_vec;
    std::vector<std::pair<size_t, size_t>> available_chunks;
};

#endif
