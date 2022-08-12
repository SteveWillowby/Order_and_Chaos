/* Implements the interface defined in graph.h
 *
 * Designed to be efficient with sparse graphs.
 *
 * Prioritizes asymptotic runtime over space efficiency.
 */

#include "graph.cpp"

#ifndef SYM__SPARSE_GRAPH_H
#define SYM__SPARSE_GRAPH_H


class SparseGraph : Graph {

public:
    SparseGraph(const bool directed, const bool has_edge_types);
    SparseGraph(const bool directed, const bool has_edge_types, size_t n);
    SparseGraph(const Graph &g);

    // size_t num_nodes() const; -- defined in graph.cpp
    // size_t num_edges() const; -- defined in graph.cpp

    const bool directed;

    // Returns the id of the new node.
    // O(1)
    int add_node();
    // Deletes node a AND if a < n-1, relabels node n-1 to have label a.
    // O(number of node a's neighbors + number of node n-1's neighbors)
    int delete_node(int a);

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

    const std::unordered_set<int> &neighbors(int a) const;
    // neighbors that node a points to
    const std::unordered_set<int> &out_neighbors(int a) const;
    // neighbors that point to node a
    const std::unordered_set<int> &in_neighbors(int a) const;

protected:
    // size_t n; -- defined in Graph
    // size_t m; -- defined in Graph

};


#endif
