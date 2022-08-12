/* Abstract graph class interface.
 *
 * Designed to handle graphs with nodes labeled 0 through n-1.
 *
 * Supports directed and undirected graphs.
 */

#ifndef SYM__GRAPH_H
#define SYM__GRAPH_H

class Graph {

public:
    Graph();

    size_t num_nodes() const;
    size_t num_edges() const;

    const bool directed;

    // Returns the id of the new node.
    virtual int add_node() = 0;
    // Deletes node a AND if a < n-1, relabels node n-1 to have label a.
    virtual int delete_node(int a) = 0;

    // If the graph has edge types, sets the edge type to 0.
    virtual void add_edge(const int a, const int b) = 0;
    // Only call on graphs with edge types.
    virtual void add_edge(const int a, const int b, const int type) = 0;
    // Returns the edge type if there is one and zero otherwise.
    //  (note that 0 is an edge type, too)
    virtual int del_edge(const int a, const int b) = 0;

    // Deletes the edge if it exists and adds it if it does not.
    //  Returns the type of the changed edge if it has one and zero otherwise.
    //  Defaults to adding edge types of zero if the graph has edge types.
    virtual int flip_edge(const int a, const int b) = 0;
    virtual int flip_edge(const int a, const int b, const int type) = 0;

    // Requires that the graph has edge types and that edge (a, b) is present.
    virtual int set_edge_type(const int a, const int b, const int type) = 0;

    virtual const std::unordered_set<int> &neighbors(int a) const = 0;
    // neighbors that node a points to
    virtual const std::unordered_set<int> &out_neighbors(int a) const = 0;
    // neighbors that point to node a
    virtual const std::unordered_set<int> &in_neighbors(int a) const = 0;

protected:
    size_t n;
    size_t m;

};

#endif
