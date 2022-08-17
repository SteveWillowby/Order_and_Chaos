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
    virtual int delete_node(const int a) = 0;

    virtual void add_edge(const int a, const int b) = 0;
    virtual void delete_edge(const int a, const int b) = 0;

    // Deletes the edge if it exists and adds it if it does not.
    virtual void flip_edge(const int a, const int b) = 0;

    virtual const std::unordered_set<int> &neighbors(const int a) const = 0;
    // neighbors that node a points to
    virtual const std::unordered_set<int> &out_neighbors(const int a) const = 0;
    // neighbors that point to node a
    virtual const std::unordered_set<int> &in_neighbors(const int a) const = 0;

protected:
    size_t n;
    size_t m;

};

#endif
