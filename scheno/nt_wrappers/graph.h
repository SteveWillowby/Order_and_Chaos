/* Abstract graph class interface.
 *
 * Handles graphs with nodes labeled 0 through n-1.
 *
 * Supports directed and undirected graphs.
 * Supports self-loops.
 */

#include<cstddef>
#include<unordered_set>

#ifndef SCHENO__GRAPH_H
#define SCHENO__GRAPH_H

class Graph {

public:
    Graph(const bool directed);

    size_t num_nodes() const;
    size_t num_edges() const;
    size_t num_loops() const;  // returns the number of self-loops

    const bool directed;

    // Returns the id of the new node.
    virtual int add_node() = 0;
    // Deletes node a AND if a < n-1, relabels node n-1 to have label a.
    //  Returns the old label of the node that is now labeled a.
    virtual int delete_node(const int a) = 0;

    // Returns true iff the edge was absent (and thus now added)
    virtual bool add_edge(const int a, const int b) = 0;
    // Returns true iff the edge was present (and thus now deleted)
    virtual bool delete_edge(const int a, const int b) = 0;

    // Deletes the edge if it exists and adds it if it does not.
    virtual void flip_edge(const int a, const int b) = 0;

    virtual bool has_edge(const int a, const int b) const = 0;

    virtual const std::unordered_set<int> &neighbors(const int a) const = 0;
    // neighbors that node a points to
    virtual const std::unordered_set<int> &out_neighbors(const int a) const = 0;
    // neighbors that point to node a
    virtual const std::unordered_set<int> &in_neighbors(const int a) const = 0;

protected:
    size_t n;
    size_t m;
    size_t num_self_loops;

};

#endif
