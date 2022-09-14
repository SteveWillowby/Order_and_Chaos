/* Implements the interface defined in graph.h
 *
 * Designed to be efficient with sparse graphs.
 *
 * Prioritizes asymptotic runtime over space efficiency.
 */

#include "graph.h"

#include<unordered_set>
#include<vector>

#ifndef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
// Comment out this next line if you want this class's input error checks to
//  be ignored by the compiler.
#define SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
#endif

#ifndef SYM__SPARSE_GRAPH_H
#define SYM__SPARSE_GRAPH_H

class SparseGraph : public Graph {

public:
    SparseGraph(const bool directed);
    SparseGraph(const bool directed, size_t n);
    SparseGraph(const Graph &g);

    // size_t num_nodes() const; -- defined in graph.cpp
    // size_t num_edges() const; -- defined in graph.cpp

    // const bool directed; -- defined in graph.h

    // Returns the id of the new node.
    // O(1)
    virtual int add_node();
    // Deletes node a AND if a < n-1, relabels node n-1 to have label a.
    //  Returns the old label of the node that is now labeled a.
    // O(number of node a's neighbors + number of node n-1's neighbors)
    virtual int delete_node(int a);

    virtual void add_edge(const int a, const int b);
    virtual void delete_edge(const int a, const int b);

    // Deletes the edge if it exists and adds it if it does not.
    virtual void flip_edge(const int a, const int b);

    virtual bool has_edge(const int a, const int b) const;

    virtual const std::unordered_set<int> &neighbors(const int a) const;
    // neighbors that node a points to
    virtual const std::unordered_set<int> &out_neighbors(const int a) const;
    // neighbors that point to node a
    virtual const std::unordered_set<int> &in_neighbors(const int a) const;

protected:
    // size_t n; -- defined in graph.h
    // size_t m; -- defined in graph.h

    std::vector<std::unordered_set<int>> _neighbors;
    std::vector<std::unordered_set<int>> _out_neighbors;
    std::vector<std::unordered_set<int>> _in_neighbors;

    inline void range_check(const int a) const;
};


#endif
