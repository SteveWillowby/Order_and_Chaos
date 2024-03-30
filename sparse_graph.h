/* Implements the interface defined in graph.h
 *
 * Designed to be efficient with sparse graphs.
 *
 * Prioritizes asymptotic runtime over space efficiency.
 */

#include "graph.h"

#include<cstddef>
#include<unordered_set>
#include<vector>

#ifndef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
// Comment out this next line if you want this class's input error checks to
//  be ignored by the compiler.
#define SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
#endif

#ifndef SCHENO__SPARSE_GRAPH_H
#define SCHENO__SPARSE_GRAPH_H

class SparseGraph : public Graph {

public:
    SparseGraph(const bool directed);
    SparseGraph(const bool directed, size_t n);
    SparseGraph(const Graph &g);
    SparseGraph(const SparseGraph &g);

    SparseGraph& operator=(const Graph& g);
    SparseGraph& operator=(const SparseGraph& g);

    // size_t num_nodes() const; -- defined in graph.cpp
    // size_t num_edges() const; -- defined in graph.cpp

    // returns the number of self-loops
    // size_t num_loops() const; -- defined in graph.cpp

    // const bool directed; -- defined in graph.h

    // Returns the id of the new node.
    // O(1)
    virtual int add_node();
    // Deletes node a AND if a < n-1, relabels node n-1 to have label a.
    //  Returns the old label of the node that is now labeled a.
    // O(number of node a's neighbors + number of node n-1's neighbors)
    virtual int delete_node(int a);

    // Returns true iff the edge was absent (and thus now added)
    virtual bool add_edge(const int a, const int b);
    // Returns true iff the edge was present (and thus now deleted)
    virtual bool delete_edge(const int a, const int b);

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
    // size_t num_self_loops; -- defined in graph.h

    std::vector<std::unordered_set<int>> _neighbors;
    std::vector<std::unordered_set<int>> _out_neighbors;
    std::vector<std::unordered_set<int>> _in_neighbors;

    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
private:
    inline void range_check(const int a) const;
    #endif
};


#endif
