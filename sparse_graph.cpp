#include "edge.h"
#include "sparse_graph.h"

#include<stdexcept>
#include<string>
#include<unordered_map>
#include<vector>


SparseGraph::SparseGraph(const bool directed) : directed(directed) {
    n = 1;
    m = 0;

    neighbors =
        std::vector<std::unordered_map<int>>(n, std::unordered_map<int>());

    if (directed) {
        out_neighbors =
            std::vector<std::unordered_map<int>>(n, std::unordered_map<int>());
        in_neighbors =
            std::vector<std::unordered_map<int>>(n, std::unordered_map<int>());
    }
}

SparseGraph::SparseGraph(const bool directed, size_t n) : directed(directed) {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    if (n == 0) {
        throw std::domain_error("Error! Cannot make a graph with 0 nodes.");
    }
    #endif
    this->n = n;
    m = 0;

    neighbors =
        std::vector<std::unordered_map<int>>(n, std::unordered_map<int>());

    if (directed) {
        out_neighbors =
            std::vector<std::unordered_map<int>>(n, std::unordered_map<int>());
        in_neighbors =
            std::vector<std::unordered_map<int>>(n, std::unordered_map<int>());
    }
}

SparseGraph::SparseGraph(const Graph &g) : directed(g.directed) {
    n = g.num_nodes();
    m = g.num_edges();

    neighbors = std::vector<std::unordered_map<int>>(n);
    for (size_t i = 0; i < n; i++) {
        neighbors[i] = std::unordered_map<int>(g.neighbors(i));
    }

    if (directed) {
        out_neighbors = std::vector<std::unordered_map<int>>(n);
        in_neighbors = std::vector<std::unordered_map<int>>(n);
        for (size_t i = 0; i < n; i++) {
            out_neighbors[i] = std::unordered_map<int>(g.out_neighbors(i));
            in_neighbors[i] = std::unordered_map<int>(g.in_neighbors(i));
        }
    }
}


int SparseGraph::add_node() {
    n++;
    neighbors.push_back(std::unordered_map<int>());
    if (directed) {
        out_neighbors.push_back(std::unordered_map<int>());
        in_neighbors.push_back(std::unordered_map<int>());
    }
    return n - 1;
}

void SparseGraph::delete_node(const int a) {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    // TODO: Implement
}


void SparseGraph::add_edge(const int a, const int b) {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    range_check(b);
    #endif

    // insert() returns a pair, the second element of which is true iff
    //  the element is new
    if (directed) {
        if (out_neighbors[a].insert(b).second) {
            in_neighbors[b].insert(a);
            m++;

            if (neighbors[a].insert(b).second) {
                neighbors[b].insert(a);
            }
        }
    } else {
        if (neighbors[a].insert(b).second) {
            neighbors[b].insert(a);
            m++;
        }
    }
}

void SparseGraph::delete_edge(const int a, const int b) {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    range_check(b);
    #endif

    // TODO: Implement
}

void SparseGraph::flip_edge(const int a, const int b) {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    range_check(b);
    #endif

    // TODO: Implement
}

const std::unordered_set<int> &SparseGraph::neighbors(const int a) const {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    return neighbors[a];
}

const std::unordered_set<int> &SparseGraph::out_neighbors(const int a) const {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    return out_neighbors[a];
}

const std::unordered_set<int> &SparseGraph::in_neighbors(const int a) const {
    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    return in_neighbors[a];
}

void SparseGraph::range_check(const int a) const {
    if (a < 0 || a >= n) {
        throw std::out_of_range("Error! Node " + std::to_string(a) +
                                " out of range - does not exist.");
    }
}
