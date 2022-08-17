#include "edge.h"
#include "sparse_graph.h"

#include<stdexcept>
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
    if (n == 0) {
        throw std::domain_error("Error! Cannot make a graph with 0 nodes.");
    }
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

int SparseGraph::delete_node(const int a) {
    
}
