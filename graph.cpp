#include "graph.h"

Graph::Graph(const bool directed) : directed(directed) {}

size_t Graph::num_nodes() const {
    return n;
}

size_t Graph::num_edges() const {
    return m;
}

size_t Graph::num_loops() const {
    return num_self_loops;
}
