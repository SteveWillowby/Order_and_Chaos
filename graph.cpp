#include "graph.h"

Graph::Graph() : directed(false) {
    n = 1;
    m = 0;
}

size_t Graph::num_nodes() const {
    return n;
}

size_t Graph::num_edges() const {
    return m;
}
