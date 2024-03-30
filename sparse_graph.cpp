#include "sparse_graph.h"

#include<stdexcept>
#include<string>
#include<unordered_set>
#include<vector>


SparseGraph::SparseGraph(const bool directed) : Graph(directed) {
    n = 1;
    m = 0;
    num_self_loops = 0;

    _neighbors =
        std::vector<std::unordered_set<int>>(n, std::unordered_set<int>());

    if (directed) {
        _out_neighbors =
            std::vector<std::unordered_set<int>>(n, std::unordered_set<int>());
        _in_neighbors =
            std::vector<std::unordered_set<int>>(n, std::unordered_set<int>());
    }
}

SparseGraph::SparseGraph(const bool directed, size_t n) : Graph(directed) {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    if (n == 0) {
        throw std::domain_error("Error! Cannot make a graph with 0 nodes.");
    }
    #endif
    this->n = n;
    m = 0;
    num_self_loops = 0;

    _neighbors =
        std::vector<std::unordered_set<int>>(n, std::unordered_set<int>());

    if (directed) {
        _out_neighbors =
            std::vector<std::unordered_set<int>>(n, std::unordered_set<int>());
        _in_neighbors =
            std::vector<std::unordered_set<int>>(n, std::unordered_set<int>());
    }
}

SparseGraph::SparseGraph(const Graph &g) : Graph(g.directed) {
    operator=(g);
}

SparseGraph::SparseGraph(const SparseGraph &g) : Graph(g.directed) {
    operator=(g);
}

SparseGraph& SparseGraph::operator=(const Graph& g) {
    if (g.directed != directed) {
        throw std::logic_error(
                std::string("Error! Cannot perform copy assignment when one") +
                " graph is directed and the other is not.");
    }

    if (this == &g) {
        return *this;
    }

    n = g.num_nodes();
    m = g.num_edges();
    num_self_loops = g.num_loops();

    _neighbors = std::vector<std::unordered_set<int>>(n);
    for (size_t i = 0; i < n; i++) {
        _neighbors[i] = std::unordered_set<int>(g.neighbors(i));
    }

    if (directed) {
        _out_neighbors = std::vector<std::unordered_set<int>>(n);
        _in_neighbors = std::vector<std::unordered_set<int>>(n);
        for (size_t i = 0; i < n; i++) {
            _out_neighbors[i] = std::unordered_set<int>(g.out_neighbors(i));
            _in_neighbors[i] = std::unordered_set<int>(g.in_neighbors(i));
        }
    }

    return *this;
}

SparseGraph& SparseGraph::operator=(const SparseGraph& g) {
    if (g.directed != directed) {
        throw std::logic_error(
                std::string("Error! Cannot perform copy assignment when one") +
                " graph is directed and the other is not.");
    }

    if (this == &g) {
        return *this;
    }

    n = g.num_nodes();
    m = g.num_edges();
    num_self_loops = g.num_loops();

    _neighbors = std::vector<std::unordered_set<int>>(n);
    for (size_t i = 0; i < n; i++) {
        _neighbors[i] = std::unordered_set<int>(g.neighbors(i));
    }

    if (directed) {
        _out_neighbors = std::vector<std::unordered_set<int>>(n);
        _in_neighbors = std::vector<std::unordered_set<int>>(n);
        for (size_t i = 0; i < n; i++) {
            _out_neighbors[i] = std::unordered_set<int>(g.out_neighbors(i));
            _in_neighbors[i] = std::unordered_set<int>(g.in_neighbors(i));
        }
    }

    return *this;
}

int SparseGraph::add_node() {
    n++;
    if (size_t(int(n - 1)) != n - 1) {
        throw std::out_of_range(
            "Error! Too many nodes to label with datatype int: " +
            std::to_string(n));
    }

    _neighbors.push_back(std::unordered_set<int>());
    if (directed) {
        _out_neighbors.push_back(std::unordered_set<int>());
        _in_neighbors.push_back(std::unordered_set<int>());
    }

    return n - 1;
}

int SparseGraph::delete_node(const int a) {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    n--;
    if (directed) {
        m -= (_out_neighbors[a].size() + _in_neighbors[a].size());
        if (_neighbors[a].find(a) != _neighbors[a].end()) {
            // One of the edges deleted was a self-loop. This would have shown
            //  up twice in the above "m -= ..." line.
            m++;
            num_self_loops--;
        }
    } else {
        m -= _neighbors[a].size();
        if (_neighbors[a].find(a) != _neighbors[a].end()) {
            num_self_loops--;
        }
    }

    int last_node = n;

    if (directed) {
        // If a points to a, this code will end up leaving a in
        //  _out_neighbors[a], but that is not a problem because _out_neighbors[a]
        //  is about to be deleted anyways.
        for (auto i = _out_neighbors[a].begin();
                i != _out_neighbors[a].end(); i++) {
            _in_neighbors[*i].erase(a);
        }
        for (auto i = _in_neighbors[a].begin();
                i != _in_neighbors[a].end(); i++) {
            _out_neighbors[*i].erase(a);
        }
    }
    for (auto i = _neighbors[a].begin();
            i != _neighbors[a].end(); i++) {
        // This check prevents a self-loop from causing iterator problems.
        //  _neighbors[a] will be deleted anyways, so it does not matter if its
        //  contents are fully erased or not.
        if (*i != a) {
            _neighbors[*i].erase(a);
        }
    }

    // At this point no nodes other than a point to a.
    //
    // We can now replace a's slot with last_node's info.
    //  This requires taking all nodes that pointed to last_node and making them
    //  point to last_node's new label (namely, a).
    if (directed && (last_node != a)) {
        _out_neighbors[a] = _out_neighbors[last_node];
        _in_neighbors[a] = _in_neighbors[last_node];

        for (auto i = _out_neighbors[a].begin();
                i != _out_neighbors[a].end(); i++) {
            if (*i == last_node) {
                continue;
            }
            _in_neighbors[*i].erase(last_node);
            _in_neighbors[*i].insert(a);
        }
        for (auto i = _in_neighbors[a].begin();
                i != _in_neighbors[a].end(); i++) {
            if (*i == last_node) {
                continue;
            }
            _out_neighbors[*i].erase(last_node);
            _out_neighbors[*i].insert(a);
        }

        // Check if there was a self-loop. If so, relabel it correctly.
        if (_out_neighbors[a].erase(last_node)) {
            _out_neighbors[a].insert(a);
            _in_neighbors[a].erase(last_node);
            _in_neighbors[a].insert(a);
        }
    }

    if (last_node != a) {
        _neighbors[a] = _neighbors[last_node];

        for (auto i = _neighbors[a].begin();
                i != _neighbors[a].end(); i++) {
            if (*i == last_node) {
                continue;
            }
            _neighbors[*i].erase(last_node);
            _neighbors[*i].insert(a);
        }

        // Check if there was a self-loop. If so, relabel it correctly.
        if (_neighbors[a].erase(last_node)) {
            _neighbors[a].insert(a);
        }
    }

    if (directed) {
        _out_neighbors.pop_back();
        _in_neighbors.pop_back();
    }
    _neighbors.pop_back();

    return last_node;
}


bool SparseGraph::add_edge(const int a, const int b) {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    range_check(b);
    #endif

    // insert() returns a pair, the second element of which is true iff
    //  the element is new
    if (directed) {
        if (_out_neighbors[a].insert(b).second) {
            _in_neighbors[b].insert(a);
            m++;

            if (_neighbors[a].insert(b).second) {
                _neighbors[b].insert(a);
            }

            if (a == b) {
                num_self_loops++;
            }

            return true;
        }
    } else {
        if (_neighbors[a].insert(b).second) {
            _neighbors[b].insert(a);
            m++;
            if (a == b) {
                num_self_loops++;
            }
            return true;
        }
    }
    return false;
}

bool SparseGraph::delete_edge(const int a, const int b) {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    range_check(b);
    #endif

    if (directed) {
        if (_out_neighbors[a].erase(b)) {
            _in_neighbors[b].erase(a);
            m--;
            if (_out_neighbors[b].find(a) == _out_neighbors[b].end()) {
                _neighbors[a].erase(b);
                _neighbors[b].erase(a);
            }
            if (a == b) {
                num_self_loops--;
            }
            return true;
        }
    } else {
        if (_neighbors[a].erase(b)) {
            _neighbors[b].erase(a);
            m--;
            if (a == b) {
                num_self_loops--;
            }
            return true;
        }
    }
    return false;
}

void SparseGraph::flip_edge(const int a, const int b) {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    range_check(b);
    #endif

    if (directed) {
        if (_out_neighbors[a].find(b) == _out_neighbors[a].end()) {
            add_edge(a, b);
        } else {
            delete_edge(a, b);
        }
    } else {
        if (_neighbors[a].find(b) == _neighbors[a].end()) {
            add_edge(a, b);
        } else {
            delete_edge(a, b);
        }
    }
}

bool SparseGraph::has_edge(const int a, const int b) const {
    if (directed) {
        return _out_neighbors[a].find(b) != _out_neighbors[a].end();
    }
    return _neighbors[a].find(b) != _neighbors[a].end();
}

const std::unordered_set<int> &SparseGraph::neighbors(const int a) const {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    return _neighbors[a];
}

const std::unordered_set<int> &SparseGraph::out_neighbors(const int a) const {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    if (directed) {
        return _out_neighbors[a];
    }
    return _neighbors[a];
}

const std::unordered_set<int> &SparseGraph::in_neighbors(const int a) const {
    #ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    range_check(a);
    #endif

    if (directed) {
        return _in_neighbors[a];
    }
    return _neighbors[a];
}

#ifdef SCHENO__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
void SparseGraph::range_check(const int a) const {
    if (a < 0 || size_t(a) >= n) {
        throw std::out_of_range("Error! Node " + std::to_string(a) +
                                " out of range - does not exist.");
    }
}
#endif
