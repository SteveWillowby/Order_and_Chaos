#include "edge_sampler.h"
#include "nt_wrappers/nauty_traces.h"

#include<cstdint>
#include<random>
#include<stdexcept>
#include<string>


EdgeSampler::EdgeSampler(const Graph& g, std::mt19937& gen) :
                           directed(g.directed), n(g.num_nodes()),
                           has_self_loops(g.num_loops() > 0),
                           generator(gen) {

    dist = std::uniform_real_distribution<double>(0.0, 1.0);

    last_op = 0;  // none

    size_t num_edges = g.num_edges();
    size_t max_num_edges = (n * size_t(has_self_loops)) +
                           (n * (n - 1)) / (size_t(!directed) + 1);

    if (n * n > SCHENO__MAX_EDGE_LABEL) {
        throw std::range_error(std::string("Error! An EdgeSampler for a ") + 
                               "graph with " + std::to_string(n) + " nodes " +
                               "requires a larger edge-int type. "
                               "Redefine SCHENO__edge_int_type " +
                               "and SCHENO__MAX_EDGE_LABEL in edge.h.");
    }

    edges_un_sampled = std::vector<SCHENO__edge_int_type>(num_edges, 0);
    non_edges_un_sampled =
        std::vector<SCHENO__edge_int_type>(max_num_edges - num_edges, 0);
    edges_sampled = std::vector<SCHENO__edge_int_type>();
    non_edges_sampled = std::vector<SCHENO__edge_int_type>();

    size_t next_edge_idx = 0;
    size_t next_non_edge_idx = 0;
    SCHENO__edge_int_type edge_id;
    for (size_t a = 0; a < n; a++) {
        for (size_t b = size_t(!directed) * a; b < n; b++) {
            if (a == b && !has_self_loops) {
                b++;
                if (b == n) {
                    break;
                }
            }
            edge_id = a * n + b;
            if (g.has_edge(a, b)) {
                edges_un_sampled[next_edge_idx] = edge_id;
                next_edge_idx++;
            } else {
                non_edges_un_sampled[next_non_edge_idx] = edge_id;
                next_non_edge_idx++;
            }
        }
    }
}

Edge EdgeSampler::sample_edge() {
    if (edges_un_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::sample_edge()")
                               + " -- All edges have already been sampled. "
                               + "Un-sample some edges first.");
    }
    last_op = 1;

    size_t idx = dist(generator) * edges_un_sampled.size();
    SCHENO__edge_int_type e = edges_un_sampled[idx];
    edges_un_sampled[idx] = edges_un_sampled[edges_un_sampled.size() - 1];
    edges_un_sampled.pop_back();

    edges_sampled.push_back(e);
    return int_to_edge(e);
}

Edge EdgeSampler::sample_non_edge() {
    if (non_edges_un_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::sample_non_")
                               + "edge() -- All non-edges have already been "
                               + "sampled. Un-sample some non-edges first.");
    }
    last_op = 2;

    size_t idx = dist(generator) * non_edges_un_sampled.size();
    SCHENO__edge_int_type e = non_edges_un_sampled[idx];
    non_edges_un_sampled[idx] =
        non_edges_un_sampled[non_edges_un_sampled.size() - 1];
    non_edges_un_sampled.pop_back();

    non_edges_sampled.push_back(e);
    return int_to_edge(e);
    
}

Edge EdgeSampler::un_sample_edge() {
    if (edges_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::un_sample_edge")
                               + "() -- No edges have been sampled. "
                               + "Sample some edges first.");
    }
    last_op = 3;

    size_t idx = dist(generator) * edges_sampled.size();
    SCHENO__edge_int_type e = edges_sampled[idx];
    edges_sampled[idx] = edges_sampled[edges_sampled.size() - 1];
    edges_sampled.pop_back();

    edges_un_sampled.push_back(e);
    return int_to_edge(e);
}

Edge EdgeSampler::un_sample_non_edge() {
    if (non_edges_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::un_sample_non_")
                               + "edge() -- No non-edges have been sampled. "
                               + "Sample some non-edges first.");
    }
    last_op = 4;

    size_t idx = dist(generator) * non_edges_sampled.size();
    SCHENO__edge_int_type e = non_edges_sampled[idx];
    non_edges_sampled[idx] = non_edges_sampled[non_edges_sampled.size() - 1];
    non_edges_sampled.pop_back();

    non_edges_un_sampled.push_back(e);
    return int_to_edge(e);
}

std::pair<Edge, Edge> EdgeSampler::swap_edge_samples() {
    if (edges_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::swap_edge_")
                               + "samples() -- No edges are sampled. "
                               + "Sample an edge first.");
    } else if (edges_un_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::swap_edge_")
                               + "samples() -- All edges are sampled. "
                               + "Un-sample an edge first.");
    }
    last_op = 5;

    size_t idx_old = dist(generator) * edges_sampled.size();
    size_t idx_new = dist(generator) * edges_un_sampled.size();

    SCHENO__edge_int_type old_edge = edges_sampled[idx_old];
    SCHENO__edge_int_type new_edge = edges_un_sampled[idx_new];

    edges_sampled[idx_old] = edges_sampled[edges_sampled.size() - 1];
    edges_sampled[edges_sampled.size() - 1] = new_edge;
    edges_un_sampled[idx_new] = edges_un_sampled[edges_un_sampled.size() - 1];
    edges_un_sampled[edges_un_sampled.size() - 1] = old_edge;

    return std::pair<Edge, Edge>(int_to_edge(new_edge), int_to_edge(old_edge));
}

std::pair<Edge, Edge> EdgeSampler::swap_non_edge_samples() {
    if (non_edges_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::swap_non_edge_")
                               + "samples() -- No non-edges are sampled. "
                               + "Sample a non-edge first.");
    } else if (non_edges_un_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::swap_non_edge_")
                               + "samples() -- All non-edges are sampled. "
                               + "Un-sample a non-edge first.");
    }
    last_op = 6;

    size_t idx_old = dist(generator) * non_edges_sampled.size();
    size_t idx_new = dist(generator) * non_edges_un_sampled.size();

    SCHENO__edge_int_type old_edge = non_edges_sampled[idx_old];
    SCHENO__edge_int_type new_edge = non_edges_un_sampled[idx_new];

    non_edges_sampled[idx_old] = non_edges_sampled[non_edges_sampled.size()-1];
    non_edges_sampled[non_edges_sampled.size() - 1] = new_edge;
    non_edges_un_sampled[idx_new] =
        non_edges_un_sampled[non_edges_un_sampled.size() - 1];
    non_edges_un_sampled[non_edges_un_sampled.size() - 1] = old_edge;

    return std::pair<Edge, Edge>(int_to_edge(new_edge), int_to_edge(old_edge));
}

void EdgeSampler::undo() {
    if (last_op == 0) {
        throw std::logic_error(std::string("Error! Cannot undo sampling from ")
                               + "this EdgeSampler object. Either no operations"
                               + " have been done or you already un-did the "
                               + "latest operation.");
    } else if (last_op == 1) {
        edges_un_sampled.push_back(*(edges_sampled.rbegin()));
        edges_sampled.pop_back();
    } else if (last_op == 2) {
        non_edges_un_sampled.push_back(*(non_edges_sampled.rbegin()));
        non_edges_sampled.pop_back();
    } else if (last_op == 3) {
        edges_sampled.push_back(*(edges_un_sampled.rbegin()));
        edges_un_sampled.pop_back();
    } else if (last_op == 4) {
        non_edges_sampled.push_back(*(non_edges_un_sampled.rbegin()));
        non_edges_un_sampled.pop_back();
    } else if (last_op == 5) {
        SCHENO__edge_int_type temp = edges_sampled[edges_sampled.size() - 1];
        edges_sampled[edges_sampled.size() - 1] =
            edges_un_sampled[edges_un_sampled.size() - 1];
        edges_un_sampled[edges_un_sampled.size() - 1] = temp;
    } else {
        SCHENO__edge_int_type temp = non_edges_sampled[non_edges_sampled.size()-1];
        non_edges_sampled[non_edges_sampled.size() - 1] =
            non_edges_un_sampled[non_edges_un_sampled.size() - 1];
        non_edges_un_sampled[non_edges_un_sampled.size() - 1] = temp;
    }

    // Mark the undo as completed.
    last_op = 0;
}

Edge EdgeSampler::int_to_edge(SCHENO__edge_int_type e) const {
    return EDGE(e / n, e % n, directed);
}

/*
SCHENO__edge_int_type EdgeSampler::edge_to_int(const Edge& e) const {
    return (e.first * n) + e.second;
}
*/
