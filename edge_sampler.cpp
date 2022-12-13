#include "edge.h"
#include "edge_sampler.h"
#include "sparse_graph.h"

#include<cstdint>
#include<random>
#include<stdexcept>
#include<string>


EdgeSampler::EdgeSampler(const SparseGraph& g, std::mt19937& gen) :
                           generator(gen), directed(g.directed),
                           has_self_loops(g.num_loops() > 0), n(g.num_nodes()) {

    dist = std::uniform_real_distribution<double>(0.0, 1.0);

    last_op = 0;  // none

    size_t num_edges = g.num_edges();
    size_t max_num_edges = (n * size_t(has_self_loops)) +
                           (n * (n - 1)) / (size_t(!directed) + 1);

    if (n * n > SYM__MAX_EDGE_LABEL) {
        throw std::range_error(std::string("Error! An EdgeSampler for a ") + 
                               "graph with " + std::to_string(n) + " nodes " +
                               "requires a larger edge-int type. "
                               "Redefine SYM__edge_int_type " +
                               "and SYM__MAX_EDGE_LABEL in edge_sampler.h.");
    }

    edges_un_sampled = std::vector<SYM__edge_int_type>(num_edges, 0);
    non_edges_un_sampled =
        std::vector<SYM__edge_int_type>(max_num_edges - num_edges, 0);
    edges_sampled = std::vector<SYM__edge_int_type>();
    non_edges_sampled = std::vector<SYM__edge_int_type>();

    size_t next_edge_idx = 0;
    size_t next_non_edge_idx = 0;
    SYM__edge_int_type edge_id;
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
    SYM__edge_int_type e = edges_un_sampled[idx];
    edges_un_sampled[idx] = edges_un_sampled[edges_un_sampled.size() - 1];
    edges_un_sampled.pop_back();

    edges_sampled.push_back(e);
    return int_to_edge(e);
}

Edge sample_non_edge() {
    if (non_edges_un_sampled.size() == 0) {
        throw std::range_error(std::string("Error! EdgeSampler::sample_non_")
                               + "edge() -- All non-edges have already been "
                               + "sampled. Un-sample some non-edges first.");
    }
    last_op = 2;

    size_t idx = dist(generator) * non_edges_un_sampled.size();
    SYM__edge_int_type e = non_edges_un_sampled[idx];
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
    SYM__edge_int_type e = edges_sampled[idx];
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
    SYM__edge_int_type e = non_edges_sampled[idx];
    non_edges_sampled[idx] = non_edges_sampled[non_edges_sampled.size() - 1];
    non_edges_sampled.pop_back();

    non_edges_un_sampled.push_back(e);
    return int_to_edge(e);
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
    } else {  // last_op == 4
        non_edges_sampled.push_back(*(non_edges_un_sampled.rbegin()));
        non_edges_un_sampled.pop_back();
    }

    // Mark the undo as completed.
    last_op = 0;
}

Edge EdgeSampler::int_to_edge(SYM__edge_int_type e) const {
    return EDGE(e / n, e % n, directed);
}

/*
SYM__edge_int_type EdgeSampler::edge_to_int(const Edge& e) const {
    return (e.first * n) + e.second;
}
*/
