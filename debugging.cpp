#include<algorithm>
#include<iostream>
#include<random>
#include<string>

#include "nt_sparse_graph.h"

#ifdef SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
bool consistency_check(const NTSparseGraph &g) {
    // First, check the SparseGraph code for consistency.
    size_t num_self_loops = 0;
    size_t num_edges_found = 0;
    size_t num_nodes = g.num_nodes();
    if (g.directed) {
        for (int i = 0; i < int(num_nodes); i++) {
            num_edges_found += g.out_neighbors(i).size();
            if (g.neighbors(i).find(i) != g.neighbors(i).end()) {
                num_self_loops++;
            }
        }
        if (num_edges_found != g.num_edges()) {
            std::cout<<"SUM of out_neighbors() sizes != g.num_edges()"
                     <<num_edges_found<<" vs. "<<g.num_edges()<<std::endl;
            return false;
        }
    } else {
        for (int i = 0; i < int(num_nodes); i++) {
            num_edges_found += g.neighbors(i).size();
            if (g.neighbors(i).find(i) != g.neighbors(i).end()) {
                num_self_loops++;
            }
        }
        if (num_edges_found !=
                ((g.num_edges() - num_self_loops) * 2 + num_self_loops)) {
            std::cout<<"SUM of neighbors() sizes != g.num_edges()"<<std::endl;
            return false;
        }
    }

    if (num_self_loops != g.num_loops()) {
        std::cout<<"Num Self Loops Found != g.num_loops()"<<std::endl;
        return false;
    }

    if (g.directed) {
        for (int i = 0; i < int(num_nodes); i++) {
            for (auto itr = g.out_neighbors(i).begin();
                      itr != g.out_neighbors(i).end(); itr++) {
                if (g.neighbors(i).find(*itr) == g.neighbors(i).end()) {
                    std::cout<<*itr<<" is in "<<i<<"'s out_neighbors but not in"
                             <<i<<"'s neighbors"<<std::endl;
                    return false;
                }
                if (g.in_neighbors(*itr).find(i) == g.in_neighbors(*itr).end()){
                    std::cout<<*itr<<" is in "<<i<<"'s out_neighbors but "<<i
                             <<"is not in "<<*itr<<"'s in_neighbors"<<std::endl;
                    return false;
                }
            }
            for (auto itr = g.in_neighbors(i).begin();
                      itr != g.in_neighbors(i).end(); itr++) {
                if (g.neighbors(i).find(*itr) == g.neighbors(i).end()) {
                    std::cout<<*itr<<" is in "<<i<<"'s in_neighbors but not in"
                             <<i<<"'s neighbors"<<std::endl;
                    return false;
                }
                if (g.out_neighbors(*itr).find(i) ==
                            g.out_neighbors(*itr).end()) {
                    std::cout<<*itr<<" is in "<<i<<"'s in_neighbors but "<<i
                             <<"isn't in "<<*itr<<"'s out_neighbors"<<std::endl;
                    return false;
                }
            }
            for (auto itr = g.neighbors(i).begin();
                      itr != g.neighbors(i).end(); itr++) {
                if (g.out_neighbors(i).find(*itr) == g.out_neighbors(i).end() &&
                      g.in_neighbors(i).find(*itr) == g.in_neighbors(i).end()) {
                    std::cout<<*itr<<" is in "<<i<<"'s neighbors but not in "<<i
                             <<"'s in_ or out_neighbors."<<std::endl;
                    return false;
                }
            }
        }
    } else {
        for (int i = 0; i < int(num_nodes); i++) {
            for (auto itr = g.neighbors(i).begin();
                      itr != g.neighbors(i).end(); itr++) {
                if (g.neighbors(*itr).find(i) == g.neighbors(*itr).end()) {
                    std::cout<<*itr<<" is in "<<i<<"'s neighbors but "<<i<<"is "
                             <<"not in "<<*itr<<"'s neighbors"<<std::endl;
                    return false;
                }
            }
        }
    }

    // Now check the Nauty/Traces data.

    if (g.out_degrees.size() != g.internal_n) {
        std::cout<<"g.out_degrees.size() != g.internal_n"<<std::endl;
        return false;
    }
    if (g.num_nodes() + g.num_edge_nodes != g.internal_n) {
        std::cout<<"g.num_nodes() + g.num_edge_nodes != g.internal_n"<<std::endl;
        return false;
    }
    if (g.node_to_startpoint.size() != g.internal_n) {
        std::cout<<"g.node_to_startpoint.size() != g.internal_n"<<std::endl;
        return false;
    }
    if (g.node_to_endpoint.size() != g.internal_n) {
        std::cout<<"g.node_to_endpoint.size() != g.internal_n"<<std::endl;
        return false;
    }
    if (g.endpoint_to_node.size() != g.internal_n) {
        std::cout<<"g.endpoint_to_node.size() != g.internal_n"
                 <<" ("<<g.endpoint_to_node.size()<<" vs. "
                 <<g.internal_n<<")"<<std::endl;
        return false;
    }
    if (g.edge_node_to_edge.size() != g.num_edge_nodes) {
        std::cout<<"g.edge_node_to_edge.size() != g.num_edge_nodes"<<std::endl;
        return false;
    }
    if (g.edge_to_edge_node.size() != g.num_edge_nodes) {
        std::cout<<"g.edge_to_edge_node.size() != g.num_edge_nodes"<<std::endl;
        return false;
    }
    if (g.edge_node_to_places.size() != g.num_edge_nodes) {
        std::cout<<"g.edge_node_to_places.size() != g.num_edge_nodes"
                 <<" ("<<g.edge_node_to_places.size()<<" vs. "
                 <<g.num_edge_nodes<<")"<<std::endl;
        return false;
    }

    // Verify that has_self_loop is correct.
    for (size_t i = 0; i < g.num_nodes(); i++) {
        if (g.has_self_loop[i] !=
                (g.neighbors(i).find(i) != g.neighbors(i).end())) {
            std::cout<<"g.has_self_loop[i] != (g.neighbors(i).find(i)"
                     <<" != g.neighbors(i).end())"<<std::endl;
            return false;
        }
    }

    // Verify that out_degrees is correct. Account for self-loops not being
    //  referenced by g.out_degrees but in a different vector.
    for (size_t i = 0; i < g.num_nodes(); i++) {
        if (g.out_degrees[i] + int(g.has_self_loop[i]) !=
                int(g.neighbors(i).size())) {
            std::cout<<"g.out_degrees["<<i<<"] + int(g.has_self_loop["<<i
                     <<"]) != g.neighbors("<<i<<").size()"<<std::endl;
            return false;
        }
        if (g.out_degrees[i] + g.node_to_startpoint[i] > g.node_to_endpoint[i]){
            std::cout<<"g.out_degrees["<<i<<"] + g.node_to_startpoint["<<i
                     <<"] > g.node_to_endpoint["<<i<<"]"<<std::endl;
            return false;
        }
    }
    for (size_t i = g.num_nodes(); i < g.internal_n; i++) {
        if (g.out_degrees[i] != 2) {
            std::cout<<"g.out_degrees["<<i<<"] != 2, but rather equals "
                     <<g.out_degrees[i]<<std::endl;
            return false;
        }
    }

    // Verify that the space of space referenced in out_neighbors_vec is the
    //  same as the size of out_neighbors_vec
    size_t space_accounted_for = 0;
    for (size_t i = 0; i < g.node_to_startpoint.size(); i++) {
        space_accounted_for += g.node_to_endpoint[i] - g.node_to_startpoint[i];
        if (g.node_to_endpoint[i] > g.out_neighbors_vec.size()) {
            std::cout<<"g.node_to_endpoint["<<i<<"] > g.out_neighbors_vec.size()"
                     <<std::endl;
            return false;
        }
    }
    auto all_extra_space_pairs = g.extra_space_and_node.all_pairs();
    for (auto itr = all_extra_space_pairs.begin();
            itr < all_extra_space_pairs.end(); itr++) {
        if (itr->second == -1) {
            space_accounted_for += itr->first;
        }
    }
    if (space_accounted_for != g.out_neighbors_vec.size()) {
        std::cout<<"space_accounted_for != g.out_neighbors_vec.size() "
                 <<space_accounted_for<<" vs "<<g.out_neighbors_vec.size()
                 <<std::endl;
        return false;
    }

    // Verify that endpoint_to_node and node_to_endpoint are consistent.
    for (size_t i = 0; i < g.internal_n; i++) {
        int endpoint = g.node_to_endpoint[i];
        auto result = g.endpoint_to_node.find(endpoint);
        if (result == g.endpoint_to_node.end()) {
            std::cout<<"g.endpoint_to_node does not contain node "<<i
                     <<"'s endpoint ("<<endpoint<<")."<<std::endl;
            return false;
        }
        if (int(i) != result->second) {
            std::cout<<"g.endpoint_to_node[node_to_endpoint["<<i<<"]] != "<<i
                     <<std::endl;
            return false;
        }
    }

    // Verify that edge_node_to_edge and edge_to_edge_node are consistent.
    for (size_t i = g.num_nodes(); i < g.internal_n; i++) {
        auto edge_result = g.edge_node_to_edge.find(i);
        if (edge_result == g.edge_node_to_edge.end()) {
            std::cout<<"g.edge_node_to_edge missing edge node "<<i<<std::endl;
            return false;
        }
        auto edge = edge_result->second;
        auto node_result = g.edge_to_edge_node.find(edge);
        if (node_result == g.edge_to_edge_node.end()) {
            std::cout<<"g.edge_to_edge_node missing edge ("<<edge.first<<", "
                     <<edge.second<<")"<<std::endl;
            return false;
        }
        if (node_result->second != int(i)) {
            std::cout<<"g.edge_to_edge_node[g.edge_node_to_edge["<<i<<"]] != "
                     <<i<<std::endl;
            return false;
        }
    }

    // Check that every spot recorded as having extra space in fact has the
    //  amount recorded.
    std::vector<std::pair<const size_t&, const int&>> extra_space =
                                            g.extra_space_and_node.all_pairs();
    for (auto itr = extra_space.begin(); itr != extra_space.end(); itr++) {
        size_t size = itr->first;
        int node = itr->second;

        if (node == -1) {
            if (g.endpoint_to_node.find(size) != g.endpoint_to_node.end()) {
                std::cout<<"Extra space is recoreded as having size "<<size
                         <<" at start of vector (i.e. 'node -1') but there is "
                         <<"an endpoint recorded at point "<<size<<" recorded "
                         <<"as being for node "
                         <<g.endpoint_to_node.find(size)->second<<std::endl;
                return false;
            }
            continue;
        }

        int space = g.node_to_endpoint[node] - g.node_to_startpoint[node];
        if (!(space >= 2 * int(g.MIN_EDGE_SPACE_PER_NODE) &&
                space >= 4 * g.out_degrees[node])) {
            std::cout<<"Node "<<node<<" is recorded as having extra space "
                     <<"when in fact it does not."<<std::endl;
            return false;
        } else if (space / 2 != int(size)) {
            std::cout<<"Node "<<node<<" is recorded as having "<<size
                     <<" extra slots when in fact it has "<<space/2<<std::endl;
            return false;
        }
    }

    // Check that every node which has extra space is recorded as having it.
    for (int node = 0; node < int(g.num_nodes()); node++) {
        int space = g.node_to_endpoint[node] - g.node_to_startpoint[node];
        if (space >= 2 * int(g.MIN_EDGE_SPACE_PER_NODE) &&
                space >= 4 * g.out_degrees[node]) {
            int capacity = space / 2;
            if (!g.extra_space_and_node.contains(capacity, node)) {
                std::cout<<"Node "<<node<<" should be recorded as having "
                         <<capacity<<" extra slots but either it's not recorded"
                         <<" as having any slots or is recorded as having the "
                         <<"the wrong capacity."<<std::endl;
                return false;
            }
        }
    }
    std::vector<size_t> startpoints = std::vector<size_t>(g.node_to_startpoint);
    std::sort(startpoints.begin(), startpoints.end());
    if (startpoints[0] > 0) {
        bool found = false;
        for (auto itr = extra_space.begin(); itr != extra_space.end(); itr++) {
            size_t size = itr->first;
            int node = itr->second;

            if (node == -1) {
                found = true;
                if (size != startpoints[0]) {
                    std::cout<<"The extra space associated with node '-1' shoul"
                             <<"d be "<<startpoints[0]<<" but it is "<<size
                             <<std::endl;
                    return false;
                }
            }
        }
        if (!found) {
            std::cout<<"There should be an extra space with node '-1' recorded"
                     <<" but there is no such recording."<<std::endl;
            return false;
        }
    }

    // Verify that edge_node_to_places is correct.
    //
    // Not only does edge node e appear at the two places that
    //  edge_node_to_places[e] references, but also those two places are within
    //  the domain of an actual node (between startpoint & startpoint + degree).
    for (size_t i = 0; i < g.internal_n; i++) {
        for (size_t idx = g.node_to_startpoint[i];
              idx < g.node_to_startpoint[i] + size_t(g.out_degrees[i]); idx++) {
            int node = g.out_neighbors_vec[idx];
            if (node < 0 || (node < int(g.num_nodes()) && i < g.num_nodes())) {
                std::cout<<"g.out_neighbors_vec["<<idx<<"] < 0 or it is"
                         <<" < g.num_nodes() while "<<idx<<" is an idx for a regular node"
                         <<std::endl;
                return false;
            }

            if (node < int(g.num_nodes())) {
                continue;
            }

            // Node is an edge node.
            auto places_itr = g.edge_node_to_places.find(node);
            if (places_itr == g.edge_node_to_places.end()) {
                std::cout<<"Edge node "<<node<<" in out_neighbors_vec but not"
                         <<" in edge_node_to_places. One or the other is wrong."
                         <<std::endl;
                return false;
            }
            if (idx != places_itr->second.first &&
                    idx != places_itr->second.second) {
                std::cout<<"Edge node "<<node<<" appearance at idx "<<idx
                         <<" is not referenced in edge_node_to_places"<<std::endl;
                return false;
            }
        }
    }

    // Verify that num_edge_nodes corresponds to the num of undirected edges.
    size_t num_undirected_edges = 0;  // excludes self-loops
    for (size_t a = 0; a < g.num_nodes(); a++) {
        for (auto b = g.neighbors(a).begin(); b != g.neighbors(a).end(); b++) {
            if (*b > int(a)) {
                num_undirected_edges++;
            }
        }
    }
    if (num_undirected_edges * (1 + int(g.directed)) != g.num_edge_nodes) {
        std::cout<<"num_undirected_edges * (1 + int(g.directed)) != "
                 <<"g.num_edge_nodes"<<std::endl;
        return false;
    }

    // Verify that edge_to_edge_node contains the actual edges.
    for (int a = 0; a < int(g.num_nodes()); a++) {
        for (auto b = g.neighbors(a).begin(); b != g.neighbors(a).end(); b++) {
            if (*b == a) {
                continue;
            }
            if (g.edge_to_edge_node.find(EDGE(a, *b, g.directed)) ==
                    g.edge_to_edge_node.end()) {
                std::cout<<"g.edge_to_edge_node missing an edge"<<std::endl;
                return false;
            }
        }
    }

    // Verify that the edges match the neighbor sets.
    if (g.directed) {
        // Directed
        for (int a = 0; a < int(g.num_nodes()); a++) {
            for (auto b = g.neighbors(a).begin();
                      b != g.neighbors(a).end(); b++) {
                if (a == *b) {
                    // A self-loop. Accounted for elsewhere.
                    continue;
                }
                if (a > *b) {
                    // We already checked this when a < *b.
                    continue;
                }

                auto en_A_itr = g.edge_to_edge_node.find(EDGE(a, *b, true));
                auto en_B_itr = g.edge_to_edge_node.find(EDGE(*b, a, true));
                int edge_node_A = en_A_itr->second; // we know it exists
                int edge_node_B = en_B_itr->second; // we know it exists

                int en_A_startpoint = g.node_to_startpoint[edge_node_A];
                int en_B_startpoint = g.node_to_startpoint[edge_node_B];

                if (g.out_neighbors_vec[en_A_startpoint] != a) {
                    std::cout<<"g.out_neighbors_vec[node_to_startpoint["
                             <<"edge_to_edge_node[(a, b)]]] != a"
                             <<" (directed graph)"<<std::endl;
                    return false;
                }
                if (g.out_neighbors_vec[en_B_startpoint] != *b) {
                    std::cout<<"g.out_neighbors_vec[node_to_startpoint["
                             <<"edge_to_edge_node[(b, a)]]] != b"
                             <<" (directed graph)"<<std::endl;
                    return false;
                }

                if (g.out_neighbors_vec[en_A_startpoint + 1] != edge_node_B) {
                    std::cout<<"g.out_neighbors_vec[node_to_startpoint["
                             <<"edge_to_edge_node[(a, b)]] + 1] != "
                             <<"edge_to_edge_node[(b, a)] (directed graph)"
                             <<std::endl;
                    return false;
                }
                if (g.out_neighbors_vec[en_B_startpoint + 1] != edge_node_A) {
                    std::cout<<"g.out_neighbors_vec[node_to_startpoint["
                             <<"edge_to_edge_node[(b, a)]] + 1] != "
                             <<"edge_to_edge_node[(a, b)] (directed graph)"
                             <<std::endl;
                    return false;
                }

                const auto &en_A_places =
                    g.edge_node_to_places.find(edge_node_A)->second;
                const auto &en_B_places =
                    g.edge_node_to_places.find(edge_node_B)->second;

                if (en_A_places.first < g.node_to_startpoint[a] ||
                        en_A_places.first >= g.node_to_startpoint[a] +
                                                  size_t(g.out_degrees[a])) {
                    std::cout<<"The place where edge_node "<<edge_node_A
                             <<" is referenced is not in "<<a<<"'s range."
                             <<std::endl;
                    return false;
                }
                if (en_B_places.first < g.node_to_startpoint[*b] ||
                        en_B_places.first >= g.node_to_startpoint[*b] +
                                                  size_t(g.out_degrees[*b])) {
                    std::cout<<"The place where edge_node "<<edge_node_B
                             <<" is referenced is not in "<<*b<<"'s range."
                             <<std::endl;
                    return false;
                }

                if (en_A_places.second !=
                        g.node_to_startpoint[edge_node_B] + 1) {
                    std::cout<<"edge_node_to_places[edge_to_edge_node[(a, b)]]"
                             <<".second != node_to_startpoint[edge_to_edge_node"
                             <<"[(b, a)]] + 1 (directed graph)"<<std::endl;
                    return false;
                }
                if (en_B_places.second !=
                        g.node_to_startpoint[edge_node_A] + 1) {
                    std::cout<<"edge_node_to_places[edge_to_edge_node[(b, a)]]"
                             <<".second != node_to_startpoint[edge_to_edge_node"
                             <<"[(a, b)]] + 1 (directed graph)"<<std::endl;
                    return false;
                }
            }
        }
    }
    else {
        // Undirected
        for (int a = 0; a < int(g.num_nodes()); a++) {
            for (auto b = g.neighbors(a).begin();
                      b != g.neighbors(a).end(); b++) {
                if (a == *b) {
                    // A self-loop. Accounted for elsewhere.
                    continue;
                }
                if (a > *b) {
                    // We already checked this when a < *b.
                    continue;
                }

                auto en_itr = g.edge_to_edge_node.find(EDGE(a, *b, false));
                int edge_node = en_itr->second; // we know it exists

                int en_startpoint = g.node_to_startpoint[edge_node];

                if (g.out_neighbors_vec[en_startpoint] != a) {
                    std::cout<<"g.out_neighbors_vec[g.edge_to_edge_node[(a, b)]"
                             <<"] != a (undirected graph)"<<std::endl;
                    return false;
                }
                if (g.out_neighbors_vec[en_startpoint + 1] != *b) {
                    std::cout<<"g.out_neighbors_vec[g.edge_to_edge_node[(a, b)]"
                             <<" + 1] != b (undirected graph)"<<std::endl;
                    return false;
                }

                auto &en_places = g.edge_node_to_places.find(edge_node)->second;
                if (en_places.first < g.node_to_startpoint[a] ||
                            en_places.first >= g.node_to_startpoint[a] +
                                                    size_t(g.out_degrees[a])) {
                    std::cout<<"g.edge_node_to_places[g.edge_to_edge_node[(a, b"
                             <<")]].first not in a's range (undirected graph)"
                             <<" | value: "<<en_places.first<<std::endl;
                    return false;
                }
                if (en_places.second < g.node_to_startpoint[*b] ||
                            en_places.second >= g.node_to_startpoint[*b] +
                                                     size_t(g.out_degrees[*b])) {
                    std::cout<<"g.edge_node_to_places[g.edge_to_edge_node[(a, b"
                             <<")]].second not in b's range (undirected graph)"
                             <<" | value: "<<en_places.second<<std::endl;
                    return false;
                }
            }
        }
    }

    // All tests are passed. Graph is fully consistent.
    return true;
}
#endif

#ifdef SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
std::vector<int> cleaned_out_N_vec(const NTSparseGraph &g) {
    std::vector<int> result = std::vector<int>(g.out_neighbors_vec.size(), 0);

    for (int i = 0; i < int(g.internal_n); i++) {
        for (size_t j = g.node_to_startpoint[i];
                  j < g.node_to_startpoint[i] + size_t(g.out_degrees[i]); j++) {
            result[j] = g.out_neighbors_vec[j];
        }
    }

    return result;
}
#endif

std::string vec_as_string(const std::vector<int> &v) {
    std::string s = "";
    for (size_t i = 0; i < v.size(); i++) {
        s += std::to_string(v[i]) + ", ";
    }
    return s;
}

#ifdef SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
std::string nt_graph_as_string(const NTSparseGraph &g) {
    std::string s = "";

    std::vector<size_t> startpoints = std::vector<size_t>(g.node_to_startpoint);
    std::vector<size_t> endpoints = std::vector<size_t>(g.node_to_endpoint);

    std::sort(startpoints.begin(), startpoints.end());
    std::sort(endpoints.begin(), endpoints.end());

    // std::cout<<"Startpoints: "<<vec_as_string(startpoints)<<std::endl;
    // std::cout<<"Endpoints: "<<vec_as_string(endpoints)<<std::endl;
    // std::cout<<"Vec Size: "<<g.out_neighbors_vec.size()<<std::endl;

    for (size_t idx = 0; idx < endpoints.size(); idx++) {
        auto node_itr = g.endpoint_to_node.find(endpoints[idx]);
        if (node_itr == g.endpoint_to_node.end()) {
            std::cout<<"Missing endpoint info for endpoint "<<endpoints[idx]<<std::endl;
            continue;
        }
        int node = node_itr->second;

        size_t endpoint = g.node_to_endpoint[node];
        if (endpoint != endpoints[idx]) {
            std::cout<<"Error in nt_graph_as_string: node_to_endpoint[endpoint_to_node[ep]] != ep"<<std::endl;
        }
        size_t startpoint = g.node_to_startpoint[node];
        if (startpoint != startpoints[idx]) {
            std::cout<<"Error in nt_graph_as_string: startpoint mis-match"<<std::endl;
        }

        s += std::to_string(node) + " @ " + std::to_string(startpoint) + ": ";
        for (size_t i = startpoint;
                    i < startpoint + size_t(g.out_degrees[node]); i++) {
            s += std::to_string(g.out_neighbors_vec[i]) + ", ";
        }
        s += "| ";
    }
    return s;
}
#endif

#ifdef SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
void print_graph(const NTSparseGraph &g) {
    std::cout<<"A Graph:"<<std::endl;
    std::cout<<"n: "<<g.num_nodes()<<"  m: "<<g.num_edges()
             <<"  internal_n: "<<g.internal_n<<"  num_edge_nodes: "
             <<g.num_edge_nodes<<std::endl;
    std::cout<<vec_as_string(cleaned_out_N_vec(g))<<std::endl;
    std::cout<<nt_graph_as_string(g)<<std::endl;
}
#endif


#ifdef SCHENO__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
// If reconstruct_frequency is greater than 0, then every
//  reconstruct_frequency iterations, the code will replace the graph with a
//  new version of it made with a constructor.
void rand_test(float add_node_prob, float delete_node_prob,
               float add_edge_prob, float delete_edge_prob,
               const bool directed, size_t iterations,
               size_t reconstruct_frequency) {

    const int initial_n = 10;

    // Normalize to sum to 1.
    float prob_sum = add_node_prob + delete_node_prob +
                     add_edge_prob + delete_edge_prob;
    add_node_prob = add_node_prob / prob_sum;
    delete_node_prob = delete_node_prob / prob_sum;
    add_edge_prob = add_edge_prob / prob_sum;
    delete_edge_prob = delete_edge_prob / prob_sum;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dist(0.0, 1.0);

    NTSparseGraph before = NTSparseGraph(directed, initial_n);
    NTSparseGraph after = NTSparseGraph(directed, initial_n);
    NTSparseGraph third = NTSparseGraph(directed, 1);

    int a, b;
    float p;

    size_t skipped_iterations = 0;
    size_t j;

    for (size_t i = 0; i < iterations + skipped_iterations; i++) {
        j = i - skipped_iterations;
        if (reconstruct_frequency != 0 && j > 0 &&
                    j % reconstruct_frequency == 0) {
            // It would be more direct to write after = before;
            //  However, keeping the intermediate object allows us to do
            //  the same thing to before below.
            third = before;
            after = third;
            if (!consistency_check(after)) {
                std::cout<<"Failed after an '=' assignment."<<std::endl;
                std::cout<<"Before"<<std::endl;
                print_graph(before);
                std::cout<<std::endl<<"After"<<std::endl;
                print_graph(after);
                return;
            }
            // I made this clunky in-between object so that
            //  thanks to determinism, after is exactly the same as
            //  before.
            before = third;
        }
        p = dist(gen);
        if (p < add_node_prob) {
            // Add node.

            after.add_node();

            if (!consistency_check(after)) {
                std::cout<<"Failed after an add_node() call."<<std::endl;
                std::cout<<"Before"<<std::endl;
                print_graph(before);
                std::cout<<std::endl<<"After"<<std::endl;
                print_graph(after);
                return;
            }

            before.add_node();

        } else if (p < add_node_prob + delete_node_prob) {
            // Delete node.

            int node = dist(gen) * after.num_nodes();
            if (node == int(after.num_nodes())) {
                node--;
            }

            after.delete_node(node);

            if (!consistency_check(after)) {
                std::cout<<"Failed after a delete_node("<<node<<") call."
                         <<std::endl;
                std::cout<<"Before"<<std::endl;
                print_graph(before);
                std::cout<<std::endl<<"After"<<std::endl;
                print_graph(after);
                std::cout<<"NOTE: In the before graph, node "<<node<<" did "
                         <<(before.neighbors(node).find(node) ==
                            before.neighbors(node).end() ? "not " : "")
                         <<"have a self-loop."<<std::endl;
                return;
            }

            before.delete_node(node);

        } else if (p < add_node_prob + delete_node_prob + add_edge_prob) {
            // Add edge.

            size_t n = after.num_nodes();
            size_t max_num_edges = n + ((n * (n - 1)) / (1 + int(!directed)));

            if (after.num_edges() == max_num_edges) {
                // Edges are already maxed out.
                skipped_iterations++;
                continue;
            }

            if (after.num_edges() >= max_num_edges / sqrt(after.num_nodes())) {
                // Add the j'th possible edge where j is random.
                // O(n)
                size_t num_could_be_added = max_num_edges - after.num_edges();
                int edge_to_add = dist(gen) * num_could_be_added;
                if (edge_to_add == int(num_could_be_added)) {
                    edge_to_add--;
                }

                a = 0;
                int covered_edges = 0;
                while (edge_to_add >= covered_edges +
                        int(after.num_nodes() - after.out_neighbors(a).size())) {
                    covered_edges += after.num_nodes() -
                                     after.out_neighbors(a).size();
                    a++;
                    if (a == int(after.num_nodes())) {
                        std::cout<<"Uh-oh a... "<<edge_to_add<<" vs. "<<covered_edges<<" vs. "<<num_could_be_added<<std::endl;
                    }
                }
                b = 0;
                for (; b < int(after.num_nodes()); b++) {
                    if (after.out_neighbors(a).find(b) == 
                               after.out_neighbors(a).end()) {
                        if (covered_edges == edge_to_add) {
                            break;
                        }
                        covered_edges++;
                    }
                }
                if (b == int(after.num_nodes())) {
                    std::cout<<"Uh-oh..."<<std::endl;
                }
            } else {
                // Add a randomly selected possible edge.
                bool done = false;
                while (!done) {
                    a = dist(gen) * after.num_nodes();
                    if (size_t(a) == after.num_nodes()) {
                        a = a - 1;  // This should never happen, but if it does it's not a problem.
                    }
                    b = dist(gen) * after.num_nodes();
                    if (size_t(b) == after.num_nodes()) {
                        b = b - 1;  // This should never happen, but if it does it's not a problem.
                    }
                    if (after.out_neighbors(a).find(b) == after.out_neighbors(a).end()) {
                        done = true;
                    }
                }
            }

            after.add_edge(a, b);

            if (!consistency_check(after)) {
                std::cout<<"Failed after an add_edge("<<a<<", "<<b<<") call."
                         <<std::endl;
                std::cout<<"Before"<<std::endl;
                print_graph(before);
                std::cout<<std::endl<<"After"<<std::endl;
                print_graph(after);
                return;
            } else if (after.num_edges() != before.num_edges() + 1) {
                std::cout<<"Failed after an add_edge("<<a<<", "<<b<<") call."
                         <<" Num edges did not increment."<<std::endl;
                std::cout<<"Before"<<std::endl;
                print_graph(before);
                std::cout<<std::endl<<"After"<<std::endl;
                print_graph(after);
                return;
            }

            before.add_edge(a, b);

        } else {
            // Delete edge.

            if (after.num_edges() == 0) {
                // No edges to delete.
                skipped_iterations++;
                continue;
            }

            // Delete the j'th edge where j is randomly drawn.
            int edge_to_delete = dist(gen) * after.num_edges();
            if (edge_to_delete == int(after.num_edges())) {
                edge_to_delete = after.num_edges() - 1;
            }

            a = 0;
            int covered_edges = 0;
            while (edge_to_delete >= covered_edges +
                                     int(after.out_neighbors(a).size())) {
                covered_edges += after.out_neighbors(a).size();
                a++;
            }

            if (after.has_self_loop[a] &&
                    dist(gen) < (1.0 / after.out_neighbors(a).size())) {
                b = a;
            } else {
                if (after.has_self_loop[a]) {
                    covered_edges++;
                }

                auto b_itr = after.out_neighbors(a).begin();
                for (int i = 0; i < edge_to_delete - covered_edges; i++) {
                    b_itr++;
                }
                b = *b_itr;
            }

            after.delete_edge(a, b);

            if (!consistency_check(after)) {
                std::cout<<"Failed after a delete_edge("<<a<<", "<<b<<") call."
                         <<std::endl;
                std::cout<<"Before"<<std::endl;
                print_graph(before);
                std::cout<<std::endl<<"After"<<std::endl;
                print_graph(after);
                return;
            } else if (after.num_edges() != before.num_edges() - 1) {
                std::cout<<"Failed after a delete_edge("<<a<<", "<<b<<") call."
                         <<" Num edges did not decrement."<<std::endl;
                std::cout<<"Before"<<std::endl;
                print_graph(before);
                std::cout<<std::endl<<"After"<<std::endl;
                print_graph(after);
                return;
            }

            before.delete_edge(a, b);

        }
    }
    std::cout<<"Finished test successfully with "<<skipped_iterations
             <<" retried iterations."<<std::endl;

    std::cout<<"Graph had "<<after.num_nodes()<<" nodes and "
             <<after.num_edges()<<" edges."<<std::endl;
    std::cout<<"That's "<<float(after.num_edges() * 100) /
           (after.num_nodes() + ((after.num_nodes() * (after.num_nodes() - 1)) /
                                                (1 + int(!after.directed))))
             <<" percent of all possible edges on that many nodes."<<std::endl;
}
#endif
