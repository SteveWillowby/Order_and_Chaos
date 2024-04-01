#include "nt_wrappers/nauty_traces.h"

#include<algorithm>
#include<cmath>
#include<iostream>
#include<string>
#include<unordered_set>
#include<utility>
#include<vector>

std::vector<int> long_to_indices(unsigned long l) {
    std::vector<int> result = {};
    for (int i = 0; i < int(sizeof(unsigned long)) * 8; i++) {
        if (l & 0x1) {
            result.push_back(i);
        }
        l >>= 1;
    }
    return result;
}

std::pair<int, int> index_to_edge(int i, int num_nodes, bool directed) {
    if (directed) {
        int a = i / (num_nodes - 1);
        int b = i % (num_nodes - 1);
        if (b >= a) {
            b++;
        }
        return {a, b};
    }

    int a = 0;
    int x = (num_nodes - 1);
    while (x <= i) {
        a++;
        x += (num_nodes - 1) - a;
    }
    x -= (num_nodes - 1) - a;
    return {a, a + (i - x) + 1};
}

std::vector<std::pair<int, int>> number_to_edges(unsigned long num, int num_nodes,
                                                 bool directed) {
    std::vector<std::pair<int, int>> result = {};
    std::vector<int> indices = std::move(long_to_indices(num));
    for (auto i_itr = indices.begin(); i_itr != indices.end(); i_itr++) {
        result.push_back(std::move(index_to_edge(*i_itr, num_nodes, directed)));
    }
    return result;
}

std::vector<int> node_map(const std::vector<int>& node_order) {
    std::vector<int> res = std::vector<int>(node_order.size(), 0);
    for (size_t i = 0; i < node_order.size(); i++) {
        res[node_order[i]] = i;
    }
    return res;
}

std::vector<std::pair<int, int>> canonical_form(
                const std::vector<std::pair<int, int>>& edges,
                const std::vector<int>& node_map, bool directed) {
    std::vector<std::pair<int, int>> res;
    for (auto e_itr = edges.begin(); e_itr != edges.end(); e_itr++) {
        int a = node_map[e_itr->first];
        int b = node_map[e_itr->second];
        if (directed || a <= b) {
            res.push_back({a, b});
        } else {
            res.push_back({b, a});
        }
    }
    std::sort(res.begin(), res.end());

    return res;
}

std::string debug_str(const std::vector<std::pair<int, int>>& edges) {
    std::string res = "";
    for (auto e_itr = edges.begin(); e_itr != edges.end(); e_itr++) {
        res += "(" + std::to_string(e_itr->first) + ", "
                   + std::to_string(e_itr->second) + ") ";
    }
    return res;
}

std::string vec_str(const std::vector<int>& vec) {
    std::string res = "";
    for (auto v_itr = vec.begin(); v_itr != vec.end(); v_itr++) {
        res += std::to_string(*v_itr) + " ";
    }
    return res;
}

std::string vec_pair_str(const std::vector<std::pair<int, int>>& vec) {
    std::string res = "";
    for (auto v_itr = vec.begin(); v_itr != vec.end(); v_itr++) {
        res += std::to_string(v_itr->first) + " " +
               std::to_string(v_itr->second) + " ";
    }
    return res;
}

int main(void) {

    const bool DIRECTED = 1;
    const int MAX_NUM_NODES = (DIRECTED ? 8 : 11);

    NautyTracesOptions nto;
    nto.get_edge_orbits = false;
    nto.get_node_orbits = false;
    nto.get_canonical_node_order = true;

    std::unordered_set<std::string> graphs;

    for (int N = 3; N <= MAX_NUM_NODES; N++) {
        unsigned long max_E = ((N * (N - 1)) / (2 - DIRECTED));
        unsigned long max_l = ((unsigned long) 1) << max_E;
        long num_graphs = 0;

        long double sum_aut = 0;

        graphs.clear();
        for (unsigned long l = 0; l < max_l; l++) {
            std::vector<std::pair<int, int>> edges =
                std::move(number_to_edges(l, N, DIRECTED));

            // TODO: Get graph info
            NTSparseGraph g(DIRECTED, N);
            for (auto e_itr = edges.begin(); e_itr != edges.end(); e_itr++) {
                g.add_edge(e_itr->first, e_itr->second);
            }
            NautyTracesResults ntr = traces(g, nto);

            std::vector<int> map = std::move(node_map(ntr.canonical_node_order));

            // std::cout<<vec_str(ntr.canonical_node_order)<<std::endl;
            // std::cout<<vec_str(map)<<std::endl;
            // std::cout<<debug_str(canonical_form(edges, map, DIRECTED))<<std::endl;

            std::string canon =
                vec_pair_str(std::move(canonical_form(edges, map, DIRECTED)));

            if (graphs.find(canon) == graphs.end()) {
                graphs.insert(canon);
                num_graphs++;

                sum_aut += ((long double) ntr.num_aut_base) *
                           (std::pow((long double) 10.0,
                                     (long double) ntr.num_aut_exponent));
            }
        }

        std::cout<<std::endl;
        std::cout<<"Number of "<<N<<"-node graphs: "<<num_graphs<<std::endl;
        std::cout<<"Number of automorphisms: "<<sum_aut<<std::endl;

        long double num_mat_thing = std::exp2((long double) max_E);
        for (size_t x = 2; x <= size_t(N); x++) {
            num_mat_thing /= (long double) x;
        }
        long double fancy_factor = num_graphs / num_mat_thing;
        long double other_factor = sum_aut / num_mat_thing;

        std::cout<<"Power: "<<(std::log2(other_factor) /
                               std::log2(fancy_factor))<<std::endl;

        // TODO: Aggregate and print graph info
    }

    return 0;
}
