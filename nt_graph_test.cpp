#include<algorithm>
#include<iostream>

#include "sparse_graph.h"
#include "nt_sparse_graph.h"

std::vector<int> cleaned_out_N_vec(const NTSparseGraph &g) {
    std::vector<int> result = std::vector<int>(g.out_neighbors_vec.size(), 0);

    for (int i = 0; i < int(g.internal_n); i++) {
        for (int j = g.node_to_startpoint[i];
                    j < g.node_to_startpoint[i] + g.out_degrees[i]; j++) {
            result[j] = g.out_neighbors_vec[j];
        }
    }

    return result;
}

std::string vec_as_string(const std::vector<int> &v) {
    std::string s = "";
    for (size_t i = 0; i < v.size(); i++) {
        s += std::to_string(v[i]) + ", ";
    }
    return s;
}

std::string nt_graph_as_string(const NTSparseGraph &g) {
    std::string s = "";

    std::vector<int> startpoints = std::vector<int>(g.node_to_startpoint);
    std::vector<int> endpoints = std::vector<int>(g.node_to_endpoint);

    std::sort(startpoints.begin(), startpoints.end());
    std::sort(endpoints.begin(), endpoints.end());

    // std::cout<<"Startpoints: "<<vec_as_string(startpoints)<<std::endl;
    // std::cout<<"Endpoints: "<<vec_as_string(endpoints)<<std::endl;
    // std::cout<<"Vec Size: "<<g.out_neighbors_vec.size()<<std::endl;

    for (size_t idx = 0; idx < startpoints.size(); idx++) {
        auto node_itr = g.endpoint_to_node.find(endpoints[idx]);
        if (node_itr == g.endpoint_to_node.end()) {
            std::cout<<"Missing endpoint info for endpoint "<<endpoints[idx]<<std::endl;
            continue;
        }
        int node = node_itr->second;

        s += std::to_string(node) + " @ " + std::to_string(startpoints[idx]) + ": ";
        for (int i = startpoints[idx]; i < endpoints[idx]; i++) {
            s += std::to_string(g.out_neighbors_vec[i]) + ", ";
        }
        s += "| ";
    }
    return s;
}

int main(void) {

    std::cout<<"All Printed Numbers Should be 1 Unless the Line Specifies Otherwise."<<std::endl<<std::endl;

    const bool directed = true;
    NTSparseGraph g1 = NTSparseGraph(directed);

    g1.add_node();
    g1.add_node();
    std::cout<<(g1.out_neighbors_vec == std::vector<int>({0,0,0,0, 0,0,0,0, 0,0,0,0}))<<std::endl;
    g1.add_edge(0, 1);
    std::cout<<(g1.out_neighbors_vec == std::vector<int>({3,0,0,0, 4,0,0,0, 0,0,0,0, 0,4, 1,3}))<<std::endl;
    g1.add_edge(2, 0);
    std::cout<<(g1.out_neighbors_vec == std::vector<int>({3,6,0,0, 4,0,0,0, 5,0,0,0, 0,4, 1,3, 2,6, 0,5}))<<std::endl;
    g1.add_edge(0, 2);
    std::cout<<(g1.out_neighbors_vec == std::vector<int>({3,6,0,0, 4,0,0,0, 5,0,0,0, 0,4, 1,3, 2,6, 0,5}))<<std::endl;
    g1.add_node();
    std::vector<int> expected = std::vector<int>({7,6,0,0, 4,0,0,0, 5,0,0,0, 0,0,0,0, 2,6, 0,5, 0,4, 1,7});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    // std::cout<<vec_as_string(expected)<<std::endl;
    // std::cout<<vec_as_string(cleaned_out_N_vec(g1))<<std::endl;
    // std::cout<<nt_graph_as_string(g1)<<std::endl;

    return 0;

    NTSparseGraph g2 = NTSparseGraph(directed);

    std::cout<<(g2.num_nodes() == 1)<<std::endl;
    std::cout<<(g2.add_node() == 1)<<std::endl;
    std::cout<<(g2.num_nodes() == 2)<<std::endl;
    g2.add_node();
    g2.add_edge(0, 1);
    g2.add_edge(1, 2);
    std::cout<<(g2.neighbors(1) == std::unordered_set<int>({0, 2}))<<std::endl;
    g2.add_node();
    g2.add_edge(2, 3);
    std::cout<<(g2.delete_node(1) == 3)<<std::endl; // TODO: delete_node unimplemented
    std::cout<<(g2.neighbors(0) == std::unordered_set<int>())<<std::endl;
    std::cout<<(g2.neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
    std::cout<<(g2.neighbors(2) == std::unordered_set<int>({1}))<<std::endl;

    std::cout<<"Hola"<<std::endl;

    if (directed) {
        std::cout<<(g2.out_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g2.in_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g2.out_neighbors(1) == std::unordered_set<int>({}))<<std::endl;
        std::cout<<(g2.in_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g2.out_neighbors(2) == std::unordered_set<int>({1}))<<std::endl;
        std::cout<<(g2.in_neighbors(2) == std::unordered_set<int>({}))<<std::endl;
    } else {
        std::cout<<(g2.out_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g2.in_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g2.out_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g2.in_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g2.out_neighbors(2) == std::unordered_set<int>({1}))<<std::endl;
        std::cout<<(g2.in_neighbors(2) == std::unordered_set<int>({1}))<<std::endl;
    }

    std::cout<<"Como Estas"<<std::endl;

    g2.flip_edge(1, 0);

    std::cout<<"Muy Bien"<<std::endl;

    g2.flip_edge(2, 2);

    std::cout<<(g2.neighbors(0) == std::unordered_set<int>({1}))<<std::endl;
    std::cout<<(g2.neighbors(1) == std::unordered_set<int>({0, 2}))<<std::endl;
    std::cout<<(g2.neighbors(2) == std::unordered_set<int>({1, 2}))<<std::endl;

    if (directed) {
        std::cout<<(g2.out_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g2.in_neighbors(0) == std::unordered_set<int>({1}))<<std::endl;
        std::cout<<(g2.out_neighbors(1) == std::unordered_set<int>({0}))<<std::endl;
        std::cout<<(g2.in_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g2.out_neighbors(2) == std::unordered_set<int>({1, 2}))<<std::endl;
        std::cout<<(g2.in_neighbors(2) == std::unordered_set<int>({2}))<<std::endl;
    }

    if (directed) {
        std::cout<<!g2.has_edge(0, 0)<<std::endl;
        std::cout<<!g2.has_edge(0, 1)<<std::endl;
        std::cout<<!g2.has_edge(0, 2)<<std::endl;
        std::cout<< g2.has_edge(1, 0)<<std::endl;
        std::cout<<!g2.has_edge(1, 1)<<std::endl;
        std::cout<<!g2.has_edge(1, 2)<<std::endl;
        std::cout<<!g2.has_edge(2, 0)<<std::endl;
        std::cout<< g2.has_edge(2, 1)<<std::endl;
        std::cout<< g2.has_edge(2, 2)<<std::endl;
    } else {
        std::cout<<!g2.has_edge(0, 0)<<std::endl;
        std::cout<< g2.has_edge(0, 1)<<std::endl;
        std::cout<<!g2.has_edge(0, 2)<<std::endl;
        std::cout<<!g2.has_edge(1, 1)<<std::endl;
        std::cout<< g2.has_edge(1, 2)<<std::endl;
        std::cout<< g2.has_edge(2, 2)<<std::endl;
    }
    std::cout<<(g2.num_edges() == 3)<<std::endl;

    SparseGraph g3 = SparseGraph(directed, 5);
    std::cout<<(g3.num_nodes() == 5)<<std::endl;
    std::cout<<(g3.num_edges() == 0)<<std::endl;

    SparseGraph g4 = SparseGraph(g2);
    if (directed) {
        std::cout<<!g4.has_edge(0, 0)<<std::endl;
        std::cout<<!g4.has_edge(0, 1)<<std::endl;
        std::cout<<!g4.has_edge(0, 2)<<std::endl;
        std::cout<< g4.has_edge(1, 0)<<std::endl;
        std::cout<<!g4.has_edge(1, 1)<<std::endl;
        std::cout<<!g4.has_edge(1, 2)<<std::endl;
        std::cout<<!g4.has_edge(2, 0)<<std::endl;
        std::cout<< g4.has_edge(2, 1)<<std::endl;
        std::cout<< g4.has_edge(2, 2)<<std::endl;
    } else {
        std::cout<<!g4.has_edge(0, 0)<<std::endl;
        std::cout<< g4.has_edge(0, 1)<<std::endl;
        std::cout<<!g4.has_edge(0, 2)<<std::endl;
        std::cout<<!g4.has_edge(1, 1)<<std::endl;
        std::cout<< g4.has_edge(1, 2)<<std::endl;
        std::cout<< g4.has_edge(2, 2)<<std::endl;
    }
    std::cout<<(g4.num_nodes() == 3)<<std::endl;
    std::cout<<(g4.num_edges() == 3)<<std::endl;

    #ifdef SYM__SPARSE_GRAPH_INCLUDE_ERROR_CHECKS
    std::cout<<"The next thing you should see is a custom error message."<<std::endl;
    #else
    std::cout<<"The next thing you should see is a generic error message."<<std::endl;
    #endif
    g4.add_edge(0, 3);

    return 0;
};
