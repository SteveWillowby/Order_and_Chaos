#include<iostream>

#include "sparse_graph.h"


int main(void) {
    const bool directed = true;
    SparseGraph g = SparseGraph(directed);

    std::cout<<"All Printed Numbers Should be 1 Unless the Line Specifies Otherwise."<<std::endl<<std::endl;

    std::cout<<(g.num_nodes() == 1)<<std::endl;
    std::cout<<(g.add_node() == 1)<<std::endl;
    std::cout<<(g.num_nodes() == 2)<<std::endl;
    g.add_node();
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    std::cout<<(g.neighbors(1) == std::unordered_set<int>({0, 2}))<<std::endl;
    g.add_node();
    g.add_edge(2, 3);
    std::cout<<(g.delete_node(1) == 3)<<std::endl;
    std::cout<<(g.neighbors(0) == std::unordered_set<int>())<<std::endl;
    std::cout<<(g.neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
    std::cout<<(g.neighbors(2) == std::unordered_set<int>({1}))<<std::endl;

    if (directed) {
        std::cout<<(g.out_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g.in_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g.out_neighbors(1) == std::unordered_set<int>({}))<<std::endl;
        std::cout<<(g.in_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g.out_neighbors(2) == std::unordered_set<int>({1}))<<std::endl;
        std::cout<<(g.in_neighbors(2) == std::unordered_set<int>({}))<<std::endl;
    } else {
        std::cout<<(g.out_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g.in_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g.out_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g.in_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g.out_neighbors(2) == std::unordered_set<int>({1}))<<std::endl;
        std::cout<<(g.in_neighbors(2) == std::unordered_set<int>({1}))<<std::endl;
    }

    g.flip_edge(1, 0);
    g.flip_edge(2, 2);
    std::cout<<(g.neighbors(0) == std::unordered_set<int>({1}))<<std::endl;
    std::cout<<(g.neighbors(1) == std::unordered_set<int>({0, 2}))<<std::endl;
    std::cout<<(g.neighbors(2) == std::unordered_set<int>({1, 2}))<<std::endl;

    if (directed) {
        std::cout<<(g.out_neighbors(0) == std::unordered_set<int>())<<std::endl;
        std::cout<<(g.in_neighbors(0) == std::unordered_set<int>({1}))<<std::endl;
        std::cout<<(g.out_neighbors(1) == std::unordered_set<int>({0}))<<std::endl;
        std::cout<<(g.in_neighbors(1) == std::unordered_set<int>({2}))<<std::endl;
        std::cout<<(g.out_neighbors(2) == std::unordered_set<int>({1, 2}))<<std::endl;
        std::cout<<(g.in_neighbors(2) == std::unordered_set<int>({2}))<<std::endl;
    }

    if (directed) {
        std::cout<<!g.has_edge(0, 0)<<std::endl;
        std::cout<<!g.has_edge(0, 1)<<std::endl;
        std::cout<<!g.has_edge(0, 2)<<std::endl;
        std::cout<< g.has_edge(1, 0)<<std::endl;
        std::cout<<!g.has_edge(1, 1)<<std::endl;
        std::cout<<!g.has_edge(1, 2)<<std::endl;
        std::cout<<!g.has_edge(2, 0)<<std::endl;
        std::cout<< g.has_edge(2, 1)<<std::endl;
        std::cout<< g.has_edge(2, 2)<<std::endl;
    } else {
        std::cout<<!g.has_edge(0, 0)<<std::endl;
        std::cout<< g.has_edge(0, 1)<<std::endl;
        std::cout<<!g.has_edge(0, 2)<<std::endl;
        std::cout<<!g.has_edge(1, 1)<<std::endl;
        std::cout<< g.has_edge(1, 2)<<std::endl;
        std::cout<< g.has_edge(2, 2)<<std::endl;
    }
    std::cout<<(g.num_edges() == 3)<<std::endl;

    SparseGraph g2 = SparseGraph(directed, 5);
    std::cout<<(g2.num_nodes() == 5)<<std::endl;
    std::cout<<(g2.num_edges() == 0)<<std::endl;

    SparseGraph g3 = SparseGraph(g);
    if (directed) {
        std::cout<<!g3.has_edge(0, 0)<<std::endl;
        std::cout<<!g3.has_edge(0, 1)<<std::endl;
        std::cout<<!g3.has_edge(0, 2)<<std::endl;
        std::cout<< g3.has_edge(1, 0)<<std::endl;
        std::cout<<!g3.has_edge(1, 1)<<std::endl;
        std::cout<<!g3.has_edge(1, 2)<<std::endl;
        std::cout<<!g3.has_edge(2, 0)<<std::endl;
        std::cout<< g3.has_edge(2, 1)<<std::endl;
        std::cout<< g3.has_edge(2, 2)<<std::endl;
    } else {
        std::cout<<!g3.has_edge(0, 0)<<std::endl;
        std::cout<< g3.has_edge(0, 1)<<std::endl;
        std::cout<<!g3.has_edge(0, 2)<<std::endl;
        std::cout<<!g3.has_edge(1, 1)<<std::endl;
        std::cout<< g3.has_edge(1, 2)<<std::endl;
        std::cout<< g3.has_edge(2, 2)<<std::endl;
    }
    std::cout<<(g3.num_nodes() == 3)<<std::endl;
    std::cout<<(g3.num_edges() == 3)<<std::endl;

    return 0;
};
