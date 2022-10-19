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
        for (int i = startpoints[idx]; i < startpoints[idx] + g.out_degrees[node]; i++) {
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
    g1.add_node();
    // Nodes                     0        1        2        3        4        7    8    5    6
    expected = std::vector<int>({7,6,0,0, 8,0,0,0, 5,0,0,0, 0,0,0,0, 0,0,0,0, 0,8, 1,7, 2,6, 0,5});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    g1.add_edge(3, 0);
    // Nodes                     0         1        2        3        4        7    8    5    6    9     10
    expected = std::vector<int>({7,6,10,0, 8,0,0,0, 5,0,0,0, 9,0,0,0, 0,0,0,0, 0,8, 1,7, 2,6, 0,5, 3,10, 0,9});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    g1.add_edge(0, 4);
    // Nodes                     0          1        2        3        4         7    8    5    6    9     10   11    12
    expected = std::vector<int>({7,6,10,11, 8,0,0,0, 5,0,0,0, 9,0,0,0, 12,0,0,0, 0,8, 1,7, 2,6, 0,5, 3,10, 0,9, 0,12, 4,11});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    g1.add_node();
    // Nodes                     0          1        2         3        4         5        13   6     9     10   11    12    7    8
    expected = std::vector<int>({7,6,10,11, 8,0,0,0, 13,0,0,0, 9,0,0,0, 12,0,0,0, 0,0,0,0, 2,6, 0,13, 3,10, 0,9, 0,12, 4,11, 0,8, 1,7});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    g1.add_edge(5, 0);
    expected = std::vector<int>(
    // Nodes -1       1        2         3        4         5         0                   11    12    7    8    13   6     9     10   14    15
            {0,0,0,0, 8,0,0,0, 13,0,0,0, 9,0,0,0, 12,0,0,0, 14,0,0,0, 7,6,10,11,15,0,0,0, 0,12, 4,11, 0,8, 1,7, 2,6, 0,13, 3,10, 0,9, 5,15, 0,14});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    g1.add_node();
    expected = std::vector<int>(
// Nodes 6        1        2         3        4         5         0                    11    12    7    8    13    16    9     10   14    15
        {0,0,0,0, 8,0,0,0, 13,0,0,0, 9,0,0,0, 12,0,0,0, 14,0,0,0, 7,16,10,11,15,0,0,0, 0,12, 4,11, 0,8, 1,7, 2,16, 0,13, 3,10, 0,9, 5,15, 0,14});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    // Now testing edge deletion.

    g1.delete_edge(5, 0);
    expected = std::vector<int>(
// Nodes 6        1        2         3        4         5        0                   11    12    7    8    13    14    9     10
        {0,0,0,0, 8,0,0,0, 13,0,0,0, 9,0,0,0, 12,0,0,0, 0,0,0,0, 7,14,10,11,0,0,0,0, 0,12, 4,11, 0,8, 1,7, 2,14, 0,13, 3,10, 0,9});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.delete_edge(1, 0); // The edge does not exist - expect unchanged vector.
    expected = std::vector<int>(
// Nodes 6        1        2         3        4         5        0                   11    12    7    8    13    14    9     10
        {0,0,0,0, 8,0,0,0, 13,0,0,0, 9,0,0,0, 12,0,0,0, 0,0,0,0, 7,14,10,11,0,0,0,0, 0,12, 4,11, 0,8, 1,7, 2,14, 0,13, 3,10, 0,9});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.delete_edge(0, 1);
    expected = std::vector<int>(
// Nodes 6        1        2        3        4         5        0                  11    12    10   9     8    7
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 9,0,0,0, 12,0,0,0, 0,0,0,0, 11,7,10,0,0,0,0,0, 0,12, 4,11, 0,9, 3,10, 2,7, 0,8});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    /*
    // A check to ensure that a node isn't added into empty space too soon.
    //  This check does not fit within the main test sequence thread --
    //  hence the return 0;
    g1.add_node();
    expected = std::vector<int>(
// Nodes 6        1        2        3        4         5        0                   7        10   9     8     13   11    12
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 9,0,0,0, 12,0,0,0, 0,0,0,0, 11,13,10,0,0,0,0,0, 0,0,0,0, 0,9, 3,10, 2,13, 0,8, 0,12, 4,11});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    return 0;
    */

    g1.delete_edge(0, 4);
    expected = std::vector<int>(
// Nodes 6        1        2        3        4        5        0                 7    8    10   9
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 9,0,0,0, 0,0,0,0, 0,0,0,0, 10,7,0,0,0,0,0,0, 0,8, 2,7, 0,9, 3,10});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.add_node();
    expected = std::vector<int>(
// Nodes 6        1        2        3        4        5        0          7        11   8     10   9
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 9,0,0,0, 0,0,0,0, 0,0,0,0, 10,11,0,0, 0,0,0,0, 0,8, 2,11, 0,9, 3,10});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.add_edge(2, 7);
    expected = std::vector<int>(
// Nodes 6        1        2         3        4        5        0          7         11   8     10   9     12    13
        {0,0,0,0, 0,0,0,0, 8,12,0,0, 9,0,0,0, 0,0,0,0, 0,0,0,0, 10,11,0,0, 13,0,0,0, 0,8, 2,11, 0,9, 3,10, 2,13, 7,12});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.delete_edge(2, 7);
    expected = std::vector<int>(
// Nodes 6        1        2        3        4        5        0          7        11   8     10   9
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 9,0,0,0, 0,0,0,0, 0,0,0,0, 10,11,0,0, 0,0,0,0, 0,8, 2,11, 0,9, 3,10});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.delete_edge(3, 0);
    expected = std::vector<int>(
// Nodes 6        1        2        3        4        5        0        7        9    8
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 9,0,0,0, 0,0,0,0, 0,8, 2,9});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.add_edge(3, 0);
    expected = std::vector<int>(
// Nodes 6        1        2        3         4        5        0         7        9    8    10    11
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 10,0,0,0, 0,0,0,0, 0,0,0,0, 9,11,0,0, 0,0,0,0, 0,8, 2,9, 3,11, 0,10});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.add_edge(0, 7);
    expected = std::vector<int>(
// Nodes 6        1        2        3         4        5        0          7         9    8    10    11    12    13
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 10,0,0,0, 0,0,0,0, 0,0,0,0, 9,11,12,0, 13,0,0,0, 0,8, 2,9, 3,11, 0,10, 0,13, 7,12});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.add_edge(5, 0);
    expected = std::vector<int>(
// Nodes 6        1        2        3         4        5         0           7         9    8    10    11    12    13    14    15
        {0,0,0,0, 0,0,0,0, 8,0,0,0, 10,0,0,0, 0,0,0,0, 14,0,0,0, 9,11,12,15, 13,0,0,0, 0,8, 2,9, 3,11, 0,10, 0,13, 7,12, 5,15, 0,14});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.add_edge(6, 0);
    expected = std::vector<int>(
// Nodes 6         1        2        3         4        5                 7         0                   
        {16,0,0,0, 0,0,0,0, 8,0,0,0, 10,0,0,0, 0,0,0,0, 14,0,0,0,0,0,0,0, 13,0,0,0, 9,11,12,15,17,0,0,0,

//ENodes 12    13    14    15    9    8    10    11    16    17
         0,13, 7,12, 5,15, 0,14, 0,8, 2,9, 3,11, 0,10, 6,17, 0,16});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    /*
    //  This check does not fit within the main test sequence thread --
    //  hence the return 0;
    g1.add_node();
    expected = std::vector<int>(
// Nodes 6         1        2         3         4        5         8        7         0                   
        {16,0,0,0, 0,0,0,0, 18,0,0,0, 10,0,0,0, 0,0,0,0, 14,0,0,0, 0,0,0,0, 13,0,0,0, 9,11,12,15,17,0,0,0,

//ENodes 12    13    14    15    9     18   10    11    16    17
         0,13, 7,12, 5,15, 0,14, 0,18, 2,9, 3,11, 0,10, 6,17, 0,16});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    return 0;
    */

    g1.add_edge(5, 4);
    expected = std::vector<int>(
// Nodes 6         1        2        3         4         5                  7         0                   
        {16,0,0,0, 0,0,0,0, 8,0,0,0, 10,0,0,0, 19,0,0,0, 14,18,0,0,0,0,0,0, 13,0,0,0, 9,11,12,15,17,0,0,0,

//ENodes 12    13    14    15    9    8    10    11    16    17    18    19
         0,13, 7,12, 5,15, 0,14, 0,8, 2,9, 3,11, 0,10, 6,17, 0,16, 5,19, 4,18});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.add_edge(5, 3);
    expected = std::vector<int>(
// Nodes 6         1        2        3          4         5                   7         0                   
        {16,0,0,0, 0,0,0,0, 8,0,0,0, 10,21,0,0, 19,0,0,0, 14,18,20,0,0,0,0,0, 13,0,0,0, 9,11,12,15,17,0,0,0,

//ENodes 12    13    14    15    9    8    10    11    16    17    18    19    20    21
         0,13, 7,12, 5,15, 0,14, 0,8, 2,9, 3,11, 0,10, 6,17, 0,16, 5,19, 4,18, 5,21, 3,20});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    /*
    //  This check does not fit within the main test sequence thread --
    //  hence the return 0;
    g1.add_node();
    expected = std::vector<int>(
// Nodes 6         1        2         3          4         5                   7         0                    8
        {16,0,0,0, 0,0,0,0, 22,0,0,0, 10,21,0,0, 19,0,0,0, 14,18,20,0,0,0,0,0, 13,0,0,0, 9,11,12,15,17,0,0,0, 0,0,0,0,

//ENodes 14    15    9     22   10    11    16    17    18    19    20    21    12    13
         5,15, 0,14, 0,22, 2,9, 3,11, 0,10, 6,17, 0,16, 5,19, 4,18, 5,21, 3,20, 0,13, 7,12});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;
    return 0;
    */

    // Both edges (2, 0) and (0, 2) are currently present.
    g1.delete_edge(2, 0);
    expected = std::vector<int>(
// Nodes 6         1        2        3          4         5                   7         0                   
        {16,0,0,0, 0,0,0,0, 8,0,0,0, 10,21,0,0, 19,0,0,0, 14,18,20,0,0,0,0,0, 13,0,0,0, 9,11,12,15,17,0,0,0,

//ENodes 12    13    14    15    9    8    10    11    16    17    18    19    20    21
         0,13, 7,12, 5,15, 0,14, 0,8, 2,9, 3,11, 0,10, 6,17, 0,16, 5,19, 4,18, 5,21, 3,20});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    g1.delete_edge(0, 2);
    expected = std::vector<int>(
// Nodes 6         1        2        3         4         5                  7         0                   
        {16,0,0,0, 0,0,0,0, 0,0,0,0, 10,9,0,0, 19,0,0,0, 14,18,8,0,0,0,0,0, 13,0,0,0, 17,11,12,15,0,0,0,0,

//ENodes 12    13    14    15    9    8    10    11    16    17    18    19
         0,13, 7,12, 5,15, 0,14, 3,8, 5,9, 3,11, 0,10, 6,17, 0,16, 5,19, 4,18});
    std::cout<<(cleaned_out_N_vec(g1) == expected)<<std::endl;

    std::cout<<vec_as_string(expected)<<std::endl;
    std::cout<<vec_as_string(cleaned_out_N_vec(g1))<<std::endl;
    std::cout<<std::endl<<std::endl;
    std::cout<<nt_graph_as_string(g1)<<std::endl;

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
