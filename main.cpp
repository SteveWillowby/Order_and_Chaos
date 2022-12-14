#include<iostream>
#include<string>

#include "file_utils.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"
#include "simulated_annealing_search.h"
#include "sparse_graph.h"

int main( void ) {
    const bool directed = false;
    const std::string edgelist_name = "simple_test_graphs/test_01_edges.txt";
    const std::string nodelist_name = "simple_test_graphs/test_01_nodes.txt";

    std::cout<<"Loading graph from files:"<<std::endl
             <<"    "<<nodelist_name<<std::endl
             <<"    "<<edgelist_name<<std::endl;
    SparseGraph g = read_graph(directed, nodelist_name, edgelist_name);

    auto result = simulated_annealing_search(g, g.num_nodes() * g.num_nodes() *
                                                g.num_nodes(),
                                             9);
    std::cout<<"  ...graph loaded."<<std::endl<<std::endl;

    int i = 0;
    for (auto result_itr = result.rbegin();
              result_itr != result.rend(); result_itr++) {
        std::cout<<"With a score of "<<result_itr->second<<" we have edges: ";

        // Make changes.
        for (auto edge_itr = result_itr->first.begin();
                  edge_itr != result_itr->first.end(); edge_itr++) {
            // g.flip_edge(edge_itr->first, edge_itr->second);
            std::cout<<"("<<edge_itr->first<<", "<<edge_itr->second<<"), ";
        }
        std::cout<<std::endl<<std::endl;


        // Flip back.
        /*
        for (auto edge_itr = result_itr->first.begin();
                  edge_itr != result_itr->first.end(); edge_itr++) {
            g.flip_edge(edge_itr->first, edge_itr->second);
        }
        */

        i++;
    }

    return 0;
}
