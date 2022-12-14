#include<iostream>
#include<string>

#include "file_utils.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"
#include "simulated_annealing_search.h"
#include "sparse_graph.h"

int main( void ) {
    // These three variables determine which graph is run on.
    const bool directed = true;
    const bool use_real_graph = true;
    const size_t graph_idx = 0;

    std::vector<std::string> fake_nodes_names =
        {"test_01_nodes.txt", "test_02_nodes.txt"
        };
    std::vector<std::string> fake_edges_names =
        {"test_01_edges.txt", "test_02_edges.txt"
        };
    std::vector<std::string> real_nodes_names =
        {"",                     ""
        };
    std::vector<std::string> real_edges_names =
        {"celegans_metabolic.g", "species_brain_1.g"
        };

    const std::string fake_prefix = "simple_test_graphs/";
    for (size_t i = 0; i < fake_nodes_names.size(); i++) {
        if (!fake_nodes_names[i].empty()) {
            fake_nodes_names[i] = fake_prefix[i] + fake_nodes_names[i];
        }
        fake_edges_names[i] = fake_prefix + fake_edges_names[i];
    }
    const std::string real_prefix = "real_world_graphs/";
    for (size_t i = 0; i < real_nodes_names.size(); i++) {
        if (!real_nodes_names[i].empty()) {
            real_nodes_names[i] = real_prefix[i] + real_nodes_names[i];
        }
        real_edges_names[i] = real_prefix + real_edges_names[i];
    }

    std::string edgelist_name;
    std::string nodelist_name;
    if (use_real_graph) {
        edgelist_name = real_edges_names[graph_idx];
        nodelist_name = real_nodes_names[graph_idx];
    } else {
        edgelist_name = fake_edges_names[graph_idx];
        nodelist_name = fake_nodes_names[graph_idx];
    }

    std::cout<<"Loading graph from files:"<<std::endl
             <<"    "<<nodelist_name<<std::endl
             <<"    "<<edgelist_name<<std::endl;
    SparseGraph g(directed);
    if (nodelist_name.empty()) {
        g = read_graph(directed, edgelist_name);
    } else {
        g = read_graph(directed, nodelist_name, edgelist_name);
    }
    std::cout<<"  ...graph loaded. It has "<<g.num_nodes()<<" nodes and "
             <<g.num_edges()<<" edges, "<<g.num_loops()<<" of which are "
             <<"self-loops."<<std::endl<<std::endl;

    size_t num_iterations = g.num_nodes() * g.num_nodes() * g.num_nodes();

    std::cout<<"Running for "<<num_iterations<<" iterations..."<<std::endl;
    auto result = simulated_annealing_search(g, num_iterations, 9);

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
