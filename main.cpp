#include<iostream>
#include<string>

#include "file_utils.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"
#include "simulated_annealing_search.h"
#include "sparse_graph.h"

#include "edge.h" // TODO: Remove
#include<unordered_set>  // TODO: Remove

int main( void ) {
    // These three variables determine which graph is run on.
    const bool directed = false;
    const bool use_real_graph = true;
    const size_t graph_idx = 2;
    // Note: Took 252 minutes for jazz collab with n^3 iterations.

    std::vector<std::string> fake_nodes_names =
        {"test_01_nodes.txt", "test_02_nodes.txt"
        };
    std::vector<std::string> fake_edges_names =
        {"test_01_edges.txt", "test_02_edges.txt"
        };
    std::vector<std::string> real_nodes_names =
        {"",                               "",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt"
        };
    // jazz_collab_mod_X.g are various man-made modifications to
    //  jazz_collaboration to try to increase the amount of symmetry and to see
    //  how the scoring function would score those modifications.
    // jazz_collab_mod_X_changes.txt contains the list of edges removed.
    std::vector<std::string> real_edges_names =
        {"celegans_metabolic.g",           "species_brain_1.g",
         "jazz_collaboration.g",           "jazz_collab_mod_1.g",
         "jazz_collab_mod_2.g",            "jazz_collab_mod_4.g",
         "jazz_collab_mod_10.g",           "jazz_collab_mod_1_changes.txt",
         "jazz_collab_mod_2_changes.txt",  "jazz_collab_mod_4_changes.txt",
         "jazz_collab_mod_10_changes.txt"
        };

    const std::string fake_prefix = "simple_test_graphs/";
    for (size_t i = 0; i < fake_nodes_names.size(); i++) {
        if (!fake_nodes_names[i].empty()) {
            fake_nodes_names[i] = fake_prefix + fake_nodes_names[i];
        }
        fake_edges_names[i] = fake_prefix + fake_edges_names[i];
    }
    const std::string real_prefix = "real_world_graphs/";
    for (size_t i = 0; i < real_nodes_names.size(); i++) {
        if (!real_nodes_names[i].empty()) {
            real_nodes_names[i] = real_prefix + real_nodes_names[i];
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

    // TODO: Remove this code.
    SparseGraph candidate_changes(directed);
    candidate_changes = read_graph(directed, "real_world_graphs/jazz_collab_nodes.txt",
                                             "real_world_graphs/jazz_collab_mod_10_changes.txt");
    std::unordered_set<Edge, EdgeHash> cc;
    for (size_t a = 0; a < candidate_changes.num_nodes(); a++) {
        for (size_t b = (!directed) * (a + 1); b < candidate_changes.num_nodes(); b++) {
            if (b == a) {
                continue;
            }
            if (candidate_changes.has_edge(a, b)) {
                cc.insert(EDGE(a, b, directed));
            }
        }
    }
    // END TODO

    size_t num_iterations = g.num_nodes() * g.num_nodes() *
                            g.num_nodes();

    std::cout<<"Running for "<<num_iterations<<" iterations..."<<std::endl;
    auto result = simulated_annealing_search(g, num_iterations, 9, cc);

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
