#include<cmath>
#include<iostream>
#include<random>
#include<string>

#include "edge.h"
#include "edge_sampler.h"
#include "file_utils.h"
#include "genetic_alg_search.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"
#include "simulated_annealing_search.h"
#include "sparse_graph.h"

int main( void ) {
    // These three variables determine which graph is run on.
    const bool DIRECTED = false;
    const bool use_real_graph = false;
    const size_t graph_idx = 5;

    const size_t top_k = 9;  // Number of candidate noise sets to keep.
    const size_t ITERS_PER_FLIP_PROB = 5;

    // const bool corrupt_original = false;
    // Only used when corrupt_original is true.
    // const size_t num_additions = 0;
    // Only used when corrupt_original is true.
    // const size_t num_removals = 0;

    const std::vector<float> FLIP_PROBS = {0.0, 0.005, 0.01, 0.025, 0.05};

    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;

    NautyTracesResults nt_results;

    std::random_device rd;  // Will provide a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    std::uniform_real_distribution<float> dist(0, 1);

    // test_01 -- a 5-chain
    // test_02 -- a 10-chain
    // test_03 -- a 15 node binary tree
    // test_04 -- a 4x5 grid
    // test_05 -- a 6x7 grid
    // test_06 -- a 31 node binary tree

    std::vector<std::string> fake_nodes_names =
        {"test_01_nodes.txt", "test_02_nodes.txt",
         "test_03_nodes.txt", "test_04_nodes.txt",
         "test_05_nodes.txt", "test_06_nodes.txt"
        };
    std::vector<std::string> fake_edges_names =
        {"test_01_edges.txt", "test_02_edges.txt",
         "test_03_edges.txt", "test_04_edges.txt",
         "test_05_edges.txt", "test_06_edges.txt"
        };
    std::vector<std::string> real_nodes_names =
        {"",                               "",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "",
         "",                               ""
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
         "jazz_collab_mod_10_changes.txt", "moreno_highschool_noweights.g",
         "roget_thesaurus.g",              "pol_blogs.g"
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
    SparseGraph g(DIRECTED);
    if (nodelist_name.empty()) {
        g = read_graph(DIRECTED, edgelist_name);
    } else {
        g = read_graph(DIRECTED, nodelist_name, edgelist_name);
    }
    std::cout<<"  ...graph loaded. It has "<<g.num_nodes()<<" nodes and "
             <<g.num_edges()<<" edges, "<<g.num_loops()<<" of which are "
             <<"self-loops."<<std::endl<<std::endl;

    bool has_self_loops = g.num_loops() > 0;

    NTSparseGraph g_info = NTSparseGraph(g);
    nt_results = traces(g_info, o);
    double log2_aut = std::log2l(nt_results.num_aut_base) +
                     ((long double)(nt_results.num_aut_exponent)) *
                                      std::log2l(10);
    std::cout<<"The original graph has log2_aut = "<<log2_aut<<std::endl;

    std::unordered_set<Edge, EdgeHash> random_deletions, random_additions;

    // For genetic alg
    const size_t NUM_ITERATIONS = 30;
    const size_t NUM_THREADS = 0;

    std::cout<<"Running for "<<NUM_ITERATIONS<<" iterations..."<<std::endl;

    for (auto flip_prob = FLIP_PROBS.begin();
              flip_prob != FLIP_PROBS.end(); flip_prob++) {
        for (size_t iter = 0; iter < ITERS_PER_FLIP_PROB; iter++) {
            std::cout<<"Iteration "<<iter<<" for flip probability "<<*flip_prob
                     <<std::endl<<std::endl;

            random_deletions = std::unordered_set<Edge, EdgeHash>();
            random_additions = std::unordered_set<Edge, EdgeHash>();

            if (*flip_prob > 0.0) {
                // Flip edges randomly
                for (size_t i = 0; i < g.num_nodes(); i++) {
                    for (size_t j = i * (!DIRECTED); j < g.num_nodes(); j++) {
                        if (i == j and !has_self_loops) {
                            continue;
                        }
                        if (dist(gen) < *flip_prob) {
                            // Flip the edge!
                            if (g.has_edge(i, j)) {
                                random_deletions.insert(EDGE(i, j, DIRECTED));
                            }
                            else {
                                random_additions.insert(EDGE(i, j, DIRECTED));
                            }
                        }
                    }
                }
                std::cout<<std::endl<<"Adding: ";
                for (auto e = random_additions.begin();
                          e != random_additions.end(); e++) {
                    std::cout<<"("<<e->first<<", "<<e->second<<"), ";
                    g.add_edge(e->first, e->second);
                }
                std::cout<<std::endl;
                std::cout<<std::endl<<"Removing: ";
                for (auto e = random_deletions.begin();
                          e != random_deletions.end(); e++) {
                    std::cout<<"("<<e->first<<", "<<e->second<<"), ";
                    g.delete_edge(e->first, e->second);
                }
                std::cout<<std::endl;
            }

            auto result = genetic_alg_search(g, NUM_ITERATIONS,
                                             top_k, NUM_THREADS);

            // Report the results
            NTSparseGraph reporter = NTSparseGraph(g);

            int i = 0;
            for (auto result_itr = result.begin();
                      result_itr != result.end(); result_itr++) {

                // Make changes.
                for (auto edge_itr = result_itr->first.begin();
                          edge_itr != result_itr->first.end(); edge_itr++) {
                    reporter.flip_edge(edge_itr->first, edge_itr->second);
                }

                nt_results = traces(reporter, o);
                double log2_aut = std::log2l(nt_results.num_aut_base) +
                                 ((long double)(nt_results.num_aut_exponent)) *
                                                  std::log2l(10);

                std::cout<<"With a score of "<<result_itr->second<<" we have "
                         <<"log2(|Aut(G_H)|) of "<<log2_aut<<std::endl;

                std::cout<<"With edges: "<<std::endl;
                // Flip back.
                for (auto edge_itr = result_itr->first.begin();
                          edge_itr != result_itr->first.end(); edge_itr++) {
                    reporter.flip_edge(edge_itr->first, edge_itr->second);
                    std::cout<<"("<<edge_itr->first<<", "<<edge_itr->second
                             <<"), ";
                }
                std::cout<<std::endl<<std::endl;

                i++;
            }

            if (*flip_prob > 0.0) {
                // Restore graph
                for (auto e = random_additions.begin();
                          e != random_additions.end(); e++) {
                    g.delete_edge(e->first, e->second);
                }
                for (auto e = random_deletions.begin();
                          e != random_deletions.end(); e++) {
                    g.add_edge(e->first, e->second);
                }
            }
        }
    }

    /*
    size_t num_iterations = g.num_nodes() * g.num_nodes() *
                            g.num_nodes();
    num_iterations = (num_iterations < 1000000 ? 1000000 : num_iterations);
    */

    // auto result = simulated_annealing_search(g, num_iterations, top_k);


    return 0;
}
