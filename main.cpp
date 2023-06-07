#include<algorithm>
#include<cmath>
#include<iostream>
#include<random>
#include<unordered_set>
#include<stdexcept>
#include<string>
#include<vector>

#include "edge.h"
#include "edge_sampler.h"
#include "file_utils.h"
#include "genetic_alg_search.h"
#include "nauty_traces.h"
#include "noise_prob_choice.h"
#include "nt_sparse_graph.h"
#include "simulated_annealing_search.h"
#include "sparse_graph.h"

std::string get_cmd_option(const std::vector<std::string>& inputs,
                           const std::string& option) {
    for (size_t i = 0; i < inputs.size(); i++) {
        if (inputs[i] == option && i < inputs.size() - 1) {
            return inputs[i + 1];
        }
    }
    throw std::invalid_argument("Error! Expected argument for flag " + option);
}

bool cmd_flag_present(const std::vector<std::string>& inputs,
                      const std::string& option) {
    for (size_t i = 0; i < inputs.size(); i++) {
        if (inputs[i] == option) {
            return true;
        }
    }
    return false;
}

int main(int argc, char* argv[]) {
    std::vector<std::string> inputs = std::vector<std::string>();
    for (int i = 1; i < argc; i++) {
        inputs.push_back(argv[i]);
    }

    if (cmd_flag_present(inputs, "-h")) {
        std::cout<<"Guide to flags:"<<std::endl
                 <<"REQUIRED:"<<std::endl<<std::endl
                 <<"-graph <arg>:\tedgelist filename (including path)"
                 <<std::endl<<std::endl
                 <<std::endl<<"OPTIONAL:"<<std::endl<<std::endl
                 <<"-nodes <arg>:\tnodelist filename (including path)"
                 <<std::endl<<"\t\t\t* use if you want to ensure that nodes "
                 <<"with no edges are"<<std::endl<<"\t\t\t\tincluded"
                 <<std::endl<<std::endl
                 <<"-d and -u:\tdirected and undirected respectively"
                 <<std::endl<<"\t\t\t* defaults to undirected"
                 <<std::endl<<std::endl
                 <<"-trials <arg>:\tnumber of times to do a full trial"
                 <<std::endl<<"\t\t\t* defaults to 1"
                 <<std::endl<<std::endl
                 <<"-topk <arg>:\tnumber of best results to return"
                 <<std::endl<<"\t\t\t* defaults to 10"
                 <<std::endl<<std::endl
                 <<"-gdepth <arg>:\tgene depth -- must be >= 1 -- defaults to 2"
                 <<std::endl<<std::endl
                 <<"-n_itr <arg>:\tnumber of iterations for genetic algorithm"
                 <<" -- must be >= 1"<<std::endl<<"\t\t\t* defaults to 80"
                 <<std::endl<<std::endl
                 <<"-noise- <arg>:\tfraction of input edges to randomly remove"
                 <<std::endl<<"\t\t\t* for example, with 0.05, randomly removes "
                 <<"5% of edges"<<std::endl<<"\t\t\t* defaults to 0"
                 <<std::endl<<std::endl
                 <<"-noise+ <arg>:\tfraction of input non-edges to randomly add"
                 <<std::endl<<"\t\t\t* IMPORTANTLY, this is scaled to the"
                 <<" number of EDGES"<<std::endl<<"\t\t\t\tin the graph, NOT"
                 <<" the number of non-edges. Thus,"<<std::endl
                 <<"\t\t\t\tfor example, a value of 0.05 means adding"
                 <<std::endl<<"\t\t\t\troughly (0.05 * num_edges) non-edges."
                 <<std::endl<<"\t\t\t\tThis means the value can be > 1. A value"
                 <<" of 2, "<<std::endl<<"\t\t\t\tfor example, means adding "
                 <<"roughly (2 * num_edges) "<<std::endl<<"\t\t\t\tnon-edges."
                 <<std::endl<<"\t\t\t* defaults to 0"
                 <<std::endl<<std::endl
                 <<"-mcf <arg>:\taka 'max change factor' -- the candidate noise"
                 <<std::endl<<"\t\t\t\tset will have at most (mcf * num_edges)"
                 <<" elements"<<std::endl<<"\t\t\t* defaults to 1.0"
                 <<std::endl<<std::endl
                 <<"-p+ <arg>:\tUsed in the scoring function as prob. edges"
                 <<" were added"<<std::endl<<"\t\t\t* default calculated based"
                 <<" on input graph"
                 <<std::endl<<"\t\t\t* if used, must also use p-"
                 <<std::endl<<std::endl
                 <<"-p- <arg>:\tUsed in the scoring function as prob. edges"
                 <<" were removed"<<std::endl<<"\t\t\t* default calculated"
                 <<" based on input graph"
                 <<std::endl<<"\t\t\t* if used, must also use p+"
                 <<std::endl<<std::endl
                 <<"-nt <arg>:\tnumber of threads -- 0 means all available"
                 <<std::endl<<"\t\t\t* defaults to 0."
                 <<std::endl<<std::endl
                 ;

        return 0;
    }

    bool directed;
    if (cmd_flag_present(inputs, "-d")) {
        if (cmd_flag_present(inputs, "-u")) {
            throw std::invalid_argument("Error! Cannot pass both -u and -d");
        }
        directed = true;
    } else if (cmd_flag_present(inputs, "-u")) {
        directed = false;
    } else {
        // std::cout<<"No (un)directed flag passed. Assuming undirected."
        //          <<std::endl;
        directed = false;
    }

    size_t nt = 0;
    if (cmd_flag_present(inputs, "-nt")) {
        nt = std::stoul(get_cmd_option(inputs, "-nt"));
    }

    size_t trials = 1;
    if (cmd_flag_present(inputs, "-trials")) {
        trials = std::stoul(get_cmd_option(inputs, "-trials"));
        if (trials == 0) {
            throw std::invalid_argument("Error! trials must be > 0");
        }
    }

    size_t top_k = 10;  // Number of candidate noise sets to keep.
    if (cmd_flag_present(inputs, "-topk")) {
        top_k = std::stoul(get_cmd_option(inputs, "-topk"));
        if (top_k == 0) {
            throw std::invalid_argument("Error! topk must be > 0");
        }
    }

    size_t gene_depth = 2;
    if (cmd_flag_present(inputs, "-gdepth")) {
        gene_depth = std::stoul(get_cmd_option(inputs, "-gdepth"));
        if (gene_depth == 0) {
            throw std::invalid_argument("Error! gdepth must be > 0");
        }
    }

    size_t num_iterations = 80;
    if (cmd_flag_present(inputs, "-n_itr")) {
        num_iterations = std::stoul(get_cmd_option(inputs, "-n_itr"));
        if (num_iterations == 0) {
            throw std::invalid_argument("Error! n_itr must be > 0");
        }
    }

    std::string edgelist_file;
    edgelist_file = get_cmd_option(inputs, "-graph");

    std::string nodelist_file = "";  // No nodelist
    if (cmd_flag_present(inputs, "-nodes")) {
        nodelist_file = get_cmd_option(inputs, "-nodes");
    }

    double noise_minus = 0.0;
    if (cmd_flag_present(inputs, "-noise-")) {
        noise_minus = std::stod(get_cmd_option(inputs, "-noise-"));
        if (noise_minus < 0.0) {
            throw std::invalid_argument("Error! noise- must be >= 0");
        } else if (noise_minus >= 1.0) {
            throw std::invalid_argument("Error! noise- must be < 1.0");
        }
    }

    double noise_plus = 0.0;
    if (cmd_flag_present(inputs, "-noise+")) {
        noise_plus = std::stod(get_cmd_option(inputs, "-noise+"));
        if (noise_plus < 0.0) {
            throw std::invalid_argument("Error! noise+ must be >= 0");
        } // Check for an overly-large noise+ later
    }

    float max_change_factor = 1.0;
    if (cmd_flag_present(inputs, "-mcf")) {
        max_change_factor = std::stof(get_cmd_option(inputs, "-mcf"));
        if (max_change_factor <= 0.0) {
            throw std::invalid_argument("Error! mcf must be > 0");
        }
    }

    bool custom_p_values = false;

    long double p_plus = 0.0;
    if (cmd_flag_present(inputs, "-p+")) {
        custom_p_values = true;
        p_plus = std::stold(get_cmd_option(inputs, "-p+"));
        if (p_plus <= 0.0) {
            throw std::invalid_argument("Error! p+ must be > 0");
        }
    }

    long double p_minus = 0.0;
    if (cmd_flag_present(inputs, "-p-")) {
        custom_p_values = true;
        p_minus = std::stold(get_cmd_option(inputs, "-p-"));
        if (p_minus <= 0.0) {
            throw std::invalid_argument("Error! p- must be > 0");
        }
    }

    if (custom_p_values && (p_plus == 0.0 || p_minus == 0.0)) {
        throw std::invalid_argument(
                "Error! Cannot provide p+ without p- and vice versa");
    }

    // const bool corrupt_original = false;
    // Only used when corrupt_original is true.
    // const size_t num_additions = 0;
    // Only used when corrupt_original is true.
    // const size_t num_removals = 0;

    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;

    NautyTracesResults nt_results;

    std::random_device rd;  // Will provide a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    std::uniform_real_distribution<float> dist(0, 1);

    // 0  test_01 -- a 5-chain
    // 1  test_02 -- a 10-chain
    // 2  test_03 -- a 15 node binary tree
    // 3  test_04 -- a 4x5 grid
    // 4  test_05 -- a 6x7 grid
    // 5  test_06 -- a 31 node binary tree
    // 6  binary_tree_127 -- a 127 node binary tree
    // 7  johnson_10_3_120 -- a 120 10-choose-3 johnson graph
    // 8  ring_128 -- a 128 node ring
    // 9  wreath_d7_128 -- a wreath of 128 nodes each with degree 7

    std::vector<std::string> fake_nodes_names =
        {"test_01_nodes.txt",         "test_02_nodes.txt",
         "test_03_nodes.txt",         "test_04_nodes.txt",
         "test_05_nodes.txt",         "test_06_nodes.txt",
         "binary_tree_127_nodes.txt", "johnson_10_3_120_nodes.txt",
         "ring_128_nodes.txt",        "wreath_d7_128_nodes.txt"
        };
    std::vector<std::string> fake_edges_names =
        {"test_01_edges.txt", "test_02_edges.txt",
         "test_03_edges.txt", "test_04_edges.txt",
         "test_05_edges.txt", "test_06_edges.txt",
         "binary_tree_127_edges.txt", "johnson_10_3_120_edges.txt",
         "ring_128_edges.txt",        "wreath_d7_128_edges.txt"
        };
    std::vector<std::string> real_nodes_names =
        {"",                               "",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "jazz_collab_nodes.txt",
         "jazz_collab_nodes.txt",          "",
         "",                               "",
         "season_0_directed_nodes.txt",    "season_0_undirected_nodes.txt",
         "season_1_directed_nodes.txt",    "season_1_undirected_nodes.txt",
         "season_2_directed_nodes.txt",    "season_2_undirected_nodes.txt",
         "season_3_directed_nodes.txt",    "season_3_undirected_nodes.txt",
         "season_4_directed_nodes.txt",    "season_4_undirected_nodes.txt",
         "season_5_directed_nodes.txt",    "season_5_undirected_nodes.txt",
         "season_6_directed_nodes.txt",    "season_6_undirected_nodes.txt",
         "season_7_directed_nodes.txt",    "season_7_undirected_nodes.txt",
         "season_8_directed_nodes.txt",    "season_8_undirected_nodes.txt",
         "season_9_directed_nodes.txt",    "season_9_undirected_nodes.txt",
         "season_10_directed_nodes.txt",   "season_10_undirected_nodes.txt"
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
         "roget_thesaurus.g",              "pol_blogs.g",
         "season_0_directed_edges.txt",    "season_0_undirected_edges.txt",
         "season_1_directed_edges.txt",    "season_1_undirected_edges.txt",
         "season_2_directed_edges.txt",    "season_2_undirected_edges.txt",
         "season_3_directed_edges.txt",    "season_3_undirected_edges.txt",
         "season_4_directed_edges.txt",    "season_4_undirected_edges.txt",
         "season_5_directed_edges.txt",    "season_5_undirected_edges.txt",
         "season_6_directed_edges.txt",    "season_6_undirected_edges.txt",
         "season_7_directed_edges.txt",    "season_7_undirected_edges.txt",
         "season_8_directed_edges.txt",    "season_8_undirected_edges.txt",
         "season_9_directed_edges.txt",    "season_9_undirected_edges.txt",
         "season_10_directed_edges.txt",   "season_10_undirected_edges.txt"
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

    std::cout<<"Loading graph from files:"<<std::endl
             <<"    "<<nodelist_file<<std::endl
             <<"    "<<edgelist_file<<std::endl;
    SparseGraph g(directed);
    if (nodelist_file.empty()) {
        g = read_graph(directed, edgelist_file);
    } else {
        g = read_graph(directed, nodelist_file, edgelist_file);
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
    std::vector<Edge> random_deletions_sorted, random_additions_sorted;
    std::vector<Edge> flip_vec;

    size_t max_possible_edges =
            (g.num_nodes() * (g.num_nodes() - 1)) / (1 + size_t(!directed)) +
            (g.num_nodes() * size_t(g.num_loops() > 0));

    size_t max_edit_factor = 7;
    size_t max_flip_or_edge = (g.num_edges() < (max_possible_edges / 2) ?
                          g.num_edges() : (max_possible_edges - g.num_edges()));
    max_flip_or_edge *= max_edit_factor;
    CombinatoricUtility comb_util(max_possible_edges, max_flip_or_edge);

    // For genetic alg

    std::cout<<"Running for "<<num_iterations<<" iterations per trial..."
             <<std::endl;

    bool has_extra_noise = noise_plus > 0.0 || noise_minus > 0.0;

    double expected_flipped_edges, new_noise_plus;

    expected_flipped_edges = (noise_plus) * (double) g.num_edges();
    new_noise_plus = expected_flipped_edges /
                        ((double) max_possible_edges - g.num_edges());
    if (new_noise_plus >= 1.0) {
        throw std::invalid_argument(std::string("Error! The value for noise+") +
                " is too large. It requires flipping on average more non-edges"+
                " than are available to be flipped.");
    }
    noise_plus = new_noise_plus;

    for (size_t trial = 0; trial < trials; trial++) {
        std::cout<<"Trial "<<trial<<" for (noise-, noise+) = ("
                 <<noise_minus<<", "<<noise_plus<<")"
                 <<std::endl<<std::endl;
        random_deletions = std::unordered_set<Edge, EdgeHash>();
        random_additions = std::unordered_set<Edge, EdgeHash>();

        if (has_extra_noise) {
            // Flip edges randomly
            for (size_t i = 0; i < g.num_nodes(); i++) {
                for (size_t j = i * (!directed); j < g.num_nodes(); j++) {
                    if (i == j and !has_self_loops) {
                        continue;
                    }
                    if (g.has_edge(i, j)) {
                        if (dist(gen) < noise_minus) {
                            random_deletions.insert(EDGE(i, j, directed));
                        }
                    }
                    else {
                        if (dist(gen) < noise_plus) {
                            random_additions.insert(EDGE(i, j, directed));
                        }
                    }
                }
            }
            random_additions_sorted= std::vector<Edge>(random_additions.begin(),
                                                       random_additions.end());
            random_deletions_sorted= std::vector<Edge>(random_deletions.begin(),
                                                       random_deletions.end());
            std::sort(random_additions_sorted.begin(),
                      random_additions_sorted.end());
            std::sort(random_deletions_sorted.begin(),
                      random_deletions_sorted.end());

            std::cout<<std::endl<<"Adding: ";
            for (auto e = random_additions_sorted.begin();
                      e != random_additions_sorted.end(); e++) {
                std::cout<<"("<<e->first<<", "<<e->second<<"), ";
                g.add_edge(e->first, e->second);
            }
            std::cout<<std::endl;
            std::cout<<std::endl<<"Removing: ";
            for (auto e = random_deletions_sorted.begin();
                      e != random_deletions_sorted.end(); e++) {
                std::cout<<"("<<e->first<<", "<<e->second<<"), ";
                g.delete_edge(e->first, e->second);
            }
            std::cout<<std::endl;
        }

        std::vector<long double> log_probs;
        if (custom_p_values) {
            log_probs = {std::log2l(p_plus),
                         std::log2l(1.0 - p_plus),
                         std::log2l(p_minus),
                         std::log2l(1.0 - p_minus)};
        } else {
            log_probs = {-1.0, -1.0, -1.0, -1.0};
            // NTSparseGraph g_nt = NTSparseGraph(g);
            // log_probs = default_log2_noise_probs(g_nt, comb_util);
        }

        auto result = genetic_alg_search(g, num_iterations,
                                         top_k, nt,
                                         gene_depth,
                                         random_deletions,
                                         random_additions,
                                         log_probs,
                                         max_change_factor);

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

            flip_vec = std::vector<Edge>();
            // Flip back.
            for (auto edge_itr = result_itr->first.begin();
                      edge_itr != result_itr->first.end(); edge_itr++) {
                flip_vec.push_back(*edge_itr);
                reporter.flip_edge(edge_itr->first, edge_itr->second);
            }

            std::sort(flip_vec.begin(), flip_vec.end());

            std::cout<<"With edges: "<<std::endl;
            for (auto edge_itr = flip_vec.begin();
                      edge_itr != flip_vec.end(); edge_itr++) {
                std::cout<<"("<<edge_itr->first<<", "<<edge_itr->second
                         <<"), ";
            }
            std::cout<<std::endl<<std::endl;

            i++;
        }

        if (has_extra_noise) {
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

    /*
    size_t num_iterations = g.num_nodes() * g.num_nodes() *
                            g.num_nodes();
    num_iterations = (num_iterations < 1000000 ? 1000000 : num_iterations);
    */

    // auto result = simulated_annealing_search(g, num_iterations, top_k);


    return 0;
}
