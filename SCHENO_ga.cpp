#include<algorithm>
#include<cmath>
#include<iostream>
#include<random>
#include<unordered_set>
#include<stdexcept>
#include<string>
#include<vector>

#include "scheno/scheno.h"

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
    if (2 * sizeof(int) > sizeof(long)) {
        std::cout<<sizeof(int)<<std::endl;
        std::cout<<"<<< Warning! 2 * sizeof(int) > sizeof(long)"
                 <<std::endl<<"\tThis may cause problems if you use a node id "
                 <<">= "<<(long(1) << ((sizeof(int) * 8) / 2))<<" >>>"
                 <<std::endl;
    }

    std::vector<std::string> inputs = std::vector<std::string>();
    for (int i = 1; i < argc; i++) {
        inputs.push_back(argv[i]);
    }

    if (cmd_flag_present(inputs, "-h")) {
        std::cout<<"Guide to flags:"<<std::endl
                 <<"REQUIRED:"<<std::endl<<std::endl
                 <<"-graph <arg>:\tedgelist filename (including path)"
                 <<std::endl<<"\t\t\t* If a nodelist is also specified, "
                 <<"then all nodes"<<std::endl<<"\t\t\t\tin the nodelist will "
                 <<"be included."<<std::endl<<"\t\t\t* If a nodelist is NOT "
                 <<"specified, then the "<<std::endl<<"\t\t\t\tgraph's nodes "
                 <<"will be 0, 1, 2, ..., X, where "<<std::endl
                 <<"\t\t\t\tX is the max node label in the edgelist."
                 <<std::endl<<std::endl
                 <<"-o <arg>:\toutput filename base -- will append .txt "
                 <<"automatically"<<std::endl<<"\t\t\t* For example, if you "
                 <<"enter -o experiments/test_results/f"<<std::endl
                 <<"\t\t\t\tthen every iteration the top noise will be written"
                 <<std::endl<<"\t\t\t\tto experiments/test_results/f_noise.txt "
                 <<"and the top"<<std::endl<<"\t\t\t\tgraph to experiments/"
                 <<"test_results/f_graph.txt"
                 <<std::endl<<std::endl
                 <<std::endl<<"OPTIONAL:"<<std::endl<<std::endl
                 <<"-nodes <arg>:\tnodelist filename (including path)"
                 <<std::endl<<std::endl
                 <<"-d and -u:\tdirected and undirected respectively"
                 <<std::endl<<"\t\t\t* defaults to undirected"
                 <<std::endl<<std::endl
                 <<"-seed <arg>:\tan edgelist filename (including path)"
                 <<std::endl<<"\t\t\t* requires that -nodes be specified"
                 <<std::endl<<"\t\t\t* used as the initial noise set"
                 <<std::endl<<std::endl
                 <<"-legal_noise <arg>:\tan edgelist filename (including path)"
                 <<std::endl<<"\t\t\t* used to constrain which edges are "
                 <<"allowed to be part"<<std::endl<<"\t\t\t\tof the noise set."
                 <<std::endl<<"\t\t\t* uses the -nodes nodelist in the same way"
                 <<" as -graph"
                 <<std::endl<<"\t\t\t* requires that -nodes be specified"
                 <<std::endl<<std::endl
                 <<"-trials <arg>:\tnumber of times to do a full trial"
                 <<std::endl<<"\t\t\t* defaults to 1"
                 <<std::endl<<std::endl
                 <<"-topk <arg>:\tnumber of best results to return"
                 <<std::endl<<"\t\t\t* defaults to 10"
                 <<std::endl<<std::endl
                 <<"-gdepth <arg>:\tgene depth -- must be >= 1 -- defaults to 1"
                 <<std::endl<<std::endl
                 <<"-n_itr <arg>:\tnumber of iterations for genetic algorithm"
                 <<" -- must be >= 1"<<std::endl<<"\t\t\t* defaults to 80"
                 <<std::endl<<std::endl
                 <<"-score_heuristic:\ttells genetic algorithm to use a"
                 <<std::endl<<"\t\t\t\theuristic score as a tiebreaker"
                 <<std::endl<<std::endl
                 <<"-no_sample_heuristic:\ttells genetic algorithm to NOT use"
                 <<std::endl<<"\t\t\t\tits edge sampling heuristic"
                 <<std::endl<<std::endl
                 <<"-noise- <arg>:\tfraction of input edges to randomly remove"
                 <<std::endl<<"\t\t\t* for example, with 0.05, randomly removes "
                 <<"5% of edges"<<std::endl<<"\t\t\t* defaults to 0"
                 <<std::endl<<"\t\t\t* If -legal_noise is specified, the "
                 <<"removed edges will be"<<std::endl<<"\t\t\t\ta subset of"
                 <<"-legal_noise."
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
                 <<std::endl<<"\t\t\t* If -legal_noise is specified, the "
                 <<"added edges will be"<<std::endl<<"\t\t\t\ta subset of"
                 <<"-legal_noise."
                 <<std::endl<<std::endl
                 // <<"-use_g_as_ses:\t'Use G as Special Edge Set' rather than"
                 // <<std::endl<<"\t\t\t\tthe edges added/remove by noise+/-"
                 // <<std::endl<<std::endl
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
                 <<"-approx_iso:\tif set, uses a fast, inexact automorphisms"
                 <<" algorithm"
                 <<std::endl<<"\t\t\t* can actually be SLOWER on some graphs"
                 <<std::endl<<"\t\t\t\t(but sometimes is much faster)"
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

    size_t gene_depth = 1;
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

    bool score_heuristic = cmd_flag_present(inputs, "-score_heuristic");
    bool sample_heuristic = !cmd_flag_present(inputs, "-no_sample_heuristic");

    std::string edgelist_file;
    edgelist_file = get_cmd_option(inputs, "-graph");

    std::string output_file;
    output_file = get_cmd_option(inputs, "-o");

    std::string nodelist_file = "";  // No nodelist
    if (cmd_flag_present(inputs, "-nodes")) {
        nodelist_file = get_cmd_option(inputs, "-nodes");
    }

    if (nodelist_file.empty()) {
        if (cmd_flag_present(inputs, "-seed") ||
                cmd_flag_present(inputs, "-legal_noise")) {
            throw std::invalid_argument(std::string("Error! Cannot use -seed") +
                                        " or -legal_noise without -nodes");
        }
    }

    std::string seed_file = "";  // Start with empty noise set
    if (cmd_flag_present(inputs, "-seed")) {
        seed_file = get_cmd_option(inputs, "-seed");
    }

    std::string legal_noise_file = "";  // No limits on noise edge choice
    if (cmd_flag_present(inputs, "-legal_noise")) {
        legal_noise_file = get_cmd_option(inputs, "-legal_noise");
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

    // bool use_g_as_ses = cmd_flag_present(inputs, "-use_g_as_ses");

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

    bool full_iso = !cmd_flag_present(inputs, "-approx_iso");

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

    if (g.num_loops() > 0) {
        std::cout<<"<<< The current code cannot handle self-loops. "
                 <<"Ignoring all. >>>"<<std::endl<<std::endl;
        for (size_t i = 0; i < g.num_nodes(); i++) {
            if (g.has_edge(i, i)) {
                g.delete_edge(i, i);
            }
        }
    }

    SparseGraph seed_noise(directed, g.num_nodes());
    bool has_seed = !seed_file.empty();
    if (has_seed) {
        seed_noise = read_graph(directed, nodelist_file, seed_file);
    }
    if (seed_noise.num_loops() > 0) {
        std::cout<<"<<< The current code cannot handle self-loops. "
                 <<"Ignoring all in seed noise. >>>"<<std::endl<<std::endl;
        for (size_t i = 0; i < seed_noise.num_nodes(); i++) {
            if (seed_noise.has_edge(i, i)) {
                seed_noise.delete_edge(i, i);
            }
        }
    }

    SparseGraph legal_noise(directed, g.num_nodes());
    bool has_legal_noise = !legal_noise_file.empty();
    if (has_legal_noise) {
        legal_noise = read_graph(directed, nodelist_file, legal_noise_file);

        if (has_seed) {
            for (size_t i = 0; i < seed_noise.num_nodes(); i++) {
                for (auto nbr = seed_noise.out_neighbors(i).begin();
                          nbr != seed_noise.out_neighbors(i).end(); nbr++) {
                    if (!legal_noise.has_edge(i, *nbr)) {
throw std::invalid_argument(std::string("Error! Seed graph has edge (") + 
                            std::to_string(i) + ", " + std::to_string(*nbr) +
                            ") which is not listed in -legal_noise");
                    }
                }
            }
        }
    }

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

    if (max_flip_or_edge < g.num_nodes() * 2) {
        max_flip_or_edge = g.num_nodes() * 2;
    }
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

    double og_noise_plus = noise_plus;
    double og_noise_minus = noise_minus;

    if (has_extra_noise && has_legal_noise) {
        double expected_edge_dels = noise_minus * g.num_edges();
        double expected_edge_adds =
            noise_plus * (max_possible_edges - g.num_edges());

        double legal_adds = 0.0;
        double legal_dels = 0.0;
        for (size_t a = 0; a < legal_noise.num_nodes(); a++) {
            for (auto b_itr = legal_noise.out_neighbors(a).begin();
                      b_itr != legal_noise.out_neighbors(a).end(); b_itr++) {
                if (!directed && *b_itr < int(a)) {
                    continue;
                }
                if (g.has_edge(a, *b_itr)) {
                    legal_dels++;
                } else {
                    legal_adds++;
                }
            }
        }

        if (legal_dels < expected_edge_dels) {
            throw std::logic_error(
std::string("Error! Due to -noise-, supposed to delete roughly ") +
std::to_string(expected_edge_dels) + " edges,\tbut only have " +
std::to_string(legal_dels) + " options due to -legal_noise");
        }
        if (legal_adds < expected_edge_adds) {
            throw std::logic_error(
std::string("Error! Due to -noise+, supposed to add roughly ") +
std::to_string(expected_edge_adds) + " edges,\nbut only have " +
std::to_string(legal_adds) + " options due to -legal_noise");
        }

        noise_plus =  expected_edge_adds / legal_adds;
        noise_minus = expected_edge_dels / legal_dels;
    }

    /*
    // Used to get a score for the "special edge set"
    //  The special edge set is usually the noise added by noise+ and noise-,
    //  but if -use_g_as_ses is specified, then the special edge set is g.
    std::unordered_set<Edge, EdgeHash> se_del;
    std::unordered_set<Edge, EdgeHash> se_add;
    */

    for (size_t trial = 0; trial < trials; trial++) {
        std::cout<<"Trial "<<trial<<" for (noise-, noise+) = ("
                 <<og_noise_minus<<", "<<og_noise_plus<<")"
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
                    if (has_legal_noise && !legal_noise.has_edge(i, j)) {
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
            // log_probs = {-1.0, -1.0, -1.0, -1.0};
            NTSparseGraph g_nt = NTSparseGraph(g);
            log_probs = default_log2_noise_probs(g_nt, comb_util);
        }

        /*
        if (use_g_as_ses) {
            se_add.clear();
            se_del.clear();
            for (size_t a = 0; a < g.num_nodes(); a++) {
                for (auto b_itr = g.out_neighbors(a).begin();
                          b_itr != g.out_neighbors(a).end(); b_itr++) {
                    se_add.insert(EDGE((int) a, *b_itr, g.directed));
                }
            }
        } else {
            se_del = random_deletions;
            se_add = random_additions;
        }
        */

        auto result = genetic_alg_search(g, num_iterations,
                                         top_k, nt,
                                         gene_depth,
                                         log_probs,
                                         max_change_factor,
                                         score_heuristic,
                                         sample_heuristic,
                                         seed_noise,
                                         legal_noise,
                                         output_file,
                                         full_iso);

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

    return 0;
}
