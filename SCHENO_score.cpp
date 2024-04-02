#include<algorithm>
#include<array>
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
                 <<"-edges <arg>:\tedgelist filename (including path)"
                 <<std::endl<<"\t\t\t* These edges are the noise to be scored."
                 <<std::endl<<std::endl
                 <<std::endl<<"OPTIONAL:"<<std::endl<<std::endl
                 <<"-nodes <arg>:\tnodelist filename (including path)"
                 <<std::endl<<std::endl
                 <<"-d and -u:\tdirected and undirected respectively"
                 <<std::endl<<"\t\t\t* defaults to undirected"
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
                 <<"-no_extra:\tReduces number of automorphism calls and"
                 <<" amount of non-score info"
                 <<std::endl<<std::endl
                 <<"-approx:\tUse an \"approximate\" isomorphism algorithm to" 
                 <<" calculate answer."
                 <<std::endl<<"\t\t\t* Scores may be too high but are often exactly correct."
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
        directed = false;
    }

    std::string graph_file;
    graph_file = get_cmd_option(inputs, "-graph");

    std::string edgelist_file;
    edgelist_file = get_cmd_option(inputs, "-edges");

    std::string nodelist_file = "";  // No nodelist
    if (cmd_flag_present(inputs, "-nodes")) {
        nodelist_file = get_cmd_option(inputs, "-nodes");
    } else {
        nodelist_file = "/tmp/__score_only_nodes.txt";
        make_nodelist(graph_file, nodelist_file, false);
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

    bool no_extra = cmd_flag_present(inputs, "-no_extra");
    bool full_iso = !cmd_flag_present(inputs, "-approx");

    NautyTracesOptions o;
    o.get_node_orbits = true;
    o.get_edge_orbits = true;
    o.get_canonical_node_order = false;

    NautyTracesResults nt_results;

    std::cout<<"Loading graph from files:"<<std::endl
             <<"    "<<nodelist_file<<std::endl
             <<"    "<<graph_file<<std::endl;
    SparseGraph g(directed);
    if (nodelist_file.empty()) {
        g = read_graph(directed, graph_file);
        std::cout<<g.num_nodes()<<std::endl;
    } else {
        g = read_graph(directed, nodelist_file, graph_file);
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

    SparseGraph noise_candidate(directed, g.num_nodes());
    noise_candidate = read_graph(directed, nodelist_file,
                                 edgelist_file);

    bool has_self_loops = g.num_loops() > 0;

    double log2_aut;
    NTSparseGraph g_info = NTSparseGraph(g);

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

    std::vector<long double> log_probs;
    if (custom_p_values) {
        log_probs = {std::log2l(p_plus),
                     std::log2l(1.0 - p_plus),
                     std::log2l(p_minus),
                     std::log2l(1.0 - p_minus)};
    } else {
        log_probs = default_log2_noise_probs(g_info, comb_util);
    }

    if (full_iso) {
        std::cout<<"Running `traces` on original graph."<<std::endl;
        nt_results = traces(g_info, o);
    } else {
        std::cout<<"Running `fake_iso` on original graph."<<std::endl;
        nt_results = fake_iso(g_info, o);
    }
    log2_aut = std::log2l(nt_results.num_aut_base) +
                 ((long double)(nt_results.num_aut_exponent)) *
                                  std::log2l(10);
    std::cout<<"The original graph has log2_aut = "<<log2_aut<<std::endl;

    std::unordered_set<Edge, EdgeHash> deletions, additions;
    std::vector<Edge> deletions_sorted, additions_sorted;
    std::vector<Edge> flip_vec;

    deletions = std::unordered_set<Edge, EdgeHash>();
    additions = std::unordered_set<Edge, EdgeHash>();

    for (size_t i = 0; i < g.num_nodes(); i++) {
        for (size_t j = i * (!directed); j < g.num_nodes(); j++) {
            if (i == j and !has_self_loops) {
                continue;
            }
            if (g.has_edge(i, j) && noise_candidate.has_edge(i, j)) {
                deletions.insert(EDGE(i, j, directed));
            }
            else if (noise_candidate.has_edge(i, j)) {
                additions.insert(EDGE(i, j, directed));
            }
        }
    }

    Coloring<Edge, EdgeHash> edge_coloring;
    for (auto x =  nt_results.edge_orbits.colors().begin();
              x != nt_results.edge_orbits.colors().end(); x++) {
        for (auto y =  nt_results.edge_orbits.cell(*x).begin();
                  y != nt_results.edge_orbits.cell(*x).end(); y++) {
            edge_coloring.set(*y, *x);
        }
    }

    std::array<long double, 6> score_info =
          score_breakdown(g, comb_util,
                          nt_results.node_orbits,
                          nt_results.edge_orbits,
                          edge_coloring,
                          additions, deletions,
                          log_probs[0], log_probs[2],
                          log_probs[1], log_probs[3],
                          max_flip_or_edge,
                          full_iso);

    if (!no_extra) {
        for (auto x = additions.begin(); x != additions.end(); x++) {
            g.flip_edge(x->first, x->second);
        }
        for (auto x = deletions.begin(); x != deletions.end(); x++) {
            g.flip_edge(x->first, x->second);
        }

        g_info = NTSparseGraph(g);

        if (full_iso) {
            nt_results = traces(g_info, o);
        } else {
            nt_results = fake_iso(g_info, o);
        }
        log2_aut   = std::log2l(nt_results.num_aut_base) +
                         ((long double)(nt_results.num_aut_exponent)) *
                                          std::log2l(10);

        std::cout<<std::endl;
        std::cout<<"The modified graph has log2_aut = "<<log2_aut<<std::endl;
    }

    std::cout<<std::endl<<"Score Breakdown:"<<std::endl;
    std::cout<<"Log(Aut from connected components)"<<std::endl<<score_info[1]
             <<std::endl;
    std::cout<<"Log(Aut from singletons)"<<std::endl<<score_info[2]<<std::endl;
    std::cout<<"Log(AO ignoring singleton swaps)"<<std::endl<<score_info[3]
             <<std::endl;
    std::cout<<"Log(AO due to singleton swaps)"<<std::endl<<score_info[4]
             <<std::endl;
    std::cout<<"Log(Noise size prob)"<<std::endl<<score_info[5]<<std::endl;
    std::cout<<"Total score:"<<std::endl<<score_info[0]<<std::endl;

    return 0;
}
