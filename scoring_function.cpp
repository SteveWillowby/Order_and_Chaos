#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"

#include<cmath>
#include<unordered_set>
#include<vector>

double score(double log2_g_aut,
             NTSparseGraph& g,
             const Coloring<int>& node_orbit_coloring,
             const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
             Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
             const std::unordered_set<Edge,EdgeHash>& edge_additions,
             const std::unordered_set<Edge,EdgeHash>& edge_removals) {

    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;

    NautyTracesResults nt_results;

    size_t n = g.num_nodes();
    size_t m = g.num_edges();
    size_t max_num_edges = (n * (n - 1)) / (int(!g.directed) + 1);

    size_t num_additions = edge_additions.size();
    size_t num_removals = edge_removals.size();

    size_t m_prime = m + num_additions - num_removals;

    double log2_stabilizer_size, log2_hypothesis_aut;

    // Edge additions will be colored with a new color (max prev color + 1).
    int addition_color = *(edge_orbit_coloring.colors.rend()) + 1;
    // Edge deletions will be colored with a color that preserves their orbit
    //  info: new color = old color + deletion_color_base.
    int deletion_color_base = addition_color + 1;

    // Perform the edge additions in actuality.
    for (auto edge_itr = edge_additions.begin();
              edge_itr != edge_additions.end(); edge_itr++) {
        g.add_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.set(*edge_itr, addition_color);
    }

    // Perform the edge deletions in color only.
    for (auto edge_itr = edge_deletions.begin();
              edge_itr != edge_deletions.end(); edge_itr++) {
        editable_edge_orbit_coloring.set(*edge_itr, deletion_color_base +
                                            edge_orbit_coloring[*edge_itr]);
    }

    // Get the size of the stabilizer set for the changes.
    //
    // Remember that the stabilizer size is the same both in g and in the
    //  hypothesis graph.
    NTPartition stabilizer_coloring =
                    g.nauty_traces_coloring(node_orbit_coloring,
                                            editable_edge_orbit_coloring);
    nt_results = traces(g, o, stabilizer_coloring);
    log2_stabilizer_size = std::log2(nt_results.num_aut_base) +
                           double(nt_results.num_aut_exponent) +
                           __combinatoric_utility.log2(10);

    // Perform the edge deletions in actuality.
    for (auto edge_itr = edge_deletions.begin();
              edge_itr != edge_deletions.end(); edge_itr++) {
        g.remove_edge(edge_itr->first, edge_itr->second);
    }

    // Get the raw auto orbit size for the hypothesis graph.
    nt_results = traces(g, o);
    log2_hypothesis_aut = std::log2(nt_results.num_aut_base) +
                          double(nt_results.num_aut_exponent) +
                          __combinatoric_utility.log2(10);


    // Restore the deleted edges.
    for (auto edge_itr = edge_deletions.begin();
              edge_itr != edge_deletions.end(); edge_itr++) {
        g.add_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.set(*edge_itr,
                                         edge_orbit_coloring[*edge_itr]);
    }

    // Un-color and un-add the added edges.
    for (auto edge_itr = edge_additions.begin();
              edge_itr != edge_additions.end(); edge_itr++) {
        g.remove_edge(edge_itr->first, edge_itr->second);
    }

    // Perform the probability calculations.
    return  2.0 * log2_hypothesis_aut -
             (log2_stabilizer_size +
              __combinatoric_utility.log2_a_choose_b(m_prime, num_additions) +
              __combinatoric_utility.log2_a_choose_b(max_num_edges - m_prime,
                                                     num_removals));
}


__CombinatoricUtility::__CombinatoricUtility() {
    // For convenience we say that 0! = 1, and thus log2(0!) = 0.
    log2_factorials = std::vector<double>(2, 0.0);
    // There is no log2 for 0, but oh well. Put a placeholder value.
    log2_s = std::vector<double>(2, 0.0);
}

double __CombinatoricUtility::log2(size_t x) {
    if (x >= log2_s.size()) {
        for (size_t i = log2_s.size(); i <= x; i++) {
            log2_s.push_back(std::log2(i));
        }
    }
    return log2_s[x];
}

double __CombinatoricUtility::log2_factorial(size_t x) {
    if (x < log2_factorials.size()) {
        return log2_factorials[x];
    }
    double v = log2_factorials[log2_factorials.size() - 1];
    while (x >= log2_factorials.size()) {
        v += log2(log2_factorials.size());
        log2_factorials.push_back(v);
    }
    return log2_factorials[x];
}

double __CombinatoricUtility::log2_a_choose_b(size_t a, size_t b) {
    double v1 = log2_factorial(a);
    return v1 - (log2_factorials[b] + log2_factorials[a - b]);
}
