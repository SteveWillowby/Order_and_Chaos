#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"

#include<cmath>
#include<stdexcept>
#include<string>
#include<unordered_set>
#include<vector>

long double score(NTSparseGraph& g,
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
    size_t max_num_edges = (n * (n - 1)) / (size_t(!g.directed) + 1);

    size_t num_additions = edge_additions.size();
    size_t num_removals = edge_removals.size();

    size_t m_prime = m + num_additions - num_removals;

    long double log2_stabilizer_size, log2_hypothesis_aut;

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
    log2_stabilizer_size = std::log2l(nt_results.num_aut_base) +
                           (long double)(nt_results.num_aut_exponent) +
                           __combinatoric_utility.log2(10);

    // Perform the edge deletions in actuality.
    for (auto edge_itr = edge_deletions.begin();
              edge_itr != edge_deletions.end(); edge_itr++) {
        g.remove_edge(edge_itr->first, edge_itr->second);
    }

    // Get the raw auto orbit size for the hypothesis graph.
    nt_results = traces(g, o);
    log2_hypothesis_aut = std::log2l(nt_results.num_aut_base) +
                          (long double)(nt_results.num_aut_exponent) +
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
    log2_factorials = std::vector<long double>(2, 0.0);
    // There is no log2 for 0, but oh well. Put a placeholder value.
    log2_s = std::vector<long double>(2, 0.0);
    max_access = 1;
}

void __CombinatoricUtility::set_max_access(size_t max_e, size_t max_f) {
    size_t edge_flip_end_1 = max_f + 1;
    size_t edge_flip_start_2 = max_e + 1 - (max_f + 1);
    size_t edge_flip_end_2 = max_e + 1;

    log2_s[0] = std::vector<long double>(max_f + 1, 0.0);
    log2_s[1] = std::vector<long double>(max_f + 1, 0.0);
    log2_factorials[0] = std::vector<long double>(max_f + 1, 0.0);
    log2_factorials[1] = std::vector<long double>(max_f + 1, 0.0);

    for (size_t i = 2; i < edge_flip_end_1; i++) {
        log2_s[0][i] = std::log2l(i);
    }
    for (size_t i = edge_flip_start_2; i < edge_flip_end_2; i++) {
        log2_s[1][i - edge_flip_start_2] = std::log2l(i);
    }

    // This is the Kahan, Babushka, and Klein sum algorithm.
    //  Copied from the pseudocode on Wikipedia (12/8/22).
    //
    // A., Klein (2006). "A generalized Kahan-Babuska-Summation-Algorithm"
    //
    //  The point of the algorithm is to do floating-point summation while
    //   minimizing floating-point rounding errors.

    long double sum = 0.0;
    long double log_value;
    long double cs = 0.0;
    long double ccs = 0.0;
    long double c = 0.0;
    long double cc = 0.0;
    volatile long double t = 0.0;  // volatile ensures that the compiler doesn't
    volatile long double z = 0.0;  //   optimize these operation orders away

    for (size_t i = 2; i < edge_flip_end_2; i++) {
        // Access memoized versions if available.
        if (i < edge_flip_end_1) {
            log_value = log2_s[0][i];
        } else if (i >= edge_flip_start_2) {
            log_value = log2_s[1][i - edge_flip_start_2];
        } else {
            log_value = std::log2l(i);
        }

        t = sum + log_value;
        // fabsl is just the long double abs() function
        if (std::fabsl(sum) >= std::fabsl(log_value)) {
            z = sum - t;
            c = z + log_value;
        } else {
            z = log_value - t;
            c = z + sum;
        }
        sum = t;
        t = cs + c;
        if (std::fabsl(cs) >= std::fabsl(c)) {
            z = cs - t;
            cc = z + c;
        } else {
            z = c - t;
            cc = z + cs;
        }
        cs = t;
        ccs = ccs + cc;

        // Store memoized versions when relevant.
        if (i < edge_flip_end_1) {
            log2_factorials[0][i] = sum + cs + ccs;
        } else if (i >= edge_flip_start_2) {
            log2_factorials[1][i - edge_flip_start_2] = sum + cs + ccs;
        }
    }

}

long double __CombinatoricUtility::log2(size_t x) {
    if (x < edge_flip_end_1) {
        return log2_s[0][x];
    }
    if (x < edge_flip_start_2 || x >= edge_flip_end_2) {
        throw std::range_error(
                std::string("Error! Each thread must have called ")
                + "__combinatoric_utility.set_max_access(max_e, max_f) with "
                + "max_e and max_f chosen according to the criterion listed "
                + "in scoring_function.h"
                + "  NOTE: This will take O(max_e) time and O(max_f) space.");
    }
    return log2_s[1][x - edge_flip_start_2];
}

long double __CombinatoricUtility::log2_factorial(size_t x) {
    if (x < edge_flip_end_1) {
        return log2_factorials[0][x];
    }
    if (x < edge_flip_start_2 || x >= edge_flip_end_2) {
        throw std::range_error(
                std::string("Error! Each thread must have called ")
                + "__combinatoric_utility.set_max_access(max_e, max_f) with "
                + "max_e and max_f chosen according to the criterion listed "
                + "in scoring_function.h"
                + "  NOTE: This will take O(max_e) time and O(max_f) space.");
    }
    return log2_factorials[1][x - edge_flip_start_2];
}

long double __CombinatoricUtility::log2_a_choose_b(size_t a, size_t b) {
    if (b > a) {
        throw std::invalid_argument("Error! Cannot do a-choose-b with b > a.");
    }
    return log2_factorial(a) - (log2_factorial(b) + log2_factorial(a - b));
}
