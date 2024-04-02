#include "wl_measures.h"

#include "nt_wrappers/nauty_traces.h"

#include "scoring_function.h"

#include<array>
#include<cmath>
#include<stdexcept>
#include<string>
#include<unordered_set>
#include<vector>

long double score(NTSparseGraph& g, const CombinatoricUtility& comb_util,
                  const Coloring<int>& node_orbit_coloring,
                  const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
                  Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals,
                  const long double log2_p_plus, const long double log2_p_minus,
                  const long double log2_1_minus_p_plus,
                  const long double log2_1_minus_p_minus,
                  const size_t max_change,
                  bool full_iso) {

    if (edge_additions.size() + edge_removals.size() > max_change) {
        return -INFINITY;  // INFINITY is defined in cmath
    }

    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;

    NautyTracesResults nt_results;

    // size_t n = g.num_nodes();
    // size_t m = g.num_edges();
    // size_t max_num_edges = (n * (n - 1)) / (size_t(!g.directed) + 1);

    size_t num_additions = edge_additions.size();
    size_t num_removals = edge_removals.size();

    // size_t m_prime = m + num_additions - num_removals;

    long double log2_stabilizer_size, log2_hypothesis_aut;

    // Edge additions will be colored with a new color (max prev color + 1).
    int addition_color = *(edge_orbit_coloring.colors().rend()) + 1;
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
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
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
    if (full_iso) {
        nt_results = traces(g, o, stabilizer_coloring);
    } else {
        nt_results = fake_iso(g, o, stabilizer_coloring);
    }
    log2_stabilizer_size = std::log2l(nt_results.num_aut_base) +
                           ((long double)(nt_results.num_aut_exponent)) *
                                            comb_util.log2(10);

    // Perform the edge deletions in actuality.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
    }

    // Get the raw auto orbit size for the hypothesis graph.
    if (full_iso) {
        nt_results = traces(g, o);
    } else {
        nt_results = fake_iso(g, o);
    }
    log2_hypothesis_aut = std::log2l(nt_results.num_aut_base) +
                          ((long double)(nt_results.num_aut_exponent)) *
                                           comb_util.log2(10);

    // Restore the deleted edges.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.add_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.set(*edge_itr,
                                         edge_orbit_coloring[*edge_itr]);
    }

    // Un-color and un-add the added edges.
    for (auto edge_itr = edge_additions.begin();
              edge_itr != edge_additions.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.erase(*edge_itr);
    }

    // Calculate the probability that this exact noise set would be chosen.

    // Note that we ultimately want to calculate the probability that a noise
    //  set of this SIZE would be chosen. However, that includes some
    //  (a choose b) terms that cancel out with other aspects of our scoring,
    //  formula, so we can do this simpler calculation instead.

    // FURTHER note that we return a value which is RELATIVE to what we would
    //  have gotten if the noise size had been zero. This reduces floating
    //  point errors.
    long double log2_sequence;
    if (log2_p_plus == -1.0 && log2_p_minus == -1.0) {
        log2_sequence = 0;
    } else {
        // Absolute Calculation
        // log2_sequence = num_additions * log2_p_minus +
        //             (m_prime - num_additions) * log2_1_minus_p_minus +
        //             num_removals * log2_p_plus +
        //        ((max_num_edges - m_prime) - num_removals) * log2_1_minus_p_plus;

        // Relative Calculation
        log2_sequence = num_additions * log2_p_minus +
                         -(num_additions * log2_1_minus_p_minus) +
                        num_removals * log2_p_plus +
                         -(num_removals * log2_1_minus_p_plus);
    }

    // Perform the probability calculations.
    return ((2.0 * log2_hypothesis_aut) - log2_stabilizer_size) + log2_sequence;
}


std::array<long double, 6>
  score_breakdown(SparseGraph& g, const CombinatoricUtility& comb_util,
                  const Coloring<int>& node_orbit_coloring,
                  const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
                  Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals,
                  const long double log2_p_plus, const long double log2_p_minus,
                  const long double log2_1_minus_p_plus,
                  const long double log2_1_minus_p_minus,
                  const size_t max_change,
                  bool full_iso) {

    if (edge_additions.size() + edge_removals.size() > max_change) {
        // INFINITY is defined in cmath
        return std::array<long double, 6>({-INFINITY, 0, 0, 0, 0, 0});
    }

    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;

    NautyTracesResults nt_results;

    // size_t n = g.num_nodes();
    // size_t m = g.num_edges();
    // size_t max_num_edges = (n * (n - 1)) / (size_t(!g.directed) + 1);

    size_t num_additions = edge_additions.size();
    size_t num_removals = edge_removals.size();

    // size_t m_prime = m + num_additions - num_removals;

    long double log2_stabilizer_size, log2_hypothesis_aut;

    long double log2_hypothesis_singleton_aut;
    long double log2_combined_singleton_aut;
    size_t num_singletons;

    // Edge additions will be colored with a new color (max prev color + 1).
    int addition_color = *(edge_orbit_coloring.colors().rend()) + 1;
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
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        editable_edge_orbit_coloring.set(*edge_itr, deletion_color_base +
                                            edge_orbit_coloring[*edge_itr]);
    }

    NTSparseGraph g_nt(g);

    // Get the size of the stabilizer set for the changes.
    //
    // Remember that the stabilizer size is the same both in g and in the
    //  hypothesis graph.
    NTPartition stabilizer_coloring =
                    g_nt.nauty_traces_coloring(node_orbit_coloring,
                                               editable_edge_orbit_coloring);
    if (full_iso) {
        nt_results = traces(g_nt, o, stabilizer_coloring);
    } else {
        nt_results = fake_iso(g_nt, o, stabilizer_coloring);
    }
    log2_stabilizer_size = std::log2l(nt_results.num_aut_base) +
                           ((long double)(nt_results.num_aut_exponent)) *
                                            comb_util.log2(10);

    num_singletons = 0;
    for (size_t n = 0; n < g.num_nodes(); n++) {
        if (g.neighbors(n).size() == 0) {
            num_singletons++;
        }
    }
    log2_combined_singleton_aut = comb_util.log2_factorial(num_singletons);

    // Perform the edge deletions in actuality.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
    }

    // Get the raw auto orbit size for the hypothesis graph.
    g_nt = NTSparseGraph(g);
    if (full_iso) {
        nt_results = traces(g_nt, o);
    } else {
        nt_results = fake_iso(g_nt, o);
    }
    log2_hypothesis_aut = std::log2l(nt_results.num_aut_base) +
                          ((long double)(nt_results.num_aut_exponent)) *
                                           comb_util.log2(10);

    num_singletons = 0;
    for (size_t n = 0; n < g.num_nodes(); n++) {
        if (g.neighbors(n).size() == 0) {
            num_singletons++;
        }
    }
    log2_hypothesis_singleton_aut = comb_util.log2_factorial(num_singletons);

    // Restore the deleted edges.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.add_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.set(*edge_itr,
                                         edge_orbit_coloring[*edge_itr]);
    }

    // Un-color and un-add the added edges.
    for (auto edge_itr = edge_additions.begin();
              edge_itr != edge_additions.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.erase(*edge_itr);
    }

    // Calculate the probability that this exact noise set would be chosen.

    // Note that we ultimately want to calculate the probability that a noise
    //  set of this SIZE would be chosen. However, that includes some
    //  (a choose b) terms that cancel out with other aspects of our scoring,
    //  formula, so we can do this simpler calculation instead.

    // FURTHER note that we return a value which is RELATIVE to what we would
    //  have gotten if the noise size had been zero. This reduces floating
    //  point errors.
    long double log2_sequence;
    if (log2_p_plus == -1.0 && log2_p_minus == -1.0) {
        log2_sequence = 0;
    } else {
        // Absolute Calculation
        // log2_sequence = num_additions * log2_p_minus +
        //             (m_prime - num_additions) * log2_1_minus_p_minus +
        //             num_removals * log2_p_plus +
        //        ((max_num_edges - m_prime) - num_removals) * log2_1_minus_p_plus;

        // Relative Calculation
        log2_sequence = num_additions * log2_p_minus +
                         -(num_additions * log2_1_minus_p_minus) +
                        num_removals * log2_p_plus +
                         -(num_removals * log2_1_minus_p_plus);
    }

    // long double log2_stabilizer_size, log2_hypothesis_aut;
    // long double log2_hypothesis_singleton_aut;
    // long double log2_combined_singleton_aut;

    long double log2_orbit_size = log2_hypothesis_aut - log2_stabilizer_size;
    long double log2_orbit_singleton_contribution =
            log2_hypothesis_singleton_aut - log2_combined_singleton_aut;

    long double total_score =
            log2_hypothesis_aut + log2_orbit_size + log2_sequence;

    return std::array<long double, 6>({
                       total_score,
                       log2_hypothesis_aut - log2_hypothesis_singleton_aut,
                       log2_hypothesis_singleton_aut,
                       log2_orbit_size - log2_orbit_singleton_contribution,
                       log2_orbit_singleton_contribution,
                       log2_sequence});
}

std::pair<long double, long double>
            score(NTSparseGraph& g, const CombinatoricUtility& comb_util,
                  const Coloring<int>& node_orbit_coloring,
                  const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
                  Coloring<Edge,EdgeHash>& editable_edge_orbit_coloring,
                  const std::unordered_set<Edge,EdgeHash>& edge_additions,
                  const std::unordered_set<Edge,EdgeHash>& edge_removals,
                  const long double log2_p_plus, const long double log2_p_minus,
                  const long double log2_1_minus_p_plus,
                  const long double log2_1_minus_p_minus,
                  const size_t max_change,
                  const double* precomputed_wl_diff,
                  double* cost_matrix,
                  ptrdiff_t* col_for_row, ptrdiff_t* row_for_col,
                  double* u, double* v, void* workspace,
                  double* pw_scores_1, double* pw_scores_2,
                  size_t* start_indices,
                  bool full_iso) {

    if (edge_additions.size() + edge_removals.size() > max_change) {
        // INFINITY is defined in cmath
        return std::pair<long double, long double>(-INFINITY, -INFINITY);
    }

    NautyTracesOptions o;
    o.get_node_orbits = false;
    o.get_edge_orbits = false;
    o.get_canonical_node_order = false;

    NautyTracesResults nt_results;

    // size_t n = g.num_nodes();
    // size_t m = g.num_edges();
    // size_t max_num_edges = (n * (n - 1)) / (size_t(!g.directed) + 1);

    size_t num_additions = edge_additions.size();
    size_t num_removals = edge_removals.size();

    // size_t m_prime = m + num_additions - num_removals;

    long double log2_stabilizer_size, log2_hypothesis_aut;

    // Edge additions will be colored with a new color (max prev color + 1).
    int addition_color = *(edge_orbit_coloring.colors().rend()) + 1;
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
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
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
    if (full_iso) {
        nt_results = traces(g, o, stabilizer_coloring);
    } else {
        nt_results = fake_iso(g, o, stabilizer_coloring);
    }
    log2_stabilizer_size = std::log2l(nt_results.num_aut_base) +
                           ((long double)(nt_results.num_aut_exponent)) *
                                            comb_util.log2(10);

    long double heuristic_log2_stabilizer_size =
                      wl_symmetry_measure(g, NULL,
                                          &editable_edge_orbit_coloring,
                                          precomputed_wl_diff,
                                          cost_matrix,
                                          col_for_row, row_for_col, u, v,
                                          workspace, pw_scores_1, pw_scores_2,
                                          start_indices);

    // Perform the edge deletions in actuality.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
    }

    // Get the raw auto orbit size for the hypothesis graph.
    if (full_iso) {
        nt_results = traces(g, o);
    } else {
        nt_results = fake_iso(g, o);
    }
    log2_hypothesis_aut = std::log2l(nt_results.num_aut_base) +
                          ((long double)(nt_results.num_aut_exponent)) *
                                           comb_util.log2(10);

    long double heuristic_log2_symmetry =
                      wl_symmetry_measure(g, NULL, NULL, NULL, cost_matrix,
                                          col_for_row, row_for_col, u, v,
                                          workspace, pw_scores_1, pw_scores_2,
                                          start_indices);

    // Restore the deleted edges.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.add_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.set(*edge_itr,
                                         edge_orbit_coloring[*edge_itr]);
    }

    // Un-color and un-add the added edges.
    for (auto edge_itr = edge_additions.begin();
              edge_itr != edge_additions.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.erase(*edge_itr);
    }

    // Calculate the probability that this exact noise set would be chosen.

    // Note that we ultimately want to calculate the probability that a noise
    //  set of this SIZE would be chosen. However, that includes some
    //  (a choose b) terms that cancel out with other aspects of our scoring,
    //  formula, so we can do this simpler calculation instead.

    // FURTHER note that we return a value which is RELATIVE to what we would
    //  have gotten if the noise size had been zero. This reduces floating
    //  point errors.
    long double log2_sequence;
    if (log2_p_plus == -1.0 && log2_p_minus == -1.0) {
        log2_sequence = 0;
    } else {
        // Absolute Calculation
        // log2_sequence = num_additions * log2_p_minus +
        //             (m_prime - num_additions) * log2_1_minus_p_minus +
        //             num_removals * log2_p_plus +
        //        ((max_num_edges - m_prime) - num_removals) * log2_1_minus_p_plus;

        // Relative Calculation
        log2_sequence = num_additions * log2_p_minus +
                         -(num_additions * log2_1_minus_p_minus) +
                        num_removals * log2_p_plus +
                         -(num_removals * log2_1_minus_p_plus);
    }

    // Perform the probability calculations.
    return std::pair<long double, long double>(
           ((2.0 * log2_hypothesis_aut) - log2_stabilizer_size) + log2_sequence,
           // heuristic_log2_symmetry);
            ((2.0 * heuristic_log2_symmetry) - heuristic_log2_stabilizer_size));
}

// Outdated
long double score(NTSparseGraph& g, const CombinatoricUtility& comb_util,
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
    int addition_color = *(edge_orbit_coloring.colors().rend()) + 1;
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
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
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
                           ((long double)(nt_results.num_aut_exponent)) *
                                            comb_util.log2(10);

    // Perform the edge deletions in actuality.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
    }

    // Get the raw auto orbit size for the hypothesis graph.
    nt_results = traces(g, o);
    log2_hypothesis_aut = std::log2l(nt_results.num_aut_base) +
                          ((long double)(nt_results.num_aut_exponent)) *
                                           comb_util.log2(10);


    // Restore the deleted edges.
    for (auto edge_itr = edge_removals.begin();
              edge_itr != edge_removals.end(); edge_itr++) {
        g.add_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.set(*edge_itr,
                                         edge_orbit_coloring[*edge_itr]);
    }

    // Un-color and un-add the added edges.
    for (auto edge_itr = edge_additions.begin();
              edge_itr != edge_additions.end(); edge_itr++) {
        g.delete_edge(edge_itr->first, edge_itr->second);
        editable_edge_orbit_coloring.erase(*edge_itr);
    }

    return  2.0 * log2_hypothesis_aut -
             (log2_stabilizer_size +
              comb_util.log2_a_choose_b(m_prime, num_additions) +
              comb_util.log2_a_choose_b(max_num_edges - m_prime,
                                        num_removals));
}

CombinatoricUtility::CombinatoricUtility(size_t max_e, size_t max_f) {
    update_max_access(max_e, max_f);
}

void CombinatoricUtility::update_max_access(size_t max_e, size_t max_f) {
    // These checks are here because we access log2(10) regardless of graph size
    if (max_e < 10) {
        max_e = 10;
    }
    if (max_f < 10) {
        max_f = 10;
    }

    // Regular assignments.
    edge_flip_end_1 = max_f + 1;
    edge_flip_start_2 = max_e + 1 - (max_f + 1);
    edge_flip_end_2 = max_e + 1;

    if (max_f * 2 >= max_e) {
        // All indices from 0 to max_e (inclusive) are to be used.
        edge_flip_end_1 = max_e + 1;
        edge_flip_start_2 = max_e + 1;
        edge_flip_end_2 = max_e + 1;
    }

    size_t second_segment_size = edge_flip_end_2 - edge_flip_start_2;

    log2_s[0] = std::vector<long double>(edge_flip_end_1, 0.0);
    log2_factorials[0] = std::vector<long double>(edge_flip_end_1, 0.0);
    if (second_segment_size > 0) {
        log2_s[1] = std::vector<long double>(second_segment_size, 0.0);
        log2_factorials[1] = std::vector<long double>(second_segment_size, 0.0);
    }

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
        if (std::fabs(sum) >= std::fabs(log_value)) {
            z = sum - t;
            c = z + log_value;
        } else {
            z = log_value - t;
            c = z + sum;
        }
        sum = t;
        t = cs + c;
        if (std::fabs(cs) >= std::fabs(c)) {
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

long double CombinatoricUtility::log2(size_t x) const {
    if (x < edge_flip_end_1) {
        return log2_s[0][x];
    }
    if (x < edge_flip_start_2 || x >= edge_flip_end_2) {
        throw std::range_error(
                std::string("Error! CombinatoricUtility object - You must ")
                + "have initialized the object or subsequently called "
                + "update_max_access() with max_e and max_f chosen according "
                + "to the criterion listed in scoring_function.h"
                + "  NOTE: This will take O(max_e) time and O(max_f) space.");
    }
    return log2_s[1][x - edge_flip_start_2];
}

long double CombinatoricUtility::log2_factorial(size_t x) const {
    if (x < edge_flip_end_1) {
        return log2_factorials[0][x];
    }
    if (x < edge_flip_start_2 || x >= edge_flip_end_2) {
        throw std::range_error(
                std::string("Error! CombinatoricUtility object - You must ")
                + "have initialized the object or subsequently called "
                + "update_max_access() with max_e and max_f chosen according "
                + "to the criterion listed in scoring_function.h"
                + "  NOTE: This will take O(max_e) time and O(max_f) space.");
    }
    return log2_factorials[1][x - edge_flip_start_2];
}

long double CombinatoricUtility::log2_a_choose_b(size_t a, size_t b) const {
    if (b > a) {
        throw std::invalid_argument("Error! Cannot do a-choose-b with b > a.");
    }
    return log2_factorial(a) - (log2_factorial(b) + log2_factorial(a - b));
}
