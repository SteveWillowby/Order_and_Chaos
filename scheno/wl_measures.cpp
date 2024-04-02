#include<cmath>
#include<cstring>
#include<map>
#include<stdexcept>
#include<unordered_set>
#include<unordered_map>
#include<utility>
#include<vector>

#include "nt_wrappers/nauty_traces.h"

#include "Jonker_Volgenant/src/assignAlgs2D.h"
#include "wl_measures.h"

// Returns the absolute value of the difference between x and y
inline double SCHENO__abs_diff(double x, double y) {
    return (x - y >= 0.0 ? x - y : -(x - y));
}

inline double SCHENO__max(double x, double y) {
    return (x > y ? x : y);
}

// Adds a pair to the set of sets (really a map of sets)
void SCHENO__pair_insert(std::unordered_map<int, std::unordered_set<int>>& pairs,
                      int x, int y) {
    auto x_itr = pairs.find(x);
    if (x_itr == pairs.end()) {
        pairs.insert(std::pair<int, std::unordered_set<int>>(
                        x, std::unordered_set<int>({y})));
    } else {
        x_itr->second.insert(y);
    }
}

void initialize_wl_diff_matrix(const Graph& g,
                               const Coloring<int>* node_coloring,
                               const Coloring<Edge, EdgeHash>* edge_coloring,
                               double* the_matrix,
                               const size_t* start_indices);

// This is called by initialize_wl_diff_matrix
void impose_coloring_constraints(const Graph& g,
                                 const Coloring<int>* node_coloring,
                                 const Coloring<Edge, EdgeHash>* edge_coloring,
                                 double* the_matrix,
                                 const size_t* start_indices);

// Perform a single iteration of WL difference computation.
//  Track which pairs still need to be calculated.
void perform_wl_diff_iteration(const Graph& g,
                               // node_coloring should be implicit in
                               //   source_pairs_to_calc
                               const Coloring<Edge, EdgeHash>* edge_coloring,
                               double* cost_matrix,
                               ptrdiff_t* col_for_row, ptrdiff_t* row_for_col,
                               double* u, double* v, void* workspace,
                               const double* source_scores, double* target_scores,
                               const size_t* start_indices,
                               std::unordered_map<int, std::unordered_set<int>>&
                                        source_pairs_to_calc,
                               std::unordered_map<int, std::unordered_set<int>>&
                                        target_pairs_to_calc);


double* wl_node_differences(const Graph& g,
                            const Coloring<int>* node_coloring,
                            const Coloring<Edge, EdgeHash>* edge_coloring,
                            double* cost_matrix,
                            ptrdiff_t* col_for_row, ptrdiff_t* row_for_col,
                            double* u, double* v, void* workspace,
                            double* pw_scores_1, double* pw_scores_2,
                            const size_t* start_indices) {

    size_t n = g.num_nodes();
    size_t matrix_size = ((n * (n - 1)) / 2) + n;

    // Initialize both matrices.
    initialize_wl_diff_matrix(g, node_coloring, edge_coloring,
                              pw_scores_1, start_indices);
    std::memcpy(pw_scores_2, pw_scores_1, matrix_size * sizeof(double));

    double* pairwise_scores[2];
    pairwise_scores[0] = pw_scores_1;
    pairwise_scores[1] = pw_scores_2;

    size_t prev_mat_idx = 0;
    size_t next_mat_idx = 1;

    // Keep track of which pairs might not have converged.
    std::unordered_map<int, std::unordered_set<int>> pairs_to_calc;
    std::unordered_map<int, std::unordered_set<int>> next_pairs_to_calc;

    int a, b;
    for (a = 0; a < (int) n; a++) {
        for (b = a + 1; b < (int) n; b++) {
            if (pw_scores_1[start_indices[a] + (b - a)] == 1.0 ||
                   g.neighbors(a).size() == 0 || g.neighbors(b).size() == 0) {
                // Distance to node is always 1.0 (or 0.0). Don't recalculate.
                continue;
            }
            SCHENO__pair_insert(pairs_to_calc, a, b);
        }
    }

    while (pairs_to_calc.size() > 0) {
        // iteration++;

        next_pairs_to_calc.clear();
        perform_wl_diff_iteration(g, edge_coloring,
                                  cost_matrix, col_for_row, row_for_col,
                                  u, v, workspace,
                                  pairwise_scores[prev_mat_idx],
                                  pairwise_scores[next_mat_idx],
                                  start_indices,
                                  pairs_to_calc, next_pairs_to_calc);


        prev_mat_idx = next_mat_idx;
        next_mat_idx = (-next_mat_idx) + 1;

        // Copy changes over lest they are not updated a second time and
        //  then accidentally revert to the old values.
        for (auto a_itr = pairs_to_calc.begin();
                  a_itr != pairs_to_calc.end(); a_itr++) {
            a = a_itr->first;
            for (auto b_itr = a_itr->second.begin();
                      b_itr != a_itr->second.end(); b_itr++) {
                b = *b_itr;
                pairwise_scores[next_mat_idx][start_indices[a] + (b - a)] =
                    pairwise_scores[prev_mat_idx][start_indices[a] + (b - a)];
            }
        }

        std::swap(pairs_to_calc, next_pairs_to_calc);
    }

    return pairwise_scores[prev_mat_idx];  // Return the final result

}

long double wl_symmetry_measure(const Graph& g,
                                const Coloring<int>* node_coloring,
                                const Coloring<Edge, EdgeHash>* edge_coloring,
                                const double* precomputed_basics,
                                double* cost_matrix,
                                ptrdiff_t* col_for_row,
                                ptrdiff_t* row_for_col,
                                double* u, double* v,
                                void* workspace,
                                double* pw_scores_1,
                                double* pw_scores_2,
                                const size_t* start_indices) {

    size_t n = g.num_nodes();
    size_t matrix_size = ((n * (n - 1)) / 2) + n;

    double *old_mat, *new_mat;

    if (precomputed_basics == NULL) {
        old_mat = wl_node_differences(g, node_coloring, edge_coloring,
                                      cost_matrix, col_for_row, row_for_col,
                                      u, v, workspace, pw_scores_1,
                                      pw_scores_2, start_indices);

        if (old_mat == pw_scores_1) {
            new_mat = pw_scores_2;
        } else {
            new_mat = pw_scores_1;
        }
    } else {
        std::memcpy(pw_scores_1, precomputed_basics,
                        matrix_size * sizeof(double));
        if (node_coloring != NULL || edge_coloring != NULL) {
            impose_coloring_constraints(g, node_coloring, edge_coloring,
                                        pw_scores_1, start_indices);

            old_mat = wl_node_differences(g, node_coloring, edge_coloring,
                                          cost_matrix,
                                          col_for_row, row_for_col,
                                          u, v, workspace, pw_scores_1,
                                          pw_scores_2, start_indices);

            if (old_mat == pw_scores_1) {
                new_mat = pw_scores_2;
            } else {
                new_mat = pw_scores_1;
            }
        } else {
            old_mat = pw_scores_1;
            new_mat = pw_scores_2;
        }
    }

    int a, b;
    std::memcpy(new_mat, old_mat, matrix_size * sizeof(double));

    // Quit when the largest orbit size is no larger than QUIT_FACTOR
    const long double QUIT_FACTOR = 1.5;

    long double value;
    long double log2_aut_approx = 0.0;
    long double max_orbit_size = 0.0;
    int node_with_max_orbit = 0;
    for (a = 0; a < int(n); a++) {
        value = wl_orbit_size(a, n, old_mat, start_indices);
        if (value > max_orbit_size) {
            max_orbit_size = value;
            node_with_max_orbit = a;
        }
    }

    std::unordered_map<int, std::unordered_set<int>> pairs_to_calc;
    std::unordered_map<int, std::unordered_set<int>> next_pairs_to_calc;

    size_t min, max;
    while (max_orbit_size > QUIT_FACTOR) {
        log2_aut_approx += std::log2(max_orbit_size);

        // Make node_with_max_orbit a singleton
        for (a = 0; a < int(n); a++) {
            if (a == node_with_max_orbit) {
                continue;
            } else if (a < node_with_max_orbit) {
                min = a;
                max = node_with_max_orbit;
            } else {
                min = node_with_max_orbit;
                max = a;
            }
            old_mat[start_indices[min] + (max - min)] = 1.0;
            new_mat[start_indices[min] + (max - min)] = 1.0;
        }

        // Figure out which pairs to recompute.
        pairs_to_calc.clear();
        for (auto nbr_itr = g.neighbors(node_with_max_orbit).begin();
                  nbr_itr != g.neighbors(node_with_max_orbit).end(); nbr_itr++){
            for (b = 0; b < int(n); b++) {
                if (b == int(node_with_max_orbit) || b == *nbr_itr ||
                        g.neighbors(b).size() == 0) {
                    continue;
                } else if (b < *nbr_itr) {
                    min = b;
                    max = *nbr_itr;
                } else {
                    min = *nbr_itr;
                    max = b;
                }
                if (old_mat[start_indices[min] + (max - min)] < 1.0) {
                    SCHENO__pair_insert(pairs_to_calc, min, max);
                }
            }
        }

        // Run until convergence
        while (pairs_to_calc.size() > 0) {
            next_pairs_to_calc.clear();
            perform_wl_diff_iteration(g, edge_coloring,
                                      cost_matrix, col_for_row, row_for_col,
                                      u, v, workspace,
                                      old_mat, new_mat,
                                      start_indices,
                                      pairs_to_calc, next_pairs_to_calc);

            // Copy all changes from new mat into old mat
            for (auto a_itr = pairs_to_calc.begin();
                      a_itr != pairs_to_calc.end(); a_itr++) {
                a = a_itr->first;
                for (auto b_itr = a_itr->second.begin();
                          b_itr != a_itr->second.end(); b_itr++) {
                    b = *b_itr;
                    old_mat[start_indices[a] + (b - a)] =
                        new_mat[start_indices[a] + (b - a)];
                }
            }

            std::swap(pairs_to_calc, next_pairs_to_calc);
        }

        // Collect results

        max_orbit_size = 0.0;
        for (a = 0; a < int(n); a++) {
            value = wl_orbit_size(a, n, old_mat, start_indices);
            if (value > max_orbit_size) {
                max_orbit_size = value;
                node_with_max_orbit = a;
            }
        }
    }

    return log2_aut_approx;
}

long double wl_orbit_size(const size_t node, const size_t num_nodes,
                          const double* pw_scores,
                          const size_t* start_indices) {
    // TODO: Experiment with this exponent value (or make it principled)
    const long double EXPONENT = 2.0;
    long double size = 0.0;
    long double similarity;
    size_t min, max;
    for (size_t a = 0; a < num_nodes; a++) {
        if (a == node) {
            size += 1.0;
            continue;
        } else if (a < node) {
            min = a;
            max = node;
        } else {
            min = node;
            max = a;
        }
        // Similarity is 1 - Difference
        similarity = 1.0 - pw_scores[start_indices[min] + (max - min)];
        if (similarity != 0.0) {
            size += std::pow(similarity, EXPONENT);
        }
    }
    return size;
}


void impose_coloring_constraints(const Graph& g,
                                 const Coloring<int>* node_coloring,
                                 const Coloring<Edge, EdgeHash>* edge_coloring,
                                 double* the_matrix,
                                 const size_t* start_indices) {

    std::unordered_map<size_t, size_t> color_counts_a =
                                        std::unordered_map<size_t, size_t>();
    std::unordered_map<size_t, size_t> color_counts_b =
                                        std::unordered_map<size_t, size_t>();

    bool can_map;
    size_t nbr;
    for (size_t a = 0; a < g.num_nodes(); a++) {
        for (size_t b = a; b < g.num_nodes(); b++) {
            if (a == b) {
                continue;
            }
            if (node_coloring != NULL &&
                    (*node_coloring)[a] != (*node_coloring)[b]) {
                the_matrix[start_indices[a] + (b - a)] = 1.0;
                continue;
            }

            // This chunk of code requires that if a has edges of color C
            //  and b has NO edges of color C, then the nodes cannot be aligned
            //  (and same with the vice versa case).
            if (edge_coloring != NULL) {
                // Reset color counts
                for (auto c_itr = edge_coloring->colors().begin();
                          c_itr != edge_coloring->colors().end(); c_itr++) {
                    color_counts_a[*c_itr] = 0;
                    color_counts_b[*c_itr] = 0;
                }
                // Collect color counts
                for (auto nbr_itr = g.out_neighbors(a).begin();
                          nbr_itr != g.out_neighbors(a).end(); nbr_itr++) {
                    nbr = *nbr_itr;
                    color_counts_a[(*edge_coloring)[EDGE(a,nbr, g.directed)]]++;
                }
                for (auto nbr_itr = g.out_neighbors(b).begin();
                          nbr_itr != g.out_neighbors(b).end(); nbr_itr++) {
                    nbr = *nbr_itr;
                    color_counts_b[(*edge_coloring)[EDGE(b,nbr, g.directed)]]++;
                }
                // Compare color existence
                can_map = true;
                for (auto ccnt_itr = color_counts_a.begin();
                          ccnt_itr != color_counts_a.end(); ccnt_itr++) {
                    if ((ccnt_itr->second == 0) !=
                            (color_counts_b[ccnt_itr->first] == 0)) {
                        can_map = false;
                        break;
                    }
                }
                if (g.directed && can_map) {
                    // Reset color counts
                    for (auto c_itr = edge_coloring->colors().begin();
                              c_itr != edge_coloring->colors().end(); c_itr++) {
                        color_counts_a[*c_itr] = 0;
                        color_counts_b[*c_itr] = 0;
                    }
                    // Collect color counts
                    for (auto nbr_itr = g.in_neighbors(a).begin();
                              nbr_itr != g.in_neighbors(a).end(); nbr_itr++) {
                        nbr = *nbr_itr;
                        color_counts_a[(*edge_coloring)[EDGE(nbr,a, g.directed)]]++;
                    }
                    for (auto nbr_itr = g.in_neighbors(b).begin();
                              nbr_itr != g.in_neighbors(b).end(); nbr_itr++) {
                        nbr = *nbr_itr;
                        color_counts_b[(*edge_coloring)[EDGE(nbr,b, g.directed)]]++;
                    }
                    // Compare color existence
                    for (auto ccnt_itr = color_counts_a.begin();
                              ccnt_itr != color_counts_a.end(); ccnt_itr++) {
                        if ((ccnt_itr->second == 0) !=
                                (color_counts_b[ccnt_itr->first] == 0)) {
                            can_map = false;
                            break;
                        }
                    }
                }

                if (!can_map) {
                    the_matrix[start_indices[a] + (b - a)] = 1.0;
                    continue;
                }
            }
        }
    }
}

void initialize_wl_diff_matrix(const Graph& g,
                               const Coloring<int>* node_coloring,
                               const Coloring<Edge, EdgeHash>* edge_coloring,
                               double* the_matrix,
                               const size_t* start_indices) {

    double aon, bon, ain, bin, o, i;
    for (size_t a = 0; a < g.num_nodes(); a++) {
        for (size_t b = a; b < g.num_nodes(); b++) {
            if (a == b) {
                the_matrix[start_indices[a]] = 0.0;
                continue;
            }

            // If g is undirected, out_neighbors = in_neighbors = neighbors
            aon = g.out_neighbors(a).size();
            bon = g.out_neighbors(b).size();
            ain =  g.in_neighbors(a).size();
            bin =  g.in_neighbors(b).size();
            if (aon == 0.0 || bon == 0.0) {
                if (aon == 0.0 && bon == 0.0) {
                    o = 0.0;  // Both singletons
                } else {
                    o = 1.0;  // One singleton
                }
            } else {
                // Both have edges
                o = SCHENO__abs_diff(aon, bon) / SCHENO__max(aon, bon);
            }
            if (ain == 0.0 || bin == 0.0) {
                if (ain  == 0.0 && bin == 0.0) {
                    i = 0.0;  // Both singletons
                } else {
                    i = 1.0;  // One singleton
                }
            } else {
                // Both have edges
                i = SCHENO__abs_diff(ain, bin) / SCHENO__max(ain, bin);
            }
            the_matrix[start_indices[a] + (b - a)] = (i + o) / 2.0;
        }
    }

    impose_coloring_constraints(g, node_coloring, edge_coloring, the_matrix,
                                start_indices);
}

void perform_wl_diff_iteration(const Graph& g,
                               // node_coloring should be implicit in
                               //   source_pairs_to_calc
                               const Coloring<Edge, EdgeHash>* edge_coloring,
                               double* cost_matrix,
                               ptrdiff_t* col_for_row, ptrdiff_t* row_for_col,
                               double* u, double* v, void* workspace,
                               const double* source_scores, double* target_scores,
                               const size_t* start_indices,
                               std::unordered_map<int, std::unordered_set<int>>&
                                        source_pairs_to_calc,
                               std::unordered_map<int, std::unordered_set<int>>&
                                        target_pairs_to_calc) {

    const double CONVERGENCE_FACTOR = 0.01;

    double prev_value, next_value, diff;
    int a, b, x, y, temp_int, a_ish, b_ish, x_ish, y_ish;
    size_t num_rows, num_cols, r_, c_, square_size;
    const std::unordered_set<int>* unmatched_a_ish;
    const std::unordered_set<int>* unmatched_b_ish;
    bool had_edges;

    // Fill in target_scores
    for (auto a_itr = source_pairs_to_calc.begin();
              a_itr != source_pairs_to_calc.end(); a_itr++) {
        a = a_itr->first;
        for (auto b_itr = a_itr->second.begin();
                  b_itr != a_itr->second.end(); b_itr++) {
            b = *b_itr;

            // Calculate difference between a and b
            next_value = 0;
            had_edges = true;
            for (size_t k = 0; k < 1 + size_t(g.directed); k++) {
                // Get copies of the neighbor sets
                //
                // a_ish will have <= neighbors than b_ish
                if (k == 0) {
                    if (g.out_neighbors(a).size() >
                            g.out_neighbors(b).size()) {
                        a_ish = b;
                        b_ish = a;
                    } else {
                        a_ish = a;
                        b_ish = b;
                    }
                    unmatched_a_ish = &g.out_neighbors(a_ish);
                    unmatched_b_ish = &g.out_neighbors(b_ish);
                } else {
                    // Go a second time for in_neighbors
                    if (g.in_neighbors(a).size() >
                            g.in_neighbors(b).size()) {
                        a_ish = b;
                        b_ish = a;
                    } else {
                        a_ish = a;
                        b_ish = b;
                    }
                    unmatched_a_ish = &g.in_neighbors(a_ish);
                    unmatched_b_ish = &g.in_neighbors(b_ish);
                }

                // Create the cost matrix
                num_rows = unmatched_a_ish->size();
                num_cols = unmatched_b_ish->size();
                if (num_rows == 0) {
                    if (!had_edges) {
                        throw std::logic_error(
"Error! Should not be calculating distance for a node with zero neighbors.");
                    }
                    had_edges = false;
                    continue;
                }
                square_size = num_cols;
                for (size_t i = 0; i < num_rows; i++) {
                    v[i] = 0.0;
                }
                for (size_t i = 0; i < num_cols; i++) {
                    u[i] = 0.0;
                }

                r_ = 0;
                for (auto x_itr = unmatched_a_ish->begin();
                          x_itr != unmatched_a_ish->end(); x_itr++) {
                    c_ = 0;
                    for (auto y_itr = unmatched_b_ish->begin();
                              y_itr != unmatched_b_ish->end(); y_itr++) {

                        x = *x_itr;
                        y = *y_itr;
                        x_ish = x;
                        y_ish = y;
                        if (x > y) {
                            x = y_ish;
                            y = x_ish;
                        }

                        diff = source_scores[start_indices[x] + (y - x)];

                        // See if edge colorings force matching (a, x) to (b, y)
                        //  to get a score of 1.0
                        if (edge_coloring != NULL && diff < 1.0 && ((k == 0 &&
                                (*edge_coloring)[EDGE(a_ish, x_ish, g.directed)] !=
                                (*edge_coloring)[EDGE(b_ish, y_ish, g.directed)]) ||
                                                                    (k == 1 &&
                                (*edge_coloring)[EDGE(x_ish, a_ish, g.directed)] !=
                                (*edge_coloring)[EDGE(y_ish, b_ish, g.directed)]))) {
                            diff = 1.0;
                        }
                        cost_matrix[r_ + square_size * c_] = diff;
                                
                        c_++;
                    }
                    r_++;
                }

                for (r_ = num_rows; r_ < square_size; r_++) {
                    for (c_ = 0; c_ < num_cols; c_++) {
                        cost_matrix[r_ + square_size * c_] = 1.0;
                    }
                }

                diff =
                    assign2DCBasic(cost_matrix, col_for_row, row_for_col,
                                   workspace, u, v, square_size, num_cols);

                if (diff == -1.0) {
                    throw std::logic_error(
                            "Assignment problem deemed infeasible.");
                }

                next_value += diff / double(square_size);

            }
            next_value /= (1.0 + double(g.directed));

            prev_value = source_scores[start_indices[a] + (b - a)];

            target_scores[start_indices[a] + (b - a)] = next_value;

            // See if the value changed by a factor of CONVERGENCE_FACTOR
            //
            // If so, we will recalculate for pairs (x, y) where x is a's
            //  neighbor and y is b's neighbor.
            if (SCHENO__abs_diff(next_value, prev_value) / next_value >
                    CONVERGENCE_FACTOR) {

                for (auto x_itr = g.neighbors(a).begin();
                          x_itr != g.neighbors(a).end(); x_itr++) {
                    for (auto y_itr = g.neighbors(b).begin();
                              y_itr != g.neighbors(b).end(); y_itr++) {
                        x = *x_itr;
                        y = *y_itr;
                        if (x == y) {
                            continue;
                        }
                        if (x > y) {
                            temp_int = x;
                            x = y;
                            y = temp_int;
                        }

                        if (target_scores[start_indices[x] + (y - x)] == 1.0) {
                            // Distance was set to max. Don't recalculate.
                            continue;
                        }

                        SCHENO__pair_insert(target_pairs_to_calc, x, y);
                    }
                }
            }
        }
    }
}
