#include<cmath>
#include<map>
#include<stdexcept>
#include<unordered_set>
#include<unordered_map>
#include<utility>
#include<vector>

#include "Jonker_Volgenant/src/assignAlgs2D.h"
#include "graph.h"

#include "wl_similarity.h"

double SYM__abs_diff(double x, double y) {
    double v = x - y;
    if (v < 0.0) {
        return -v;
    }
    return v;
}

void SYM__pair_insert(std::unordered_map<int, std::unordered_set<int>>& pairs,
                      int x, int y) {
    auto x_itr = pairs.find(x);
    if (x_itr == pairs.end()) {
        pairs.insert(std::pair<int, std::unordered_set<int>>(
                        x, std::unordered_set<int>({y})));
    } else {
        x_itr->second.insert(y);
    }
}


double* wl_similarity_measure(bool sum_result, long double* neg_sum,
                              const Graph& g, const size_t node,
                              double* cost_matrix,
                              ptrdiff_t* col_for_row, ptrdiff_t* row_for_col,
                              double* u, double* v, void* workspace,
                              double* pw_scores_1, double* pw_scores_2,
                              size_t* start_indices) {

    const double CONVERGENCE_FACTOR = 0.05;
    const size_t MAX_ITERATIONS = -1;  // Go until convergence
    size_t iteration = 0;

    size_t n = g.num_nodes();
    // size_t m = g.num_edges();
    // bool has_self_loops = g.num_loops() > 0;
    bool directed = g.directed;
    // For simplicity, always leave space for self-loops
    size_t matrix_size = ((n * (n - 1)) / 2) + n;

    double score_bound = ((double) n);
    double gap_cost =       score_bound / 2.0;
    double scaling_factor = gap_cost;
    if (directed) {
        // Sums can get twice as large, so divide by twice as much
        scaling_factor = gap_cost * 2.0;
    }

    double* pairwise_scores[2];
    pairwise_scores[0] = pw_scores_1;
    pairwise_scores[1] = pw_scores_2;

    if (directed) {
        for (size_t a = 0; a < n; a++) {
            for (size_t b = a; b < n; b++) {
                if (a == b) {
                    pairwise_scores[0][start_indices[a]] = 0.0;
                    pairwise_scores[1][start_indices[a]] = 0.0;
                    continue;
                }
                if (a == node || b == node) {
                    pairwise_scores[0][start_indices[a] + (b - a)] = n;
                    pairwise_scores[1][start_indices[a] + (b - a)] = n;
                    continue;
                }
                pairwise_scores[0][start_indices[a] + (b - a)] =
                    ((SYM__abs_diff((double) g.out_neighbors(a).size(),
                                    (double) g.out_neighbors(b).size()) +
                      SYM__abs_diff((double) g.in_neighbors(a).size(),
                                    (double) g.in_neighbors(b).size())) / 2.0);
            }
        }
    } else {
        for (size_t a = 0; a < n; a++) {
            for (size_t b = a; b < n; b++) {
                if (a == b) {
                    pairwise_scores[0][start_indices[a]] = 0.0;
                    pairwise_scores[1][start_indices[a]] = 0.0;
                    continue;
                }
                if (a == node || b == node) {
                    pairwise_scores[0][start_indices[a] + (b - a)] = n;
                    pairwise_scores[1][start_indices[a] + (b - a)] = n;
                    continue;
                }
                pairwise_scores[0][start_indices[a] + (b - a)] =
                    SYM__abs_diff((double) g.neighbors(a).size(),
                                  (double) g.neighbors(b).size());
            }
        }
    }

    size_t prev_mat_idx = 0;
    size_t next_mat_idx = 1;

    double prev_value, next_value, diff;
    int a, b, x, y, temp_int;
    size_t num_rows, num_cols, r_, c_, square_size;
    const std::unordered_set<int>* unmatched_1;
    const std::unordered_set<int>* unmatched_2;

    // Keep track of which pairs might not have converged.
    std::unordered_map<int, std::unordered_set<int>> pairs_to_calc;
    std::unordered_map<int, std::unordered_set<int>> next_pairs_to_calc;

    for (a = 0; a < (int) n; a++) {
        for (b = a + 1; b < (int) n; b++) {
            if ((size_t) a == node || (size_t) b == node) {
                // Distance to node is always n. Don't recalculate.
                continue;
            }
            SYM__pair_insert(pairs_to_calc, a, b);
        }
    }

    while (pairs_to_calc.size() > 0 && iteration < MAX_ITERATIONS) {
        iteration++;

        next_pairs_to_calc.clear();

        // Fill in pairwise_scores[next_mat_idx]
        for (auto a_itr = pairs_to_calc.begin();
                  a_itr != pairs_to_calc.end(); a_itr++) {
            a = a_itr->first;
            for (auto b_itr = a_itr->second.begin();
                      b_itr != a_itr->second.end(); b_itr++) {
                b = *b_itr;

                // Calculate difference between a and b
                next_value = 0;
                for (size_t k = 0; k < 1 + size_t(directed); k++) {
                    // Get copies of the neighbor sets
                    if (k == 0) {
                        if (g.out_neighbors(a).size() >
                                g.out_neighbors(b).size()) {
                            unmatched_2 = &g.out_neighbors(a);
                            unmatched_1 = &g.out_neighbors(b);
                        } else {
                            unmatched_1 = &g.out_neighbors(a);
                            unmatched_2 = &g.out_neighbors(b);
                        }
                    } else {
                        // Go a second time for in_neighbors
                        if (g.in_neighbors(a).size() >
                                g.in_neighbors(b).size()) {
                            unmatched_2 = &g.in_neighbors(a);
                            unmatched_1 = &g.in_neighbors(b);
                        } else {
                            unmatched_1 = &g.in_neighbors(a);
                            unmatched_2 = &g.in_neighbors(b);
                        }
                    }

                    // Create the cost matrix
                    num_rows = unmatched_1->size();
                    num_cols = unmatched_2->size();
                    square_size = num_cols;
                    for (size_t i = 0; i < num_rows; i++) {
                        v[i] = 0.0;
                    }
                    for (size_t i = 0; i < num_cols; i++) {
                        u[i] = 0.0;
                    }

                    r_ = 0;
                    for (auto x_itr = unmatched_1->begin();
                              x_itr != unmatched_1->end(); x_itr++) {
                        c_ = 0;
                        for (auto y_itr = unmatched_2->begin();
                                  y_itr != unmatched_2->end(); y_itr++) {

                            x = *x_itr;
                            y = *y_itr;
                            if (x > y) {
                                temp_int = x;
                                x = y;
                                y = temp_int;
                            }

                            diff = pairwise_scores[prev_mat_idx]
                                                  [start_indices[x] + (y - x)];
                            cost_matrix[r_ + square_size * c_] = diff;
                                
                            c_++;
                        }
                        r_++;
                    }

                    for (r_ = num_rows; r_ < square_size; r_++) {
                        for (c_ = 0; c_ < num_cols; c_++) {
                            cost_matrix[r_ + square_size * c_] = gap_cost;
                        }
                    }

                    diff =
                        assign2DCBasic(cost_matrix, col_for_row, row_for_col,
                                       workspace, u, v, square_size, num_cols);

                    if (diff == -1.0) {
                        throw std::logic_error(
                                "Assignment problem deemed infeasible.");
                    }
                    next_value += diff;

                }
                next_value /= scaling_factor;

                prev_value = pairwise_scores[prev_mat_idx]
                                            [start_indices[a] + (b - a)];


                pairwise_scores[next_mat_idx]
                               [start_indices[a] + (b - a)] = next_value;

                // See if the value changed by a factor of CONVERGENCE_FACTOR
                //
                // If so, we will recalculate for pairs (x, y) where x is a's
                //  neighbor and y is b's neighbor.
                if (SYM__abs_diff(next_value, prev_value) / next_value >
                        CONVERGENCE_FACTOR) {

                    for (auto x_itr = g.neighbors(a).begin();
                              x_itr != g.neighbors(a).end(); x_itr++) {
                        if ((size_t) *x_itr == node) {
                            // Distance to node is always n. Don't recalculate.
                            continue;
                        }
                        for (auto y_itr = g.neighbors(b).begin();
                                  y_itr != g.neighbors(b).end(); y_itr++) {
                            if ((size_t) *y_itr == node) {
                                // Distance to node is always n.
                                //  Don't recalculate.
                                continue;
                            }
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
                            SYM__pair_insert(next_pairs_to_calc, x, y);
                        }
                    }
                }
            }
        }

        prev_mat_idx = next_mat_idx;
        next_mat_idx = (-next_mat_idx) + 1;

        // Copy changes over in case they are not updated a second time and
        //  then revert to the old values.
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

    if (!sum_result) {
        return pairwise_scores[prev_mat_idx];  // Return the final result
    }

    // Return the (negated) total sum of pairwise_scores[prev_mat_idx]

    // This is the Kahan, Babushka, and Klein sum algorithm.
    //  Copied from the pseudocode on Wikipedia (12/8/22).
    //
    // A., Klein (2006). "A generalized Kahan-Babuska-Summation-Algorithm"
    //
    //  The point of the algorithm is to do floating-point summation while
    //   minimizing floating-point rounding errors.

    long double value;
    long double sum = 0.0;
    long double cs = 0.0;
    long double ccs = 0.0;
    long double c = 0.0;
    long double cc = 0.0;
    volatile long double t = 0.0;  // volatile ensures that the compiler doesn't
    volatile long double z = 0.0;  //   optimize these operation orders away

    for (size_t i = 0; i < matrix_size; i++) {
        // Access memoized versions if available.
        value = pairwise_scores[prev_mat_idx][i];

        t = sum + value;
        if (std::fabs(sum) >= std::fabs(value)) {
            z = sum - t;
            c = z + value;
        } else {
            z = value - t;
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
    }

    (*neg_sum) = -(sum + cs + ccs);
    return pairwise_scores[prev_mat_idx];  // Return the final result
}
