#include<cmath>
#include<iostream>  // TODO: remove
#include<map>
#include<stdexcept>
#include<unordered_set>
#include<unordered_map>
#include<utility>
#include<vector>

#include "graph.h"

#include "scoring_heuristic.h"

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

long double wl_symmetry_measure(const Graph& g) {

    const double CONVERGENCE_FACTOR = 0.05;

    size_t n = g.num_nodes();
    // size_t m = g.num_edges();
    // bool has_self_loops = g.num_loops() > 0;
    bool directed = g.directed;
    // For simplicity, always leave space for self-loops
    size_t matrix_size = ((n * (n - 1)) / 2) + n;

    double score_bound = ((double) n);
    double scaling_factor = score_bound / 2.0;
    double gap_cost =       score_bound / 2.0;
    if (directed) {
        // Sums can get twice as large, so divide by twice as much
        scaling_factor = scaling_factor * 2.0;
    }

    double* pairwise_scores[3];
    pairwise_scores[0] = new double[matrix_size];
    pairwise_scores[1] = new double[matrix_size];
    pairwise_scores[2] = new double[matrix_size];  // Stores degree differences

    size_t* start_indices = new size_t[n];
    size_t start_idx = 0;
    for (size_t i = 0; i < n; i++) {
        start_indices[i] = start_idx;
        start_idx += n - i;
    }

    if (directed) {
        for (size_t a = 0; a < n; a++) {
            for (size_t b = a; b < n; b++) {
                if (a == b) {
                    pairwise_scores[0][start_indices[a]] = 0.0;
                    pairwise_scores[1][start_indices[a]] = 0.0;
                    pairwise_scores[2][start_indices[a]] = 0.0;
                    continue;
                }
                pairwise_scores[2][start_indices[a] + (b - a)] =
                    ((SYM__abs_diff((double) g.out_neighbors(a).size(),
                                    (double) g.out_neighbors(b).size()) +
                      SYM__abs_diff((double) g.in_neighbors(a).size(),
                                    (double) g.in_neighbors(b).size())) / 2.0);
                pairwise_scores[0][start_indices[a] + (b - a)] = 
                    pairwise_scores[1][start_indices[a] + (b - a)] = 
                        pairwise_scores[2][start_indices[a] + (b - a)];
            }
        }
    } else {
        for (size_t a = 0; a < n; a++) {
            for (size_t b = a; b < n; b++) {
                if (a == b) {
                    pairwise_scores[0][start_indices[a]] = 0.0;
                    pairwise_scores[1][start_indices[a]] = 0.0;
                    pairwise_scores[2][start_indices[a]] = 0.0;
                    continue;
                }
                pairwise_scores[2][start_indices[a] + (b - a)] =
                    SYM__abs_diff((double) g.neighbors(a).size(),
                                  (double) g.neighbors(b).size());

                pairwise_scores[0][start_indices[a] + (b - a)] = 
                    pairwise_scores[1][start_indices[a] + (b - a)] = 
                        pairwise_scores[2][start_indices[a] + (b - a)];
            }
        }
    }

    size_t prev_mat_idx = 2;
    size_t next_mat_idx = 1;

    double prev_value, next_value, diff;
    int a, b, x, y, temp_int;
    std::unordered_set<int> unmatched_1, unmatched_2;
    std::vector<int> shared = std::vector<int>();
    std::map<double, std::vector<std::pair<int, int>>> edge_ranks =
                        std::map<double, std::vector<std::pair<int, int>>>();

    // Keep track of which pairs of values might not have converged.
    std::unordered_map<int, std::unordered_set<int>> pairs_to_calc;
    std::unordered_map<int, std::unordered_set<int>> next_pairs_to_calc;

    for (a = 0; a < (int) n; a++) {
        for (b = a + 1; b < (int) n; b++) {
            SYM__pair_insert(pairs_to_calc, a, b);
        }
    }

    bool first_unfinished;

    std::cout<<"Started..."<<std::endl;
    while (pairs_to_calc.size() > 0) {
        std::cout<<"\t"<<pairs_to_calc.size()<<std::endl;

        next_pairs_to_calc.clear();

        first_unfinished = false;

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
                            unmatched_2 = g.out_neighbors(a);
                            unmatched_1 = g.out_neighbors(b);
                        } else {
                            unmatched_1 = g.out_neighbors(a);
                            unmatched_2 = g.out_neighbors(b);
                        }
                    } else {
                        // Go a second time for in_neighbors
                        if (g.in_neighbors(a).size() >
                                g.in_neighbors(b).size()) {
                            unmatched_2 = g.in_neighbors(a);
                            unmatched_1 = g.in_neighbors(b);
                        } else {
                            unmatched_1 = g.in_neighbors(a);
                            unmatched_2 = g.in_neighbors(b);
                        }
                    }

                    // Remove all shared nodes
                    shared.clear();
                    for (auto x_itr = unmatched_1.begin();
                              x_itr != unmatched_1.end(); x_itr++) {
                        auto y_itr = unmatched_2.find(*x_itr);
                        if (y_itr != unmatched_2.end()) {
                            shared.push_back(*x_itr);
                            unmatched_2.erase(y_itr);
                        }
                    }
                    for (auto x_itr = shared.begin();
                              x_itr != shared.end(); x_itr++) {
                        unmatched_1.erase(*x_itr);
                    }

                    // For the remaining nodes, store their pairwise costs
                    edge_ranks.clear();
                    for (auto x_itr = unmatched_1.begin();
                              x_itr != unmatched_1.end(); x_itr++) {
                        for (auto y_itr = unmatched_2.begin();
                                  y_itr != unmatched_2.end(); y_itr++) {
                            x = *x_itr;
                            y = *y_itr;
                            if (x > y) {
                                temp_int = x;
                                x = y;
                                y = temp_int;
                            }
                            // Here we just want x to be the smaller value
                            diff = pairwise_scores[prev_mat_idx]
                                                  [start_indices[x] + (y - x)];
                            auto er_itr = edge_ranks.find(diff);
                            if (er_itr == edge_ranks.end()) {
                                // Order of x vs. y matters here
                                edge_ranks.insert(std::pair<double,
                                                   std::vector<std::pair<int,int>>>(
                                                    diff, {{*x_itr, *y_itr}}));
                            } else {
                                // Order of x vs. y matters here
                                er_itr->second.push_back({*x_itr, *y_itr});
                            }

                        }
                    }

                    // Pick optimal pairs greedily, which is guaranteed to give
                    //  the optimal score due to triangle inequality held by the
                    //  overall scoring function.
                    for (auto er_itr = edge_ranks.begin();
                              er_itr != edge_ranks.end(); er_itr++) {
                        for (auto pair_itr = er_itr->second.begin();
                                  pair_itr != er_itr->second.end(); pair_itr++) {
                            x = pair_itr->first;
                            y = pair_itr->second;
                            auto x_itr = unmatched_1.find(x);
                            if (x_itr == unmatched_1.end()) {
                                continue;
                            }
                            auto y_itr = unmatched_2.find(y);
                            if (y_itr == unmatched_2.end()) {
                                continue;
                            }

                            // Both endpoints are unclaimed. Claim this edge.
                            unmatched_1.erase(x_itr);
                            unmatched_2.erase(y_itr);

                            next_value += er_itr->first;
                            if (unmatched_1.size() == 0) {
                                break;
                            }
                        }
                        if (unmatched_1.size() == 0) {
                            break;
                        }
                    }

                    if (unmatched_1.size() > 0) {
                        // TODO: Remove this check
                        throw std::logic_error("Wat?????????????");
                    }

                    // In the undirected case, this is the degree difference.
                    //
                    // In the directed case, it is the average of the two degree
                    //  differences -- since this loops twice, it works out.
                    next_value += pairwise_scores[2][start_indices[a] + (b - a)]
                                        * gap_cost;
                }
                next_value /= scaling_factor;

                prev_value = pairwise_scores[prev_mat_idx]
                                            [start_indices[a] + (b - a)];

                if (next_value < prev_value) {
                    std::cout<<"Wat????"<<std::endl;
                    std::cout<<"Old: "<<prev_value<<" vs. New: "<<next_value<<std::endl;
                    for (auto q_itr = g.out_neighbors(a).begin();
                              q_itr != g.out_neighbors(a).end(); q_itr++) {
                        for (auto z_itr = g.out_neighbors(b).begin();
                                  z_itr != g.out_neighbors(b).end(); z_itr++) {
                            std::cout<<"("<<*q_itr<<", "<<*z_itr<<"): "
                                     <<pairwise_scores[prev_mat_idx][start_indices[*q_itr] + (*z_itr - *q_itr)]<<", ";
                        }
                    }
                    std::cout<<std::endl;
                    throw std::logic_error("");
                }

                pairwise_scores[next_mat_idx]
                               [start_indices[a] + (b - a)] = next_value;

                // See if the value changed by a factor of CONVERGENCE_FACTOR
                //
                // If so, we will recalculate for pairs (x, y) where x is a's
                //  neighbor and y is b's neighbor.
                if (SYM__abs_diff(next_value, prev_value) / next_value >
                        CONVERGENCE_FACTOR) {

                    if (!first_unfinished) {
                        std::cout<<"\t\t"<<a<<", "<<b<<std::endl;
                        std::cout<<"\t\t"<<next_value<<" vs. "<<prev_value<<std::endl;
                        std::cout<<(g.neighbors(a).find(b) != g.neighbors(a).end() ?
                                    "\t\tNeighbors" : "\t\tNot Neighbors")<<std::endl;
                        std::cout<<std::endl;
                        first_unfinished = true;
                    }

                    for (auto x_itr = g.neighbors(a).begin();
                              x_itr != g.neighbors(a).end(); x_itr++) {
                        for (auto y_itr = g.neighbors(b).begin();
                                  y_itr != g.neighbors(b).end(); y_itr++) {
                            x = *x_itr;
                            y = *y_itr;
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

    delete pairwise_scores[0];
    delete pairwise_scores[1];
    delete pairwise_scores[2];
    delete start_indices;

    std::cout<<"...Done."<<std::endl;

    return -(sum + cs + ccs);
}
