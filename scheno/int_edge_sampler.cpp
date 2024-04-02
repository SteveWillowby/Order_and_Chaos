#include "int_edge_sampler.h"
#include "thread_pool_wl_sim.h"

#include "nt_wrappers/nauty_traces.h"

#include<random>
#include<stdexcept>
#include<unordered_set>
#include<vector>

IntEdgeConverterAndSampler::IntEdgeConverterAndSampler(const Graph& g,
                                                       bool weight_samples) :
        IntEdgeConverterAndSampler(g, SparseGraph(g.directed, g.num_nodes()),
                                   weight_samples) {}

IntEdgeConverterAndSampler::IntEdgeConverterAndSampler(
              const Graph& g, const Graph& legal_edges, bool weight_samples) :
        directed(g.directed), n(g.num_nodes()), self_loops(g.num_loops() > 0) {

    edges = std::unordered_set<SCHENO__edge_int_type>();
    for (size_t a = 0; a < n; a++) {
        const auto& nbrs = g.out_neighbors(a);
        for (auto b = nbrs.begin(); b != nbrs.end(); b++) {
            if (!directed && ((size_t) *b) < a) {
                continue;
            }
            edges.insert((a * n) + *b);
        }
    }

    const std::vector<std::vector<double>> *fuzzy_orbit_sizes;
    if (weight_samples) {
        ThreadPoolWLSim tpwls(0, n);
        std::vector<std::pair<const Graph*, size_t>> tasks = {{&g, -1}};
        for (size_t a = 0; a < n; a++) {
            tasks.push_back({&g, a});
        }
        fuzzy_orbit_sizes = tpwls.get_fuzzy_orbit_sizes(&tasks);
    }

    heuristic_scores = std::vector<long double>(n*n, 0.0);
    bool self_loops = g.num_loops() > 0;
    size_t edge_int;
    long double score;
    long double total;

    // This uses the Kahan, Babushka, and Klein sum algorithm.
    //  Copied from the pseudocode on Wikipedia (12/8/22).
    //
    // A., Klein (2006). "A generalized Kahan-Babuska-Summation-Algorithm"
    //
    //  The point of the algorithm is to do floating-point summation while
    //   minimizing floating-point rounding errors.

    long double sum = 0.0;
    long double cs = 0.0;
    long double ccs = 0.0;
    long double c = 0.0;
    long double cc = 0.0;
    volatile long double t = 0.0;  // volatile ensures that the compiler doesn't
    volatile long double z = 0.0;  //   optimize these operation orders away

    for (size_t a = 0; a < n; a++) {
        for (size_t b = ((size_t) !directed) * a; b < n; b++) {
            if (!self_loops && a == b) {
                continue;
            }
            edge_int = a * n + b;
            if (legal_edges.num_edges() == 0 || legal_edges.has_edge(a, b)) {
                if (weight_samples) {
                    if (a == b) {
                        score = 2.0 * (*fuzzy_orbit_sizes)[0][a];
                    } else {
                        score = (((long double) (*fuzzy_orbit_sizes)[0][a]) * 
                                  ((long double) (*fuzzy_orbit_sizes)[a+1][b]) +
                                 ((long double) (*fuzzy_orbit_sizes)[0][b]) *
                                  ((long double) (*fuzzy_orbit_sizes)[b+1][a]));
                    }
                    // TODO: Experiment with this further.
                    score = 1.0 / std::sqrt(score);
                } else {
                    score = 1.0;
                }
            } else {
                score = 0.0;
            }

            heuristic_scores[edge_int] = score;

            t = sum + score;
            if (std::fabs(sum) >= std::fabs(score)) {
                z = sum - t;
                c = z + score;
            } else {
                z = score - t;
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
    }
    total = sum + cs + ccs;

    // Now normalize and get cumulative sums.
    cumulative_sums = std::vector<long double>(n*n, 0.0);
    sum = 0.0;
    cs = 0.0;
    ccs = 0.0;
    c = 0.0;
    cc = 0.0;
    t = 0.0;
    z = 0.0;
    for (size_t i = 0; i < n*n; i++) {
        score = heuristic_scores[i] / total;
        heuristic_scores[i] = score;

        t = sum + score;
        if (std::fabs(sum) >= std::fabs(score)) {
            z = sum - t;
            c = z + score;
        } else {
            z = score - t;
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

        cumulative_sums[i] = sum + cs + ccs;
        if (score > 0.0 && i > 0 &&
                cumulative_sums[i] == cumulative_sums[i - 1]) {
            // TODO: Verify if it's even possible that this error occur
            throw std::logic_error(
"Error! Due to lack of floating point bits, some edges cannot be sampled.");
        }
    }

    total_sum = cumulative_sums[n*n - 1];
    if (total_sum > 1.0) {
        // Hopefully this will never occur.
        //  TODO: Consider making this impossible.
        throw std::logic_error(
"Error! Due to rounding accumulation, some edges will never be sampled.");
    }
}

bool IntEdgeConverterAndSampler::is_edge(SCHENO__edge_int_type e) const {
    return edges.find(e) != edges.end();
}

Edge IntEdgeConverterAndSampler::edge(SCHENO__edge_int_type e) const {
    if (e >= n*n) {
        throw std::domain_error("Error! Edge int too large (i.e. >= n * n )");
    }
    return EDGE(e / n, e % n, directed);
}

SCHENO__edge_int_type IntEdgeConverterAndSampler::edge(const Edge& e) const {
    return (e.first * n) + e.second;
}

SCHENO__edge_int_type IntEdgeConverterAndSampler::sample(
                std::mt19937& gen,
                std::uniform_real_distribution<long double>& dist) const {

    // Sample a cumulative sum
    long double sample = dist(gen);
    while (sample >= total_sum) {
        sample = dist(gen);
    }

    // Use binary search to find the leftmost index i such that
    //  cumulative_sums[i] > sample
    size_t low = 0;
    size_t high = cumulative_sums.size();
    size_t mid = ((high - low) / 2) + low;
    while (mid > 0) {
        if (cumulative_sums[mid] <= sample) {
            low = mid;
        } else if (cumulative_sums[mid - 1] <= sample) {
            break;
        } else {
            high = mid;
        }
        mid = ((high - low) / 2) + low;
    }
    return mid;
}

SCHENO__edge_int_type IntEdgeConverterAndSampler::simple_sample(
                std::mt19937& gen,
                std::uniform_int_distribution<SCHENO__edge_int_type>& dist) const {

    SCHENO__edge_int_type x, a, b;
    while (true) {
        x = dist(gen);  // Decompose as x = a * n + b
        a = x / n;
        b = x % n;
        if (a == b && !self_loops) {
            continue;
        }
        if (!directed && a > b) {
            return (b * n) + a;
        }
        return x;
    }
}

const std::vector<long double>&
        IntEdgeConverterAndSampler::get_heuristic_scores() const {
    return heuristic_scores;
}
