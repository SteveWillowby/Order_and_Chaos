#include "int_edge_sampler.h"
#include "thread_pool_scorer.h"

#include "nt_wrappers/nauty_traces.h"

#include<map>
#include<memory>
#include<mutex>
#include<random>
#include<string>
#include<unordered_set>
#include<utility>
#include<vector>

#ifndef SCHENO__GENETIC_ALG_SEARCH_H
#define SCHENO__GENETIC_ALG_SEARCH_H

// Given a graph `g`, does a search to find good candidates for noise.
//
// Lists the top `k` candidates with their associated scores.
//
// `nt` is the number of threads to be used by the code. Pass 0 to use
//      a default value of std::threads::hardware_concurrency();
//
// `num_iterations` is the number of mate, mutate, and score steps.
//
// `gene_depth` is how many layers of sub-genes there are.
//      A value of 1 means a gene is an edge set
//      A value of 2 means a gene is a set of edge sets
//      A value of 3 means a gene is a set of sets of edge sets
//          Etc.
//
//  `log_probs` should contain {log2(p+), log2(1-p+), log2(p-), log2(1-p-)}
//      where p+ is the prob. a noise edge was added to get `g`, and p- is the
//      prob. an edge was randomly removed to get `g`.
//
//  `max_change_factor` is used to limit how large the solution sets can be.
//      A set of changes with more than `max_change_factor` * g.num_edges()
//          will receive an infinitely bad score.
//
//  `scoring_heuristic` determines whether or not noise sets get a heuristic
//      score as a tiebreaker
//
//  `sampling_heuristic` determines whether or not random edge sampling is
//      weighted by a heuristic
//
//  `seed_noise` is the initial noise that the genetic algorithm considers.
//      If you want to start with no noise, just make it a graph with no edges.
//
//  `legal_edges` is used if you want to constrain which edges can be part of
//      the noise. If you want to include all edges as an option, make
//      `legal_edges` an empty graph on the same number of nodes as `g`.
//
//  `file_base` says what to write the outputs to
//
//  `full_iso` means to use an algorithm like `traces`
//      if `full_iso` is set to false, then an approximate algorithm is used
//      for automorphism calculations
std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                 genetic_alg_search(const Graph& g,
                                    size_t num_iterations,
                                    size_t k,
                                    size_t nt,
                                    size_t gene_depth,
                                    const std::vector<long double>& log_probs,
                                    float max_change_factor,
                                    bool scoring_heuristic,
                                    bool sampling_heuristic,
                                    const Graph& seed_noise,
                                    const Graph& legal_edges,
                                    const std::string& file_base,
                                    bool full_iso);

class GenePool;
class GeneEdgeSetPair;

// Immutable object -- all operations are const
class Gene {
public:
    // Basic constructors

    Gene();  // The empty, depth-0 gene
    Gene(SCHENO__edge_int_type elt);
    // Assumes elts is in sorted order
    Gene(const std::vector<SCHENO__edge_int_type>& elts);
    Gene(Gene* elt);
    // Assumes elts is in sorted order
    // Assumes each elt is of the same depth
    Gene(const std::vector<Gene*>& elts);

    // 0 means it contains edge ints
    //  k > 0 means it's k steps above edge ints
    size_t depth() const;
    // Number of relevant elements.
    size_t size() const;
    // Equals size() + sum over children of child.weight()
    //  Used to help figure out where to mutate
    size_t weight() const;

    const std::vector<SCHENO__edge_int_type>& edge_ints() const;
    std::vector<SCHENO__edge_int_type> sub_edge_ints();
    // Throws an error if d == 0
    const std::vector<Gene*>& sub_genes() const;

    size_t hash_value() const;

    size_t edge_int_hash_value();

protected:
    const size_t d;
    size_t w;     // Never changes - not a const to make constructor code easier
    size_t hash;  // Never changes - not a const to make constructor code easier
    size_t e_hash;  // Assigned a value when sub_edge_ints() is called

    // Keep in sorted order for efficiency
    std::vector<SCHENO__edge_int_type> e;
    // Keep in sorted order for efficiency
    //  Used only if d > 0
    std::vector<Gene*> sub_g;

    // Note: modifies lists
    static std::vector<SCHENO__edge_int_type> merged(
        std::vector<std::vector<SCHENO__edge_int_type>>& lists);

    static std::vector<SCHENO__edge_int_type> merged(
            const std::vector<SCHENO__edge_int_type>& a,
            const std::vector<SCHENO__edge_int_type>& b);
};

// Contains a record of each distinct gene and a hash ID for them
//
// Runs a thread-pool
class GenePool {
public:
    // Keeps a list of the top k elements ever found
    //
    // Creates an initial population by starting with a single gene and
    //  continuously mutating it until we get to pop size.
    GenePool(const Graph& g, const IntEdgeConverterAndSampler& iecas,
             const Graph& seed_noise,
             size_t gene_depth, size_t pop_size, size_t num_results,
             bool scoring_heuristic, bool sampling_heuristic);

    // Grows the population by 10x
    //  (creates 3x matings and 6x mutations)
    // Then scores the new entries
    // Lastly culls the pop back down to pop_size
    void evolve(ThreadPoolScorer& tps);

    const std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                                long double>>& top_k_results() const;

protected:

    // Mating "constructor"
    //  Returns pointer to a Gene allocated on the heap
    //  Keeps every sub-gene that's in both a and b
    //      Every sub-gene that's just in a or just in b is kept with prob 1/2
    //
    // Requires random generators (and dists) to be passed into it for
    //  thread-safety purposes
    // dist should a std::uniform_int_distribution<uint_8t>(0, 1)
    //
    // Returns a pair (n, g). n is true iff g is a new gene.
    //  If the operation fails entirely (new gene made but had hash collision)
    //  then g will be NULL
    std::pair<bool, Gene*> mated(const Gene& a, const Gene& b,
                std::mt19937& gen,
                std::uniform_int_distribution<uint8_t>& dist);

    // Mutate "constructor"
    //  Returns pointer to a Gene allocated on the heap
    //  Creates a new gene identical to g but where a gene or sub-gene
    //      is added or removed.
    //
    // Requires random generators (and dists) to be passed into it for
    //  thread-safety purposes
    // disti should be
    //  std::uniform_int_distribution<SCHENO__edge_int_type>(0, (n * n) - 1)
    // distl should be std::uniform_real_distribution<double>(0, 1)
    // distll should be std::uniform_real_distribution<long double>(0, 1)
    //
    // Returns a pair (n, g). n is true iff g is a new gene.
    //  If the operation fails entirely (new gene made but had hash collision)
    //  then g will be NULL
    //
    // Note to self: make sure to not try to add something when nothing
    //  can be added or remove something when nothing can be removed.
    std::pair<bool, Gene*> mutated(const Gene& g,
                  std::mt19937& gen,
                  std::uniform_int_distribution<SCHENO__edge_int_type>& disti,
                  std::uniform_real_distribution<double>& distl,
                  std::uniform_real_distribution<long double>& distll);

    // Adds the gene to the database
    //
    // Returns (a, b) where:
    //  a is true iff gene is added
    //  b is gene iff gene is added
    //  b is a canonical Gene iff an identical version of gene already exists
    //  b is NULL iff gene has a hash collision with a different Gene
    std::pair<bool, Gene*> add(Gene* gene);

    // Takes a gene and flags it's (sub-) occurrences in `used`.
    void take_census(Gene* gene);
    // Scans the top-level population with
    //  census() in order to determine which sub-genes are actually used.
    // Then goes through all sub-genes and removes all unused ones.
    void cull();

    const IntEdgeConverterAndSampler iecas;

    const size_t depth;
    // Number of elements to keep at the top level.
    const size_t pop_size;
    const size_t k;
    const size_t n;
    const bool scoring_heuristic;
    const bool sampling_heuristic;

    // For each depth level, stores a list of each Gene
    std::vector<std::vector<std::unique_ptr<Gene>>> pool_vec;
    // Used by cull() to know when to delete lower-level genes
    std::vector<std::vector<bool>> used;

    // Scores for top-level genes - maps a score to a map of
    //  edge-int-hashes to gene-hashes.
    std::map<std::pair<long double, long double>,
             std::unordered_map<size_t, size_t>> scores;

    // Same as above but for heuristic scores getting the priority
    std::map<std::pair<long double, long double>,
             std::unordered_map<size_t, size_t>> h_scores;

    // For each depth level, maps a hash of a Gene to its index in pool_vec
    std::vector<std::unordered_map<size_t, size_t>> pool_map;

    // For each depth level, the number of threads reading the pool map
    // std::vector<size_t> num_reading; // TODO: remove

    // A mutex for each depth level
    std::vector<std::mutex> pool_locks;

    // The indices of all top-level genes in pool_vec awaiting scores.
    // std::vector<size_t> awaiting_scores; TODO: Remove

    std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                          long double>> top_k;
};

// Make sure that SCHENO__HASH_PRIME * SCHENO__HASH_FACTOR < 0xFFFFFFFFFFFFFFFF
#define SCHENO__HASH_PRIME 87178291199
#define SCHENO__HASH_FACTOR 700001
// More primes: 233, 1597, 28657, 514229, 700001, 2971215073,
//              87178291199, 3314192745739

// We use this in order to save RAM
class GeneEdgeSetPair : public EdgeSetPair {
public:
    GeneEdgeSetPair(const IntEdgeConverterAndSampler& iecas,
                    Gene& g);

    std::pair<std::unordered_set<Edge, EdgeHash>,
              std::unordered_set<Edge, EdgeHash>>
                                 edges_and_non_edges();

    // Do not call until after edges_and_non_edges()
    long double heuristic_score() const;

protected:
    const IntEdgeConverterAndSampler& iecas;
    Gene& g;
    long double h_score;
};

#endif
