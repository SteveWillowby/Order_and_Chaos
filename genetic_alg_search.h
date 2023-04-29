#include "edge.h"
#include "graph.h"
#include "thread_pool_scorer.h"

#include<memory>
#include<mutex>
#include<unordered_set>
#include<utility>
#include<vector>

#ifndef SYM__GENETIC_ALG_SEARCH_H
#define SYM__GENETIC_ALG_SEARCH_H

// Given a graph `g`, does a search to find good candidates for noise.
//
// Lists the top `k` candidates with their associated scores.
//
// `nt` is the number of threads to be used by the code. Pass 0 to use
//      a default value of std::threads::hardware_concurrency();
std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                                 genetic_alg_search(const Graph& g,
                                                    size_t num_iterations,
                                                    size_t k,
                                                    size_t nt);

class GenePool;
class IntToEdgeConverter;
class GeneEdgeSet : public EdgeSet;

// Immutable object -- all operations are const
class Gene {
public:
    // Basic constructors

    Gene(size_t depth);  // The empty gene
    Gene(SYM__edge_int_type elt);
    // Assumes elts is in sorted order
    Gene(const std::vector<SYM__edge_int_type>& elts);
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

    const std::vector<SYM__edge_int_type>& edge_ints() const;
    // Throws an error if d == 0
    const std::vector<Gene*>& sub_genes() const;

    size_t hash_value() const;

protected:
    const size_t d;
    size_t w;     // Never changes - not a const to make constructor code easier
    size_t hash;  // Never changes - not a const to make constructor code easier

    // Keep in sorted order for efficiency
    std::vector<SYM__edge_int_type> e;
    // Keep in sorted order for efficiency
    //  Used only if d > 0
    std::vector<Gene*> sub_g;

    // Note: modifies lists
    static std::vector<SYM__edge_int_type> merged(
        std::vector<std::vector<SYM__edge_int_type>>& lists);

    static std::vector<SYM__edge_int_type> merged(
            const std::vector<SYM__edge_int_type>& a,
            const std::vector<SYM__edge_int_type>& b);
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
    GenePool(size_t gene_depth, size_t pop_size, size_t k, size_t n);

    // Grows the population by 10x
    //  (creates 5x matings and 4x mutations)
    // Then scores the new entries
    // Lastly culls the pop back down to pop_size
    void evolve(const IntEdgeConverterAndSampler& iecas,
                ThreadPoolScorer& tps);

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
    //  std::uniform_int_distribution<SYM__edge_int_type>(0, n * n)
    // distl should be std::uniform_real_distribution<double>(0, 1)
    //
    // Returns a pair (n, g). n is true iff g is a new gene.
    //  If the operation fails entirely (new gene made but had hash collision)
    //  then g will be NULL
    //
    // Note to self: make sure to not try to add something when nothing
    //  can be added or remove something when nothing can be removed.
    std::pair<bool, Gene*> mutated(const Gene& g,
                  const IntEdgeConverterAndSampler& iecas,
                  std::mt19937& gen,
                  std::uniform_real_distribution<SYM__edge_int_type>& disti,
                  std::uniform_real_distribution<double>& distl);

    // Adds the gene to the database
    //
    // Returns (a, b) where:
    //  a is true iff gene is added
    //  b is gene iff gene is added
    //  b is a canonical Gene iff an identical version of gene already exists
    //  b is NULL iff gene has a hash collision with a different Gene
    std::pair<bool, Gene*> add(Gene* gene);

    // Takes a gene and increments all counts of it's (sub-) occurrences.
    void take_census(Gene*);
    // Brings the top-level population down and then scans what remains with
    //  census() in order to determine which sub-genes are actually used.
    // Then goes through all sub-genes and removes all unused ones.
    void cull();

    const size_t depth;
    // Number of elements to keep at the top level.
    const size_t pop_size;
    const size_t n;
    const size_t k;

    // For each depth level (except 0, which will be empty),
    //  stores a list of each Gene
    std::vector<std::vector<std::unique_ptr<Gene>>> pool_vec;
    // Used by cull() to know when to delete lower-level genes
    std::vector<std::vector<bool>> used;
    // Scores for top-level genes
    std::vector<long double> scores;
    // For each depth level, maps a hash of a Gene to its index in pool_vec
    std::vector<std::unordered_map<size_t, size_t>> pool_map;
    // For each depth level, the number of threads reading the pool map
    std::vector<size_t> num_reading;
    // A mutex for each depth level
    std::vector<std::mutex> pool_locks;

    // The indices of all top-level genes in pool_vec awaiting scores.
    // std::vector<size_t> awaiting_scores; TODO: Remove

    std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                          long double>> top_k;
};

class IntEdgeConverterAndSampler {
public:
    IntEdgeConverterAndSampler(const Graph& g);

    // dist should be
    //  std::uniform_int_distribution<SYM__edge_int_type>(0, n * n)
    SYM__edge_int_type sample(std::mt19937& gen,
                std::uniform_int_distribution<SYM__edge_int_type>& dist) const;

    bool is_edge(SYM__edge_int_type e) const;
    Edge edge(SYM__edge_int_type e) const;

protected:
    const bool directed;
    const size_t n;
    std::unordered_set<SYM__edge_int_type> edges;
};

#define GENE_MODE_EDGES 0
#define GENE_MODE_NON_EDGES 1

// We use this in order to save RAM
class GeneEdgeSet : public EdgeSet {
public:
    GeneEdgeSet(const IntEdgeConverterAndSampler& iecas,
                const Gene& g, bool mode);

    const std::unordered_set<Edge, EdgeHash>& edges();

protected:
    // mode determines whether this returns edges or non-edges
    bool mode;
    const IntEdgeConverterAndSampler& iecas;
    const Gene& g;
};

#endif
