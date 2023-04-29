#include "edge.h"
#include "genetic_alg_search.h"
#include "graph.h"
#include "thread_pool_scorer.h"

#include<rand>
#include<stdexcept>
#include<unordered_set>
#include<vector>

// Given a graph `g`, does a search to find good candidates for noise.
//
// Lists the top `k` candidates with their associated scores.
//
// `nt` is the number of threads to be used by the code.
std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                                 genetic_alg_search(const Graph& g,
                                                    size_t num_iterations,
                                                    size_t k,
                                                    size_t nt) {
    // Initialize Basics

    NTSparseGraph g_nt(g);
    bool directed = g_nt.directed();
    size_t num_nodes = g_nt.num_nodes();
    size_t num_edges = g_nt.num_edges();

    // Genetic Alg Variables

    const size_t DEPTH = 2;  // How many layers of gene heirarchy are there?
    size_t POP_SIZE = num_nodes * 10;
    // The top TOP_POP_SIZE performers are guaranteed to be kept.
    //  The remaining (POP_SIZE - TOP_POP_SIZE) are selected randomly from
    //      the broad population.
    size_t TOP_POP_SIZE = POP_SIZE / 2;

    size_t MATE_SIZE_A = (POP_SIZE * POP_SIZE) / 1000;
    size_t MATE_SIZE_B = POP_SIZE * 10;
    size_t MATE_SIZE = (MATE_SIZE_A < MATE_SIZE_B ? MATE_SIZE_B : MATE_SIZE_A);
    size_t MUTATE_SIZE = POP_SIZE / 10;

    // Initialization of Scoring Thread Pool
    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed));
    size_t max_flip_or_edge = num_edges * 5;

    CombinatoricUtility comb_util(max_possible_edges, max_flip_or_edge);

    NautyTracesOptions nt_options;
    nt_options.get_node_orbits = true;
    nt_options.get_edge_orbits = true;
    nt_options.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g_nt, nt_options);

    // Calculate the noise probability at which the full graph is equally
    //  likely to be the noise as it is to be the structure.
    long double alpha_plus =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                                (long double)(num_edges));
    long double log2_p_plus = std::log2l(alpha_plus) -
                              std::log2l(1.0 + alpha_plus);
    long double log2_1_minus_p_plus = -std::log2l(1.0 + alpha_plus);

    long double alpha_minus =
        std::exp2l((2.0 * ((std::log2l(nt_result.num_aut_base) +
                            (std::log2l(10) * nt_result.num_aut_exponent)) -
                             comb_util.log2_factorial(num_nodes))) /
                        (long double)(max_possible_edges - num_edges));
    long double log2_p_minus = std::log2l(alpha_minus) -
                               std::log2l(1.0 + alpha_minus);
    long double log2_1_minus_p_minus = -std::log2l(1.0 + alpha_minus);

    ThreadPoolScorer TPS(nt, g_nt, comb_util,
                         nt_result.node_orbits, nt_result.edge_orbits,
                         log2_p_plus, log2_p_minus,
                         log2_1_minus_p_plus, log2_1_minus_p_minus);


    // Initialization of Gene Population

    auto top_results =
       std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>();

    // TODO: At lower levels, make edges represented as ints, and represent
    //  genes as ints as well. That way we only need to convert to actual edges
    //  at the top level.

    // TODO: Keep track of when a gene or meta-gene first appears so that we can
    //  give them unique numeric IDs. That way, we can convert the population to
    //  an actual edge set only at the time of evaluation. (space savings)

    // TODO: Make a recursive merge function that handles conflicts.
}

Gene::Gene(size_t depth) : d(depth) {
    w = 0;
    hash = 0;
}

Gene::Gene(SYM__edge_int_type elt) : d(0) {
    w = 1;
    hash = std::hash<SYM__edge_int_type>(elt);
    e = std::vector<SYM__edge_int_type>(1, elt);
}

// Assumes elts is in sorted order
Gene::Gene(const std::vector<SYM__edge_int_type>& elts) : d(0) {
    w = elts.size();
    hash = 0;
    for (size_t i = 0; i < elts.size(); i++) {
        hash ^= std::hash<SYM__edge_int_type>(elts[i]);
    }
    e = elts;
}

Gene::Gene(Gene* elt) : d(elt->depth() + 1) {
    w = elt->weight() + 1;
    hash = std::hash<Gene*>(elt);
    sub_g = std::vector<Gene*>(1, elt);
}

// Assumes elts is in sorted order
Gene::Gene(const std::vector<Gene*>& elts) : d(elts[0]->depth() + 1) {
    w = elts.size();
    hash = 0;
    for (size_t i = 0; i < elts.size(); i++) {
        w += elts[i]->weight();
        hash ^= std::hash<Gene*>(elts[i]);
    }
    sub_g = elts;
}

// 0 means it contains edge ints
//  k > 0 means it's k steps above edge ints
size_t Gene::depth() const {
    return d;
}

// Number of relevant elements.
size_t Gene::size() const {
    if (d == 0) {
        return e.size();
    }
    return sub_g.size();
}

// Equals size() + sum over children of child.weight()
//  Used to help figure out where to mutate
size_t Gene::weight() const {
    return w;
}

const Gene::std::vector<SYM__edge_int_type>& edge_ints() const {
    if (d == 0) {
        return e;
    }

    auto sub_vecs = std::vector<std::vector<SYM__edge_int_type>>(sub_g.size());
    for (size_t i = 0; i < sub_g.size(); i++) {
        sub_vecs[i] = sub_g[i]->edge_ints();
    }

    return merged(sub_vecs);
}

// Throws an error if d == 0
const Gene::std::vector<Gene*>& sub_genes() const {
    if (d == 0) {
        throw std::domain_error("Error! No sub-genes on a gene of depth 0.");
    }
    return sub_g;
}

size_t Gene::hash_value() const {
    return hash;
}

// Note: modifies lists
std::vector<SYM__edge_int_type> Gene::merged(
        std::vector<std::vector<SYM__edge_int_type>>& lists) {
    // Aggregates the lists in a binary tree fashion
    size_t step_size = 1;
    while (step_size < lists.size()) {
        for (size_t i = 0; i + step_size < lists.size(); i += (2 * step_size)) {
            lists[i] = merged(lists[i], lists[i + step_size]);
        }
        step_size << 1;  // *= 2
    }
    return lists[0];
}

std::vector<SYM__edge_int_type> Gene::merged(
        const std::vector<SYM__edge_int_type>& a,
        const std::vector<SYM__edge_int_type>& b) {
    auto result = std::vector<SYM__edge_int_type>();
    size_t i = 0;
    size_t j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) {
            result.push_back(a[i]);
            i++;
            j++;
        } else if (a[i] < b[j]) {
            result.push_back(a[i]);
            i++;
        } else {
            result.push_back(b[j]);
            j++;
        }
    }
    if (i < a.size()) {
        while (i < a.size) {
            result.push_back(a[i]);
            i++;
        }
    } else if (j < b.size()) {
        while (j < b.size()) {
            result.push_back(b[j]);
            j++;
        }
    }
    return result;
}

// Keeps a list of the top k elements ever found
//
// Creates an initial population by starting with a single gene and
//  continuously mutating it until we get to pop size.
GenePool::GenePool(size_t gene_depth, size_t pop_size, size_t k, size_t n) :
            depth(gene_depth), pop_size(pop_size), k(k), n(n) {

    if (depth < 2) {
        throw std::domain_error("Error! Cannot make gene pool with depth < 2");
    }

    // For each depth level (except 0, which will be empty),
    //  stores a list of each Gene
    pool_vec = std::vector<std::vector<std::unique_ptr<Gene>>>();
    // Stores how many times a gene is referenced
    //  Used to know when to delete lower-level genes
    pool_counts = std::vector<std::vector<size_t>>();
    // Scores for top-level genes
    scores = std::vector<long double>();
    // For each depth level, maps a hash of a Gene to its index in pool_vec
    pool_map = std::vector<std::unordered_map<size_t, size_t>>();
    // For each depth level, the number of threads reading the pool map
    num_reading = std::vector<size_t>(depth, 0);
    // A mutex for each depth level
    pool_locks = std::vector<std::mutex>(depth);

    for (size_t i = 0; i < depth; i++) {
        if (i < depth - 1) {
            pool_counts.push_back(std::vector<size_t>());
        }
        pool_map.push_back(std::unordered_map<size_t, size_t>());
        pool_vec.push_back(std::vector<std::unique_ptr<Gene>>());
    }

    top_k = std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                                  long double>>();

    // awaiting_scores = std::vector<size_t>(); TODO: remove
}

// Grows the population by 10x
//  (creates 5x matings and 4x mutations)
// Then scores the new entries
// Lastly culls the pop back down to pop_size
void GenePool::evolve(const IntEdgeConverterAndSampler& iecas,
                      ThreadPoolScorer& tps) {

}

const std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                            long double>>& GenePool::top_k_results() const {
    return top_k;
}

// Mating "constructor"
//  Returns pointer to a Gene allocated on the heap
//  Keeps every sub-gene that's in both a and b
//      Every sub-gene that's just in a or just in b is kept with prob 1/2
//
// Requires random generators (and dists) to be passed into it for
//  thread-safety purposes
// dist should a std::uniform_int_distribution<uint_8t>(0, 1)
//
// Requires that a or b is non-empty
//
// Returns a pair (n, g). n is true iff g is a new gene.
//  If the operation fails entirely (new gene made but had hash collision)
//  then g will be NULL
std::pair<bool, Gene*> GenePool::mated(const Gene& a, const Gene& b,
                                       std::mt19937& gen,
                            std::uniform_int_distribution<uint8_t>& dist) {
    // Assumes a.depth == b.depth and a.depth > 0
    std::vector<Gene*> sub_genes = std::vector<Gene*>();
    const std::vector<Gene*>& a_sub_g = a.sub_genes();
    const std::vector<Gene*>& b_sub_g = b.sub_genes();

    size_t i, j;
    while (sub_genes.size() == 0) {
        i = 0;
        j = 0;
        while (i < a_sub_g.size() && j < b_sub_g.size()) {
            if (a_sub_g[i] == b_sub_g[j]) {
                sub_genes.push_back(a_sub_g[i]);
                i++;
                j++;
            } else if (a_sub_g[i] < b_sub_g[i]) {
                if (dist(gen)) {
                    sub_genes.push_back(a_sub_g[i]);
                }
                i++;
            } else {
                if (dist(gen)) {
                    sub_genes.push_back(b_sub_g[j]);
                }
                j++;
            }
        }
        if (i < a_sub_g.size()) {
            while (i < a_sub_g.size()) {
                if (dist(gen)) {
                    sub_genes.push_back(a_sub_g[i]);
                }
                i++;
            }
        } else {
            while (j < b_sub_g.size()) {
                if (dist(gen)) {
                    sub_genes.push_back(b_sub_g[j]);
                }
                j++;
            }
        }
    }
    // We now have a new, non-empty vector of sub-genes.

    Gene* made = new Gene(sub_genes);
    return add(made);
}

// Assumes that when add() is called, nothing is being deleted from the pool.
//  This is an important assumption for locking choices.
std::pair<bool, Gene*> GenePool::add(Gene* gene) {
    size_t hash = gene->hash_value();
    size_t d = gene->depth();

    std::unique_lock<std::mutex> l(pool_locks[depth - 1]);
    auto exists = pool_map[d].find(hash);

    size_t new_idx;

    if (exists == pool_map[d].end()) {
        // New
        new_idx = pool_vec[d].size();
        pool_vec[d].push_back(std::unique_ptr<Gene>(gene));
        score_vec[d].push_back(0);
        // if (d == depth - 1) { TODO: remove
        //     awaiting_scores.push_back(new_idx);
        // }
        pool_map[d].insert(std::pair<size_t, size_t>(hash, new_idx));
        l.unlock();

        return std::pair<bool, Gene*>(true, gene);
    }

    Gene* hash_match = pool_vec[d][*exists].get();
    l.unlock();

    if (d == 0) {
        const std::vector<SYM__edge_int_type>& edge_ints = gene->edge_ints();
        const std::vector<SYM__edge_int_type>& hm_ei = hash_match->edge_ints();
        if (hm_ei.size() != edge_ints.size()) {
            // New but hash matches
            std::cout<<"Hash collision at depth "<<d<<" with hash value "<<hash
                     <<std::endl;

            delete gene;
            return std::pair<bool, Gene*>(false, NULL);
        }
        for (size_t i = 0; i < edge_ints.size(); i++) {
            if (hm_ei[i] != edge_ints[i]) {
                // New but hash matches and of same size
                std::cout<<"Hash collision at depth "<<d<<" with hash value "
                         <<hash<<std::endl;
                delete gene;
                return std::pair<bool, Gene*>(false, NULL);
            }
        }
        // Already exists
        delete gene;
        return std::pair<bool, Gene*>(false, hash_match);
    }

    const std::vector<Gene*>& sub_genes = gene->sub_genes();
    const std::vector<Gene*>& hm_sg = hash_match->sub_genes();

    if (hm_sg.size() != sub_genes.size()) {
        // New but hash matches
        std::cout<<"Hash collision at depth "<<d<<" with hash value "<<hash
                 <<std::endl;

        delete gene;
        return std::pair<bool, Gene*>(false, NULL);
    }
    for (size_t i = 0; i < sub_genes.size(); i++) {
        if (hm_sg[i] != sub_genes[i]) {
            // New but hash matches and of same size
            std::cout<<"Hash collision at depth "<<d<<" with hash value "<<hash
                     <<std::endl;
            delete gene;
            return std::pair<bool, Gene*>(false, NULL);
        }
    }
    // Already exists
    delete gene;
    return std::pair<bool, Gene*>(false, hash_match);
}

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
std::pair<bool, Gene*> GenePool::mutated(const Gene& g,
                    const IntEdgeConverterAndSampler& iecas,
                    std::mt19937& gen,
                    std::uniform_real_distribution<SYM__edge_int_type>& disti,
                    std::uniform_real_distribution<double>& distl) {

    if (g.depth() == 0) {
        const std::vector<SYM__edge_int_type>& prev = g.edge_ints();
        std::vector<SYM__edge_int_type> next;
        if (g.size() > 1 && distl(gen) < 0.5) {
            // Remove edge int
            next.reserve(prev.size() - 1); // Prevent re-allocations for speed
            size_t to_remove = distl(gen) * ((double) prev.size());
            for (size_t i = 0; i < to_remove; i++) {
                next.push_back(prev[i]);
            }
            for (size_t i = to_remove + 1; i < prev.size(); i++) {
                next.push_back(prev[i]);
            }
        } else {
            // Add edge int
            next = std::vector<SYM__edge_int_type>(prev.size() + 1, 0);
            SYM__edge_int_type addition;
            bool done;  // Keep generating until it's new.
            size_t spot;
            while (true) {
                done = true;
                addition = iecas.sample(gen, disti);
                for (spot = 0; spot < prev.size(); spot++) {
                    if (addition > prev[spot]) {
                        next[spot] = prev[spot];
                    } else if (addition == prev[spot]) {
                        done = false;
                        break;
                    } else {
                        break;
                    }
                }
                if (!done) {
                    continue;
                }
                next[spot] = addition;
                for (; spot < prev.size(); spot++) {
                    next[spot + 1] = prev[spot];
                }
                break;
            }
        }
        Gene* made = new Gene(next);
        return add(made);
    }

    // Make 1 / (depth + 1)
    double p_mutate_here = 1.0 / ((double) g.depth() + 1.0);

    std::unique_lock<std::mutex> l(pool_locks[g.depth() - 1]);
    size_t sub_options = pool_vec[g.depth() - 1].size();
    l.unlock();

    if (sub_options > 1 && distl(gen) < p_mutate_here) {
        // Same idea as above, but with Gene* instead.

        const std::vector<Gene*>& prev = g.sub_genes();
        std::vector<Gene*> next;
        if ((g.size() > 1 && distl(gen) < 0.5) || (g.size() == sub_options)) {
            // Remove a sub-gene
            next.reserve(prev.size() - 1); // Prevent re-allocations for speed
            size_t to_remove = distl(gen) * ((double) prev.size());
            for (size_t i = 0; i < to_remove; i++) {
                next.push_back(prev[i]);
            }
            for (size_t i = to_remove + 1; i < prev.size(); i++) {
                next.push_back(prev[i]);
            }
        } else {
            // Add a sub-gene
            next = std::vector<Gene*>(prev.size() + 1, 0);
            Gene* addition;
            bool done;  // Keep generating until it's new.
            size_t spot, selection;
            while (true) {
                done = true;
                // Need to lock the pool briefly
                l.lock();
                selection = distl(gen) * pool_vec[g.depth() - 1].size();
                addition = pool_vec[g.depth() - 1][selection];
                l.unlock();
                for (spot = 0; spot < prev.size(); spot++) {
                    if (addition > prev[spot]) {
                        next[spot] = prev[spot];
                    } else if (addition == prev[spot]) {
                        done = false;
                        break;
                    } else {
                        break;
                    }
                }
                if (!done) {
                    continue;
                }
                next[spot] = addition;
                for (; spot < prev.size(); spot++) {
                    next[spot + 1] = prev[spot];
                }
                break;
            }
        }

        Gene* made = new Gene(next);
        return add(made);
    }

    // Figure out which sub-gene to mutate based on weight
    double sub_weight = g.weight() - g.size();
    double choice_weight = distl(gen) * sub_weight;
    double cumulative_weight = 0;
    const std::vector<Gene*>& sub_genes = g.sub_genes();
    size_t choice_idx;
    for (choice_idx = 0; choice_idx != sub_genes.size(); choice_idx++) {
        cumulative_weight += sub_genes[choice_idx]->weight();
        if (choice_weight < cumulative_weight) {
            break;
        }
    }
    if (choice_idx == sub_genes.size()) {  // This should never happen
        choice_idx--;
    }
    Gene* old_sub_gene = sub_genes[choice_idx];
    std::pair<bool, Gene*> x = mutated(*new_sub_gene, iecas, gen, disti, distl);
    Gene* new_sub_gene = x.second;

    if (new_sub_gene == NULL) {
        // Failed to mutate (had a hash collision)
        return x;
    }

    // The new sub-gene might be entirely new or might pre-exist.

    // TODO: Double-check this loop logic
    std::vector<Gene*> next = std::vector<Gene*>(sub_genes.size(), NULL);
    size_t old_offset = 0;
    bool added = false;
    for (size_t i = 0; i < sub_genes.size(); next++) {
        if (i == choice_idx) {
            old_offset++;
        }

        if (!added && sub_genes[i + old_offset] > new_sub_gene) {
            added = true;
            next[i] = new_sub_gene;
            old_offset--;
        } else {
            next[i] = sub_genes[i + old_offset];
        }
    }

    Gene* made = new Gene(next);
    return add(made);
}

IntEdgeConverterAndSampler::IntEdgeConverterAndSampler(const Graph& g) :
        directed(g.directed), n(g.num_nodes()), self_loops(g.num_loops() > 0),
        max_ne(((g.num_nodes() * (g.num_nodes() - 1)) / 
                (1 + size_t(!g.directed()))) +
               size_t(g.num_loops() > 0) * g.num_nodes()) {

    edges = std::unordered_set<SYM__edge_int_type>();
    for (size_t a = 0; a < n; a++) {
        const auto& nbrs = g.out_neighbors(a);
        for (auto b = nbrs.begin(); b != nbrs.end(); b++) {
            if (!directed && *b < a) {
                continue;
            }
            edges.insert((a * n) + *b);
        }
    }
}

bool IntEdgeConverterAndSampler::is_edge(SYM__edge_int_type e) const {
    return edges.find(e) != edges.end();
}

Edge IntEdgeConverterAndSampler::edge(SYM__edge_int_type e) const {
    if (e >= n*n) {
        throw std::domain_error("Error! Edge int too large (i.e. >= n * n )");
    }
    return EDGE(e / n, e % n, directed);
}

SYM__edge_int_type IntEdgeConverterAndSampler::sample(std::mt19937& gen,
                std::uniform_int_distribution<SYM__edge_int_type>& dist) const {

    SYM__edge_int_type x, a, b;
    while (true) {
        x = dist(gen);  // Decompose as x = a * n + b
        a = x / n;
        b = x % n;
        if (a == b) {
            continue;
        }
        if (!directed && a > b) {
            return b * n + a;
        }
        return x;
    }
}

GeneEdgeSet::GeneEdgeSet(const IntEdgeConverterAndSampler& iecas,
                         const Gene& g, bool mode) :
                            g(g), mode(mode), iecas(iecas) {}

const std::unordered_set<Edge, EdgeHash>& GeneEdgeSet::edges() {
    result = std::unordered_set<Edge, EdgeHash>();

    std::vector<SYM__edge_int_type> e_ints = g.edge_ints();
    size_t s = e_ints.size();
    SYM__edge_int_type e;

    if (mode == GENE_MODE_EDGES) {
        for (size_t i = 0; i < s; i++) {
            e = e_ints[i];
            if (iecas.is_edge(e)) {
                result.push_back(iecas.edge(e));
            }
        }
    } else {
        for (size_t i = 0; i < s; i++) {
            e = e_ints[i];
            if (!iecas.is_edge(e)) {
                result.push_back(iecas.edge(e));
            }
        }
    }

    return result;
}
