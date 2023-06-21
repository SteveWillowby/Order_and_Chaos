#include "basic_edge_set.h"
#include "edge.h"
#include "file_utils.h"
#include "genetic_alg_search.h"
#include "graph.h"
#include "int_edge_sampler.h"
#include "nauty_traces.h"
#include "noise_prob_choice.h"
#include "thread_pool_scorer.h"

#include<algorithm>
#include<iostream>
#include<map>
#include<memory>
#include<mutex>
#include<stdexcept>
#include<string>
#include<unordered_set>
#include<utility>
#include<vector>

std::vector<std::pair<std::unordered_set<Edge,EdgeHash>, long double>>
                 genetic_alg_search(const Graph& g,
                                    size_t num_iterations,
                                    size_t k,
                                    size_t nt,
                                    size_t gene_depth,
                                    std::unordered_set<Edge, EdgeHash> add,
                                    std::unordered_set<Edge, EdgeHash> del,
                                    const std::vector<long double>& log_probs,
                                    float max_change_factor,
                                    bool use_heuristic,
                                    const Graph& legal_edges,
                                    const std::string& file_base) {
    // Initialize Basics

    NTSparseGraph g_nt(g);
    bool directed = g.directed;
    size_t num_nodes = g_nt.num_nodes();
    size_t num_edges = g_nt.num_edges();

    // Initialization of Scoring Thread Pool
    size_t max_possible_edges =
            (num_nodes * (num_nodes - 1)) / (1 + size_t(!directed)) +
            (num_nodes * size_t(g.num_loops() > 0));

    if (num_edges == 0) {
        throw std::domain_error("Error! Cannot run on an empty graph.");
    } else if (num_edges == max_possible_edges) {
        throw std::domain_error("Error! Cannot run on a clique graph.");
    }

    size_t max_flip_or_edge = (num_edges < (max_possible_edges / 2) ?
                               num_edges : (max_possible_edges - num_edges));
    size_t max_change_size = ((double) max_flip_or_edge) * max_change_factor;
    max_flip_or_edge = ((double) (max_flip_or_edge)) *
                       (1.2 + max_change_factor);  // The 0.2 is a safety factor

    if (max_flip_or_edge < num_nodes * 2) {
        max_flip_or_edge = num_nodes * 2.0;
    }

    CombinatoricUtility comb_util(max_possible_edges, max_flip_or_edge);

    NautyTracesOptions nt_options;
    nt_options.get_node_orbits = true;
    nt_options.get_edge_orbits = true;
    nt_options.get_canonical_node_order = false;
    NautyTracesResults nt_result = traces(g_nt, nt_options);

    long double log2_p_plus = log_probs[0];
    long double log2_1_minus_p_plus = log_probs[1];
    long double log2_p_minus = log_probs[2];
    long double log2_1_minus_p_minus = log_probs[3];

    std::cout<<"log2_p_plus:          "<<log2_p_plus<<std::endl;
    std::cout<<"log2_1_minus_p_plus:  "<<log2_1_minus_p_plus<<std::endl;
    std::cout<<"log2_p_minus:         "<<log2_p_minus<<std::endl;
    std::cout<<"log2_1_minus_p_minus: "<<log2_1_minus_p_minus<<std::endl;
    std::cout<<"p_plus:  "<<std::exp2l(log2_p_plus)<<std::endl;
    std::cout<<"p_minus: "<<std::exp2l(log2_p_minus)<<std::endl;
    std::cout<<std::endl;

    std::cout<<"  Beginning edge heuristic pre-computation..."<<std::endl;
    IntEdgeConverterAndSampler iecas(g, legal_edges);
    std::cout<<"  ...Finished edge heuristic pre-computation."<<std::endl;

    ThreadPoolScorer tps(nt, g_nt, comb_util,
                         nt_result.node_orbits, nt_result.edge_orbits,
                         log2_p_plus, log2_p_minus,
                         log2_1_minus_p_plus, log2_1_minus_p_minus,
                         max_change_size, use_heuristic);

    std::vector<std::unique_ptr<EdgeSetPair>> start_task =
        std::vector<std::unique_ptr<EdgeSetPair>>();
    start_task.push_back(std::unique_ptr<EdgeSetPair>(
                            new BasicEdgeSetPair(del, add)));
    std::vector<std::pair<long double, long double>> start_score =
                    tps.get_scores(&start_task);
    std::cout<<"The special edge set gets a score of "
             <<start_score[0].first<<std::endl<<std::endl;

    // Initialization of Gene Population

    // size_t pop_size = (g.num_nodes() / 4 + 1) * 10 * // Note the 10
    //                   (g.num_nodes() < 200 ? 200 : g.num_nodes());
    size_t pop_size = 10 * size_t(std::pow(g.num_nodes() + 400, 1.4));
    GenePool gp(g, iecas, gene_depth, pop_size, k, use_heuristic);

    std::string output_graph_file = file_base + "_graph.txt";
    std::string output_noise_file = file_base + "_noise.txt";
    std::string output_nodes_file = file_base + "_nodes.txt";

    SparseGraph noise_graph(g.directed, g.num_nodes());
    SparseGraph modified_graph(g);

    for (size_t i = 0; i < num_iterations; i++) {
        std::cout<<"Beginning Iteration "<<(i + 1)<<"..."<<std::endl;

        gp.evolve(tps);

        const std::vector<std::pair<std::unordered_set<Edge,EdgeHash>,
                                    long double>>& top_k = gp.top_k_results();

        std::cout<<"...Finished iteration with a best score of "
                 <<(top_k[0].second)<<std::endl<<std::endl;

        // Write out the current best solutions
        const std::unordered_set<Edge, EdgeHash>& noise = top_k[0].first;
        for (auto edge_itr=noise.begin(); edge_itr != noise.end(); edge_itr++) {
            noise_graph.add_edge(edge_itr->first, edge_itr->second);
            modified_graph.flip_edge(edge_itr->first, edge_itr->second);
        }
        write_graph(noise_graph, output_nodes_file, output_noise_file);
        write_graph(modified_graph, output_nodes_file, output_graph_file);
        for (auto edge_itr=noise.begin(); edge_itr != noise.end(); edge_itr++) {
            noise_graph.delete_edge(edge_itr->first, edge_itr->second);
            modified_graph.flip_edge(edge_itr->first, edge_itr->second);
        }

    }

    return gp.top_k_results();
}

Gene::Gene() : d(0) {
    w = 0;
    hash = 0;
    e = std::vector<SYM__edge_int_type>();
}

Gene::Gene(SYM__edge_int_type elt) : d(0) {
    w = 1;
    hash = elt % SYM__HASH_PRIME;
    e = std::vector<SYM__edge_int_type>(1, elt);
}

// Assumes elts is in sorted order
Gene::Gene(const std::vector<SYM__edge_int_type>& elts) : d(0) {
    w = elts.size();
    hash = 0;
    for (size_t i = 0; i < elts.size(); i++) {
        // This works because the elts are in the same order each time
        //  If they could come in a different order, the hash would be different
        hash *= SYM__HASH_FACTOR;
        hash = hash % SYM__HASH_PRIME;
        hash += elts[i] % SYM__HASH_PRIME;
        hash = hash % SYM__HASH_PRIME;
    }
    e = elts;
}

Gene::Gene(Gene* elt) : d(elt->depth() + 1) {
    w = elt->weight() + 1;
    hash = ((size_t) elt) % SYM__HASH_PRIME;
    sub_g = std::vector<Gene*>(1, elt);
    e = std::vector<SYM__edge_int_type>();
}

// Assumes elts is in sorted order
Gene::Gene(const std::vector<Gene*>& elts) : d(elts[0]->depth() + 1) {
    w = elts.size();
    hash = 0;
    for (size_t i = 0; i < elts.size(); i++) {
        w += elts[i]->weight();
        // This works because the elts are in the same order each time
        //  If they could come in a different order, the hash would be different
        hash *= SYM__HASH_FACTOR;
        hash = hash % SYM__HASH_PRIME;
        hash += ((size_t) elts[i]) % SYM__HASH_PRIME;
        hash = hash % SYM__HASH_PRIME;
    }
    sub_g = elts;
    e = std::vector<SYM__edge_int_type>();
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

const std::vector<SYM__edge_int_type>& Gene::edge_ints() const {
    if (d != 0) {
        throw std::logic_error("Error! Call edge_ints() only when d = 0");
    }
    return e;
}

std::vector<SYM__edge_int_type> Gene::sub_edge_ints() {
    if (d == 0) {
        throw std::logic_error("Error! Call sub_edge_ints() only when d > 0");
    }

    if (e.size() > 0) {
        return e;
    }

    auto sub_vecs = std::vector<std::vector<SYM__edge_int_type>>(sub_g.size());
    for (size_t i = 0; i < sub_g.size(); i++) {
        if (d == 1) {
            sub_vecs[i] = sub_g[i]->edge_ints();
        } else {
            sub_vecs[i] = sub_g[i]->sub_edge_ints();
        }
    }

    e = merged(sub_vecs);

    e_hash = 0;
    for (size_t i = 0; i < e.size(); i++) {
        // This works because the elts are in the same order each time
        //  If they could come in a different order, the hash would be different
        e_hash *= SYM__HASH_FACTOR;
        e_hash = e_hash % SYM__HASH_PRIME;
        e_hash += e[i] % SYM__HASH_PRIME;
        e_hash = e_hash % SYM__HASH_PRIME;
    }
    return e;
}

// Throws an error if d == 0
const std::vector<Gene*>& Gene::sub_genes() const {
    if (d == 0) {
        throw std::domain_error("Error! No sub-genes on a gene of depth 0.");
    }
    return sub_g;
}

size_t Gene::hash_value() const {
    return hash;
}

size_t Gene::edge_int_hash_value() {
    if (d == 0) {
        return hash;
    }
    if (e.size() == 0) {
        sub_edge_ints();
    }
    return e_hash;
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
        step_size *= 2;
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
        while (i < a.size()) {
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
GenePool::GenePool(const Graph& g, const IntEdgeConverterAndSampler& iecas,
                   size_t gene_depth, size_t pop_size, size_t num_results,
                   bool use_heuristic) :
            iecas(iecas), depth(gene_depth), pop_size(pop_size), k(num_results),
            n(g.num_nodes()), use_heuristic(use_heuristic) {

    if (depth < 1) {
        throw std::domain_error("Error! Cannot make gene pool with depth < 1");
    }

    // For each depth level, stores a list of each Gene
    pool_vec = std::vector<std::vector<std::unique_ptr<Gene>>>();

    used = std::vector<std::vector<bool>>();

    // Scores for top-level genes - maps a score to a map of
    //  edge-int-hashes to gene-hashes.
    scores = std::map<std::pair<long double, long double>,
                      std::unordered_map<size_t, size_t>>();

    if (use_heuristic) {
        h_scores = std::map<std::pair<long double, long double>,
                            std::unordered_map<size_t, size_t>>();
    }

    // For each depth level, maps a hash of a Gene to its index in pool_vec
    pool_map = std::vector<std::unordered_map<size_t, size_t>>();

    // A mutex for each depth level
    pool_locks = std::vector<std::mutex>(depth);

    // Number of noise sets that should be scored in the next round
    in_need_of_scores = 1;

    // Scores for noise sets that were scored in the last generation
    carry_over_scores = std::vector<std::pair<long double, long double>>(
                            pop_size, {-INFINITY, -INFINITY});  // Dummy scores

    // Probabilities when selecting who mates/mutates
    cumulative_probs = std::vector<long double>(pop_size, 0.0);

    Gene* gene = NULL;
    for (size_t i = 0; i < depth; i++) {
        used.push_back(std::vector<bool>());
        pool_map.push_back(std::unordered_map<size_t, size_t>());
        pool_vec.push_back(std::vector<std::unique_ptr<Gene>>());

        // Build the empty gene stack
        if (i == 0) {
            gene = new Gene();
        } else {
            gene = new Gene(gene);
        }
        add(gene);
    }

    top_k = std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                                  long double>>();
}

// Score the un-scored population
//  Use scores to select who mutates and mates
//  Delete the parent generation (unless a child equals a parent)
//      (And make sure to keep the top scoring noise set)
void GenePool::evolve(ThreadPoolScorer& tps) {

    size_t prev_generation_size = pool_vec[depth - 1].size();

    std::cout<<"\tScoring"<<std::endl;

    // Get scores for new members, if any
    if (in_need_of_scores > 0) {
        size_t score_start = prev_generation_size - in_need_of_scores;

        std::vector<std::unique_ptr<EdgeSetPair>> tasks =
                    std::vector<std::unique_ptr<EdgeSetPair>>();

        for (size_t i = score_start; i < prev_generation_size; i++) {
            tasks.push_back(std::unique_ptr<EdgeSetPair>(
                    new GeneEdgeSetPair(iecas, *(pool_vec[depth - 1][i]))));
        }

        // Get scores for new members
        const std::vector<std::pair<long double, long double>>& new_scores =
                    tps.get_scores(&tasks);

        for (size_t i = 0; i < in_need_of_scores; i++) {
            carry_over_scores[score_start + i] = new_scores[i];
        }
    }

    // Now all scores are in carry_over_scores

    std::vector<bool> to_keep = std::vector<bool>(pop_size, false);
    std::vector<size_t> keeper_indices = std::vector<size_t>();
    size_t num_kept = 1;  // Always keep the top scorer
    size_t top_score_idx = 0;

    const size_t NUM_MUTATIONS = pop_size / 4;
    const size_t NUM_MATINGS = pop_size - NUM_MUTATIONS;

    // Get Selection Probabilities -- and save top k scorers

    // NOTE: Right now the code simply ignores the heuristic values
    std::set<std::pair<std::pair<long double, long double>, size_t>> top_ks;
    top_ks.insert({carry_over_scores[0], 0});

    long double max_score = carry_over_scores[0].first;
    long double min_score = max_score;
    for (size_t i = 1; i < prev_generation_size; i++) {
        if (carry_over_scores[i].first > max_score) {
            max_score = carry_over_scores[i].first;
            top_score_idx = i;
        } else if (carry_over_scores[i].first < min_score) {
            min_score = carry_over_scores[i].first;
        }

        if (top_ks.size() < k) {
            top_ks.insert({carry_over_scores[i], i});
        } else if (top_ks.begin()->first < carry_over_scores[i]) {
            top_ks.erase(top_ks.begin());
            top_ks.insert({carry_over_scores[i], i});
        }
    }
    to_keep[top_score_idx] = true;
    keeper_indices.push_back(top_score_idx);

    // Save the top k to the vector
    top_k.clear();
    for (auto tk_itr = top_ks.rbegin(); tk_itr != top_ks.rend(); tk_itr++) {
        GeneEdgeSetPair gesp(iecas, *pool_vec[depth - 1][tk_itr->second]);
        std::pair<std::unordered_set<Edge, EdgeHash>,
                  std::unordered_set<Edge, EdgeHash>> e_ne =
                                        gesp.edges_and_non_edges();
        for (auto e = e_ne.second.begin();
                  e != e_ne.second.end(); e++) {
            e_ne.first.insert(*e);
        }
        top_k.push_back(std::pair<std::unordered_set<Edge,EdgeHash>,
                                  long double>(e_ne.first,tk_itr->first.first));
    }

    long double gap = (max_score - min_score);
    if (gap == 0.0) {
        gap = 1.0;
    }

    // The best-scoring object will be
    //      (NORMALIZER + 1) / NORMALIZER times as likely to be
    //  selected. (e.g. (.25 + 1) / .25 = 5)
    //    (normalizer = (1 / k) --> (k + 1) times as likely)
    const long double NORMALIZER = 0.0625;

    for (size_t i = 0; i < prev_generation_size; i++) {
        cumulative_probs[i] = ((carry_over_scores[i].first - min_score) / gap);
        cumulative_probs[i] = cumulative_probs[i] * cumulative_probs[i]
                            + NORMALIZER;
    }
    for (size_t i = 1; i < prev_generation_size; i++) {
        cumulative_probs[i] += cumulative_probs[i - 1];
    }
    for (size_t i = 0; i < prev_generation_size - 1; i++) {
        cumulative_probs[i] /= cumulative_probs[prev_generation_size - 1];
    }
    cumulative_probs[prev_generation_size - 1] = 1.0; // Just in case x / x != 1

    // Now we have probabilities that we can sample with

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<SYM__edge_int_type> disti(0, n * n - 1);
    std::uniform_real_distribution<double> distl(0.0, 1.0);
    std::uniform_real_distribution<long double> distll(0.0, 1.0);
    std::uniform_int_distribution<uint8_t> distbool(0, 1);
    std::pair<bool, Gene*> made;

    // Perform Mutations
    std::cout<<"\tMutating"<<std::endl;

    const size_t FAILURE_MAX = 12;

    size_t selected_idx, placement_idx, second_sel_idx;

    size_t current_pop_size = num_kept;
    size_t failure_count = 0;
    while (current_pop_size < NUM_MUTATIONS) {
        // Select a gene to mutate
        selected_idx = weighted_sample(cumulative_probs, prev_generation_size,
                                        gen, distll);


        // See if we get a new gene
        made = mutated(*pool_vec[depth - 1][selected_idx],
                       gen, disti, distl, distll);

        if (made.first) {
            // New
            current_pop_size++;
            failure_count = 0;
        } else if (made.second != NULL) {
            placement_idx = pool_map[depth - 1][made.second->hash_value()];
            if (placement_idx >= prev_generation_size) {
                // Re-created one that's new from this round.
                failure_count++;
            } else if (to_keep[placement_idx]) {
                // Already re-created this once. Thus a failure.
                failure_count++;
            } else {
                to_keep[placement_idx] = true;
                keeper_indices.push_back(placement_idx);
                current_pop_size++;
                num_kept++;
                failure_count = 0;
            }
        } else {
            failure_count++;
        }

        if (failure_count == FAILURE_MAX) {
            // Too hard to mutate -- quit early.
            std::cout<<"Ran out of mutations."<<std::endl;
            break;
        }
    }

    std::cout<<"\tMating"<<std::endl;

    // Perform Matings iff we have a pop to work with.
    if (prev_generation_size > 1) {
        // Mate

        // Accounts for the possibility that some mutations failed.
        size_t target_pop_size = current_pop_size + NUM_MATINGS;

        failure_count = 0;
        while (current_pop_size < target_pop_size) {
            selected_idx = weighted_sample(cumulative_probs,
                                           prev_generation_size,
                                           gen, distll);
            second_sel_idx = selected_idx;
            while (second_sel_idx == selected_idx) {
                second_sel_idx = weighted_sample(cumulative_probs,
                                                 prev_generation_size,
                                                 gen, distll);
            }

            // See if we get a new gene
            made = mated(*pool_vec[depth - 1][selected_idx],
                         *pool_vec[depth - 1][second_sel_idx], gen, distbool);

            if (made.first) {
                // New
                current_pop_size++;
                failure_count = 0;
            } else if (made.second != NULL) {
                placement_idx = pool_map[depth - 1][made.second->hash_value()];
                if (placement_idx >= prev_generation_size) {
                    // Re-created one that's new from this round.
                    failure_count++;
                } else if (to_keep[placement_idx]) {
                    // Already re-created this once. Thus a failure.
                    failure_count++;
                } else {
                    to_keep[placement_idx] = true;
                    keeper_indices.push_back(placement_idx);
                    current_pop_size++;
                    num_kept++;
                    failure_count = 0;
                }
            } else {
                failure_count++;
            }

            if (failure_count == FAILURE_MAX) {
                // Too hard to mate -- quit early.
                std::cout<<"Ran out of matings."<<std::endl;
                break;
            }
        }
    }

    // Relevant for the next time evolve() is called.
    in_need_of_scores = current_pop_size - num_kept;
    std::cout<<"\t\t"<<(100.0 * double(in_need_of_scores) / current_pop_size)
             <<"\% of genes are new."<<std::endl;
    size_t size_total = 0;
    for (size_t x = 0; x < prev_generation_size; x++) {
        if (depth == 1) {
            size_total += pool_vec[depth - 1][x]->edge_ints().size();
        } else {
            size_total += pool_vec[depth - 1][x]->sub_edge_ints().size();
        }
    }
    std::cout<<"\t\tAverage gene size: "
             <<(double(size_total) / prev_generation_size)<<std::endl;

    // Delete old generation
    size_t i = 0;
    size_t j;
    size_t i_hash, j_hash;


    if (keeper_indices.size() > 0) {
        std::sort(keeper_indices.begin(), keeper_indices.end());
    }

    // Put all the old keepers toward the front
    for (auto idx_itr = keeper_indices.begin();
              idx_itr != keeper_indices.end(); idx_itr++) {
        j = *idx_itr;
        if (i == j) {
            i++;
            continue;
        }

        carry_over_scores[i] = carry_over_scores[j];
        i_hash = pool_vec[depth - 1][i]->hash_value();
        j_hash = pool_vec[depth - 1][j]->hash_value();
        std::swap(pool_vec[depth - 1][i], pool_vec[depth - 1][j]);
        pool_map[depth - 1].erase(i_hash);
        pool_map[depth - 1].erase(j_hash);
        pool_map[depth - 1].insert(std::pair<size_t,size_t>(j_hash, i));

        i++;
    }

    // Now move the new elements into the remaining slots
    while (pool_vec[depth - 1].size() > current_pop_size) {

        // Move the last element to overwrite idx i
        j = pool_vec[depth - 1].size() - 1;
        
        i_hash = pool_vec[depth - 1][i]->hash_value();
        j_hash = pool_vec[depth - 1][j]->hash_value();
        std::swap(pool_vec[depth - 1][i], pool_vec[depth - 1][j]);
        pool_vec[depth - 1].pop_back();
        pool_map[depth - 1].erase(i_hash);
        pool_map[depth - 1].erase(j_hash);
        pool_map[depth - 1].insert(std::pair<size_t,size_t>(j_hash, i));

        i++;
    }

    if (depth > 1) {
        // Take census and cull lower-level genes if necessary.
        cull();
    }

    //////////////////// OLD CODE /////////////////////

    /*

    // Figure out how many mutations to perform.
    const size_t mutate_factor = 6;
    size_t mutate_start_size = pool_vec[depth - 1].size();
    size_t mutate_end_size = pop_size * (mutate_factor + 1);
    if (mutate_start_size == 1) {
        if (n > 16) {
            mutate_end_size = (n * n) / 16;
        } else {
            mutate_end_size = (n * n) / 2;
        }
    } else if (mutate_start_size < pop_size) {
        mutate_end_size = mutate_start_size * mutate_start_size;
        if (mutate_end_size > pop_size * (mutate_factor + 1)) {
            mutate_end_size = pop_size * (mutate_factor + 1);
        }
    }

    // Figure out how many matings to perform.
    const size_t mate_factor = 3;
    size_t mate_start_size = mutate_end_size;
    size_t mate_end_size = mate_start_size * (mate_factor + 1);
    if (mate_start_size > pop_size) {
        mate_start_size = pop_size;
        mate_end_size = pop_size * (mate_factor + mutate_factor + 1);
    }

    // Useful variables
    size_t i, j;
    size_t quit_counter = 0;
    size_t max_quit = 0;
    // Quit early after quit_factor consecutive failed attempts.
    const size_t quit_factor = 50;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<SYM__edge_int_type> disti(0, n * n - 1);
    std::uniform_real_distribution<double> distl(0.0, 1.0);
    std::uniform_real_distribution<long double> distll(0.0, 1.0);
    std::uniform_int_distribution<uint8_t> distbool(0, 1);
    std::pair<bool, Gene*> made;

    std::cout<<"\tMutating"<<std::endl;

    // Perform mutations
    while (pool_vec[depth - 1].size() < mutate_end_size) {
        i = distl(gen) * mutate_start_size;
        made = mutated(*pool_vec[depth - 1][i], gen, disti, distl, distll);
        if (!made.first || made.second == NULL) {
            quit_counter++;
            if (quit_counter > max_quit) {
                max_quit++;
                // std::cout<<"Max quit now "<<max_quit<<" for mutation."
                //          <<std::endl;
            }
            if (quit_counter == quit_factor) {
                break;
            }
        } else if (quit_counter > 0) {
            quit_counter = 0;
        }
    }

    // If we quit early, update mate_start_size
    if (pool_vec[depth - 1].size() < mutate_end_size) {
        mate_start_size = pool_vec[depth - 1].size();
    }

    std::cout<<"\tMating"<<std::endl;

    // Perform matings
    quit_counter = 0;
    max_quit = 0;
    while (pool_vec[depth - 1].size() < mate_end_size) {
        i = distl(gen) * mate_start_size;
        j = distl(gen) * mate_start_size;
        if (i == j) {
            continue;
        }
        made = mated(*pool_vec[depth - 1][i],
                     *pool_vec[depth - 1][j], gen, distbool);

        if (!made.first || made.second == NULL) {
            quit_counter++;
            if (quit_counter > max_quit) {
                max_quit++;
                // std::cout<<"Max quit now "<<max_quit<<" for mating."
                //          <<std::endl;
            }
            if (quit_counter == quit_factor) {
                break;
            }
        } else if (quit_counter > 0) {
            quit_counter = 0;
        }
    }

    std::cout<<"\tScoring"<<std::endl;

    // Score new genes.
    size_t num_already_scored = 0;
    for (auto x = scores.begin(); x != scores.end(); x++) {
        num_already_scored += x->second.size();
    }
    size_t num_regular_scored = num_already_scored;
    if (use_heuristic) {
        for (auto x = h_scores.begin(); x != h_scores.end(); x++) {
            num_already_scored += x->second.size();
        }
    }
    size_t num_heuristic_scored = num_already_scored - num_regular_scored;

    std::vector<std::unique_ptr<EdgeSetPair>> tasks =
                std::vector<std::unique_ptr<EdgeSetPair>>();

    for (i = num_already_scored; i < pool_vec[depth - 1].size(); i++) {
        tasks.push_back(std::unique_ptr<EdgeSetPair>(
                        new GeneEdgeSetPair(iecas, *(pool_vec[depth - 1][i]))));
    }

    if (tasks.size() > 0) {
        // Get scores for new members
        const std::vector<std::pair<long double, long double>>& new_scores =
                    tps.get_scores(&tasks);

        // std::cout<<"\tPost-processing"<<std::endl;

        std::vector<size_t> keep_sizes = {pop_size};
        std::vector<std::map<std::pair<long double, long double>,
                             std::unordered_map<size_t, size_t>>*> score_maps;
        score_maps = {&scores};

        if (use_heuristic) {
            size_t to_keep = (pop_size * 3) / 4;
            keep_sizes = {to_keep + num_heuristic_scored, pop_size};
            score_maps = {&scores, &h_scores};
        }

        // Use the score map like a min-heap to keep the top population members
        std::pair<long double, long double> score;
        std::pair<long double, long double> min_score;
        size_t hash_value, e_hash_value;
        size_t orig_num_already_scored = num_already_scored;
        for (size_t smi = 0; smi < 1 + size_t(use_heuristic); smi++) {

            auto score_map = score_maps[smi];
            size_t keep_size = keep_sizes[smi];

            for (i = 0; i < tasks.size(); i++) {
                if (smi == 0) {
                    // Regular scores
                    score = new_scores[i];
                } else {
                    // Heuristic first, regular second
                    score = std::pair<long double, long double>(
                              new_scores[i].second, new_scores[i].first);
                }
                min_score = score_map->begin()->first;
                if (num_already_scored == keep_size && score < min_score) {
                    continue;  // New element not good enough to keep.
                }
                if (smi == 1 && scores.find(new_scores[i]) != scores.end()) {
                    // This score pair is already covered by scores.
                    //  Don't record it in h_scores in case we get a duplicate.

                    // TODO: Consider making this check more precise, to check
                    //  if this EXACT edge set is already taken.
                    //  However, the current setup might perform better...
                    continue;
                }

                hash_value =
                       pool_vec[depth - 1]
                               [i + orig_num_already_scored]->hash_value();
                e_hash_value =
                   pool_vec[depth - 1]
                           [i + orig_num_already_scored]->edge_int_hash_value();

                auto x = score_map->find(score);
                if (x == score_map->end()) {
                    // New score
                    score_map->insert(std::pair<std::pair<long double,
                                                          long double>,
                                        std::unordered_map<size_t, size_t>>(
                               score, {{e_hash_value, hash_value}}));
                } else {
                    // Old score -- edge set might be redundant
                    auto y = x->second.find(e_hash_value);
                    if (y == x->second.end()) {
                        // Not redundant
                        x->second[e_hash_value] = hash_value;
                    } else {
                        // Redundant -- replace with 10% prob
                        if (distl(gen) < 0.1) {
                            x->second[e_hash_value] = hash_value;
                        }
                        num_already_scored--;
                    }
                }
                num_already_scored++;

                if (num_already_scored > keep_size) {
                    num_already_scored--;
                    // Replace a smallest element.
                    auto smallest = score_map->begin();
                    if (smallest->second.size() == 1) {
                        score_map->erase(smallest);
                    } else {
                        // Remove one element from the vector.
                        // TODO: Consider making this random.
                        smallest->second.erase(smallest->second.begin());
                    }
                }
            }
        }

        // Now that we have hash values for the top `pop_size` members,
        //  keep only those top-level genes.
        top_k.clear();
        j = 0;
        size_t j_hash, i_hash;
        bool show_score;
        for (size_t smi = 0; smi < 1 + size_t(use_heuristic); smi++) {
            auto score_map = score_maps[smi];
            show_score = true;
            for (auto x = score_map->rbegin(); x != score_map->rend(); x++) {
                const std::unordered_map<size_t, size_t>& hash_values =
                                                                    x->second;
                if (show_score && use_heuristic) {
                    show_score = false;
                    std::cout<<"\tTop Score & Heuristic Values: ";
                    if (smi == 0) {
                        std::cout<<x->first.first<<", "
                                 <<x->first.second<<std::endl;
                    } else {
                        std::cout<<x->first.second<<", "
                                 <<x->first.first<<std::endl;
                    }
                }
                for (auto h = hash_values.begin(); h != hash_values.end(); h++){
                    // Swap elements i and j
                    i_hash = h->second;
                    i = pool_map[depth - 1].find(i_hash)->second;
                    if (i != j) {
                        j_hash = pool_vec[depth - 1][j]->hash_value();
                        std::swap(pool_vec[depth-1][i], pool_vec[depth-1][j]);
                        pool_map[depth - 1].erase(i_hash);
                        pool_map[depth - 1].erase(j_hash);
                        pool_map[depth - 1].insert(
                                    std::pair<size_t,size_t>(i_hash, j));
                        pool_map[depth - 1].insert(
                                    std::pair<size_t,size_t>(j_hash, i));
                    }
                    if (j < k) {
                        // Save one of the top k.
                        GeneEdgeSetPair gesp(iecas, *pool_vec[depth - 1][j]);
                        std::pair<std::unordered_set<Edge, EdgeHash>,
                                  std::unordered_set<Edge, EdgeHash>> e_ne =
                                        gesp.edges_and_non_edges();
                        for (auto e = e_ne.second.begin();
                                      e != e_ne.second.end(); e++) {
                            e_ne.first.insert(*e);
                        }
                        top_k.push_back(
                            std::pair<std::unordered_set<Edge,EdgeHash>,
                                      long double>(e_ne.first, x->first.first));
                    }
                    j++;
                }
            }
        }
        while (pool_vec[depth - 1].size() > num_already_scored) {
            std::unique_ptr<Gene> to_del =
                std::move(pool_vec[depth - 1][pool_vec[depth - 1].size() - 1]);
            pool_vec[depth - 1].pop_back();
            pool_map[depth - 1].erase(to_del->hash_value());
        }

        // Take census and cull lower-level genes if necessary.
        cull();
    } else {
        std::cout<<"Unexpected! No new genes created by evolve()!"<<std::endl;
    }
    */
}

// Takes a gene and flags it's (sub-) occurrences in `used`.
void GenePool::take_census(Gene* gene) {
    size_t d = gene->depth();

    // TODO: Add locks if this is to be made multi-threaded
    size_t idx = pool_map[d][gene->hash_value()];

    if (used[d][idx]) {
        return;
    }
    used[d][idx] = true;
    if (d == 0) {
        return;
    }

    const std::vector<Gene*>& sub_genes = gene->sub_genes();
    for (auto g_sub = sub_genes.begin(); g_sub != sub_genes.end(); g_sub++) {
        take_census(*g_sub);
    }
}

// Scans the top-level population with
//  census() in order to determine which sub-genes are actually used.
// Then goes through all sub-genes and removes all unused ones.
void GenePool::cull() {
    if (depth == 1) {
        return;
    }

    // Prepare the `used` vectors
    for (size_t i = 0; i < depth; i++) {
        used[i] = std::vector<bool>(pool_vec[i].size(), false);
    }

    for (auto gene = pool_vec[depth - 1].begin();
             gene != pool_vec[depth - 1].end(); gene++) {
        take_census(gene->get());
    }

    size_t j;
    size_t last;
    for (size_t i = 0; i < depth - 1; i++) {
        j = 0;
        last = pool_vec[i].size() - 1;
        while (j < pool_vec[i].size()) {
            if (used[i][j]) {
                j++;
                continue;
            }

            while (!used[i][last]) {
                // At least one value must be used, so there is no danger of
                //  underflow.
                pool_map[i].erase(pool_vec[i][last]->hash_value());
                pool_vec[i].pop_back();
                last--;
            }
            // since !used[i][j], we know that last will be > j or < j
            if (last < j) {
                break;
            }

            used[i][j] = used[i][last];
            pool_map[i].erase(pool_vec[i][j]->hash_value());

            if (j < last) {
                pool_vec[i][j] =
                    std::move(pool_vec[i][last]);
                pool_map[i][pool_vec[i][j]->hash_value()] = j;
            }
            pool_vec[i].pop_back();
            last--;
        }
    }
}

size_t GenePool::weighted_sample(const std::vector<long double>& cum_probs,
                                 size_t high_idx, std::mt19937& gen,
                    std::uniform_real_distribution<long double>& dist) const {

    // Sample a cumulative sum
    long double sample = dist(gen);

    // Use binary search to find the leftmost index i such that
    //  cum_probs[i] > sample
    size_t low = 0;
    size_t high = high_idx;
    size_t mid = ((high - low) / 2) + low;
    while (mid > 0) {
        if (cum_probs[mid] <= sample) {
            low = mid;
        } else if (cum_probs[mid - 1] <= sample) {
            break;
        } else {
            high = mid;
        }
        mid = ((high - low) / 2) + low;
    }
    return mid;
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

    if (a.depth() == 0) {
        std::vector<SYM__edge_int_type> elts =
            std::vector<SYM__edge_int_type>();
        const std::vector<SYM__edge_int_type>& a_elts = a.edge_ints();
        const std::vector<SYM__edge_int_type>& b_elts = b.edge_ints();

        size_t i, j;
        while (elts.size() == 0) {
            i = 0;
            j = 0;
            while (i < a_elts.size() && j < b_elts.size()) {
                if (a_elts[i] == b_elts[j]) {
                    elts.push_back(a_elts[i]);
                    i++;
                    j++;
                } else if (a_elts[i] < b_elts[j]) {
                    if (dist(gen)) {
                        elts.push_back(a_elts[i]);
                    }
                    i++;
                } else {
                    if (dist(gen)) {
                        elts.push_back(b_elts[j]);
                    }
                    j++;
                }
            }
            if (i < a_elts.size()) {
                while (i < a_elts.size()) {
                    if (dist(gen)) {
                        elts.push_back(a_elts[i]);
                    }
                    i++;
                }
            } else {
                while (j < b_elts.size()) {
                    if (dist(gen)) {
                        elts.push_back(b_elts[j]);
                    }
                    j++;
                }
            }
        }
        // We now have a new, non-empty vector of elements.

        Gene* made = new Gene(elts);
        return add(made);
    }

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
            } else if (a_sub_g[i] < b_sub_g[j]) {
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

    std::unique_lock<std::mutex> l(pool_locks[d]);
    auto exists = pool_map[d].find(hash);

    size_t new_idx;

    if (exists == pool_map[d].end()) {
        // New
        new_idx = pool_vec[d].size();
        pool_vec[d].push_back(std::unique_ptr<Gene>(gene));
        pool_map[d].insert(std::pair<size_t, size_t>(hash, new_idx));
        l.unlock();

        return std::pair<bool, Gene*>(true, gene);
    }

    Gene* hash_match = pool_vec[d][exists->second].get();
    l.unlock();

    if (d == 0) {
        const std::vector<SYM__edge_int_type>& edge_ints =
                gene->edge_ints();
        const std::vector<SYM__edge_int_type>& hm_ei =
                hash_match->edge_ints();
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
                // for (size_t j = 0; j < edge_ints.size(); j++) {
                //     std::cout<<"("<<hm_ei[j]<<" vs. "<<edge_ints[j]<<"), ";
                // }
                // std::cout<<std::endl;
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
                 <<" (case A)"<<std::endl;

        delete gene;
        return std::pair<bool, Gene*>(false, NULL);
    }
    for (size_t i = 0; i < sub_genes.size(); i++) {
        if (hm_sg[i] != sub_genes[i]) {
            // New but hash matches and of same size
            std::cout<<"Hash collision at depth "<<d<<" with hash value "<<hash
                     <<" (case B)"<<std::endl;
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
//      is added, removed, or replaced.
//
// Requires random generators (and dists) to be passed into it for
//  thread-safety purposes
// disti should be
//  std::uniform_int_distribution<SYM__edge_int_type>(0, (n * n) - 1)
// distl should be std::uniform_real_distribution<double>(0, 1)
// distll should be std::uniform_real_distribution<long double>(0, 1)
//
// Returns a pair (n, g). n is true iff g is a new gene.
//  If the operation fails entirely (new gene made but had hash collision)
//  then g will be NULL
//
// Note to self: make sure to not try to add something when nothing
//  can be added or remove something when nothing can be removed.
std::pair<bool, Gene*> GenePool::mutated(const Gene& g,
                    std::mt19937& gen,
                    std::uniform_int_distribution<SYM__edge_int_type>& disti,
                    std::uniform_real_distribution<double>& distl,
                    std::uniform_real_distribution<long double>& distll) {

    double rand;

    if (g.depth() == 0) {
        const std::vector<SYM__edge_int_type>& prev = g.edge_ints();

        SYM__edge_int_type addition;
        std::vector<SYM__edge_int_type> next;
        rand = distl(gen);

        if (g.size() > 1 && rand < 0.6) {
            // Replace edge int
            size_t remove_idx = distl(gen) * ((double) prev.size());

            bool done = false;
            size_t spot;

            while (!done) {
                done = true;
                next = std::vector<SYM__edge_int_type>(prev);
                addition = iecas.weighted_sample(gen, distll);
                for (spot = 0; spot < remove_idx; spot++) {
                    if (addition > next[spot]) {
                        continue;
                    } else if (addition == next[spot]) {
                        done = false;
                        break;
                    } else {
                        break;  // We found where it goes
                    }
                }
                if (done) {
                    if (spot < remove_idx) {
                        for (size_t i = remove_idx; i > spot; i--) {
                            next[i] = next[i - 1];
                        }
                        next[spot] = addition;
                    } else if (next[spot] == addition) {
                        // The replacement is the removed object
                        done = false;
                    } else {
                        // Still searching for where addition goes.
                        //  It's at or after remove_idx, so start shifting data
                        for (; spot < prev.size() - 1; spot++) {
                            if (addition == next[spot + 1]) {
                                done = false;
                                break;
                            } else if (addition < next[spot + 1]) {
                                break;
                            }
                            next[spot] = next[spot + 1];
                        }
                        if (next[spot] == addition) {
                            done = false;
                        } else {
                            next[spot] = addition;
                        }
                    }
                }
            }
        } else if (g.size() > 1 && rand < 0.8) {
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
            bool done = false;  // Keep generating until it's new.
            size_t spot;
            while (!done) {
                done = true;
                addition = iecas.weighted_sample(gen, distll);
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
            }
            next[spot] = addition;
            for (; spot < prev.size(); spot++) {
                next[spot + 1] = prev[spot];
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
        if (g.size() > 1 && distl(gen) < 0.6 &&
                (g.size() <= (sub_options * 2) / 3)) {
            // ^^ Make sure there are enough replacements that we have a
            //  decent chance of randomly selecting one.

            // Replace a subgene
            size_t remove_idx = distl(gen) * ((double) prev.size());

            Gene* addition;
            bool done = false;
            size_t spot, selection;

            while (!done) {
                done = true;

                next = std::vector<Gene*>(prev);
                l.lock();
                selection = distl(gen) * pool_vec[g.depth() - 1].size();
                addition = pool_vec[g.depth() - 1][selection].get();
                l.unlock();
                for (spot = 0; spot < remove_idx; spot++) {
                    if (addition > next[spot]) {
                        continue;
                    } else if (addition == next[spot]) {
                        done = false;
                        break;
                    } else {
                        break;  // We found where it goes
                    }
                }
                if (done) {
                    if (spot < remove_idx) {
                        for (size_t i = remove_idx; i > spot; i--) {
                            next[i] = next[i - 1];
                        }
                        next[spot] = addition;
                    } else if (next[spot] == addition) {
                        done = false;
                    } else {
                        // Still searching for where addition goes.
                        //  It's after or at remove_idx, so start shifting data
                        for (; spot < prev.size() - 1; spot++) {
                            if (addition == next[spot + 1]) {
                                done = false;
                                break;
                            } else if (addition < next[spot + 1]) {
                                break;
                            }
                            next[spot] = next[spot + 1];
                        }
                        if (next[spot] == addition) {
                            done = false;
                        } else {
                            next[spot] = addition;
                        }
                    }
                }
            }
        } else if (g.size() > 1 && ((g.size() >= (sub_options * 2) / 3)
                                    || distl(gen) < 0.5)) {
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
            bool done = false;  // Keep generating until it's new.
            size_t spot, selection;
            while (!done) {
                done = true;
                // Need to lock the pool briefly
                l.lock();
                selection = distl(gen) * pool_vec[g.depth() - 1].size();
                addition = pool_vec[g.depth() - 1][selection].get();
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
            }
            next[spot] = addition;
            for (; spot < prev.size(); spot++) {
                next[spot + 1] = prev[spot];
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
    std::pair<bool, Gene*> x = mutated(*old_sub_gene, gen, disti, distl,distll);
    Gene* new_sub_gene = x.second;

    if (new_sub_gene == NULL) {
        // Failed to mutate (had a hash collision)
        return x;
    }

    // The new sub-gene might be entirely new or might pre-exist.

    std::vector<Gene*> next = std::vector<Gene*>(sub_genes.size(), NULL);
    size_t old_offset = 0;
    bool added = false;
    for (size_t i = 0; i < sub_genes.size(); i++) {
        if (i == choice_idx) {
            old_offset++;
        }

        if (!added && (i + old_offset == sub_genes.size() || 
                       sub_genes[i + old_offset] > new_sub_gene)) {
            added = true;
            next[i] = new_sub_gene;
            old_offset--;
        } else {
            if (sub_genes[i + old_offset] == new_sub_gene) {
                // Failed to mutate. The new gene is the same as an old one.
                return std::pair<bool, Gene*>(false, NULL);
            }
            next[i] = sub_genes[i + old_offset];
        }
    }

    Gene* made = new Gene(next);
    return add(made);
}


GeneEdgeSetPair::GeneEdgeSetPair(const IntEdgeConverterAndSampler& iecas,
                                 Gene& g) :
                                    iecas(iecas), g(g) {h_score = 0.0;}

std::pair<std::unordered_set<Edge, EdgeHash>,
          std::unordered_set<Edge, EdgeHash>>
                         GeneEdgeSetPair::edges_and_non_edges() {

    auto result = std::pair<std::unordered_set<Edge, EdgeHash>,
                            std::unordered_set<Edge, EdgeHash>>(
                            std::unordered_set<Edge, EdgeHash>(),
                            std::unordered_set<Edge, EdgeHash>());

    std::vector<SYM__edge_int_type> e_ints;
    if (g.depth() == 0) {
        e_ints = g.edge_ints();
    } else {
        e_ints = g.sub_edge_ints();
    }

    size_t s = e_ints.size();

    if (s == 0) {
        h_score = 1.0;
    } else {
        SYM__edge_int_type e;
        h_score = 0.0;
        const std::vector<long double>& h_scores = iecas.get_heuristic_scores();

        for (size_t i = 0; i < s; i++) {
            e = e_ints[i];
            if (iecas.is_edge(e)) {
                result.first.insert(iecas.edge(e));
            } else {
                result.second.insert(iecas.edge(e));
            }
            h_score += h_scores[e];
        }
        // TODO: Consider making this the exact average (i.e. remove the + 1)
        h_score /= (long double) (e_ints.size() + 1);
    }

    return result;
}


long double GeneEdgeSetPair::heuristic_score() const {
    if (h_score == 0.0) {
        throw std::logic_error(
"Error! Must call edges_and_non_edges() before calling heuristic_score()");
    }
    return h_score;
}
