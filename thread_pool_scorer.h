#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"

#include<unordered_set>
#include<vector>

#ifndef SYM__THREAD_POOL_SCORER_H
#define SYM__THREAD_POOL_SCORER_H

class ThreadPoolScorer {
public:
    ThreadPoolScorer(size_t num_threads, const Graph& base_graph,
                     const CombinatoricUtility& comb_util,
                     const Coloring<int>& node_orbit_coloring,
                     const Coloring<Edge,EdgeHash>& edge_orbit_coloring,
                     double log2_p_plus, double log2_p_minus,
                     double log2_1_minus_p_plus, double log2_1_minus_p_minus);

    ~ThreadPoolScorer();

    const std::vector<double>& get_scores(
            std::vector<std::pair<
                     std::unordered_set<Edge, EdgeHash>,
                     std::unordered_set<Edge, EdgeHash>>> *tasks);

    void terminate();

protected:
    const size_t num_threads;

    std::vector<std::thread> pool;
    std::vector<NTSparseGraph> graphs;  // One copy of the graph per thread.
    // One editable coloring per thread.
    std::vector<Coloring<Edge, EdgeHash>> colorings;

    std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                          std::unordered_set<Edge, EdgeHash>>> *task_vec;
    std::vector<double> score_vec;

    bool terminate_pool;

    size_t tasks_begun, tasks_finished, num_tasks;
    size_t threads_launched;

    const CombinatoricUtility& comb_util;
    const Coloring<int>& node_orbit_coloring;
    const Colorint<Edge,EdgeHash>& edge_orbit_coloring;
    const double log2_p_plus;
    const double log2_p_minus;
    const double log2_1_minus_p_plus;
    const double log2_1_minus_p_minus;

    std::mutes m_queue;
    std::mutex m_meta;
    std::condition_variable wait_signal;

    void run();
};

#endif
