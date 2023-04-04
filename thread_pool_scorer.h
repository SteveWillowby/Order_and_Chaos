#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"

#include<condition_variable>
#include<mutex>
#include<thread>
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
                     long double log2_p_plus, long double log2_p_minus,
                     long double log2_1_minus_p_plus,
                     long double log2_1_minus_p_minus);

    ~ThreadPoolScorer();

    // IMPORTANT: The vector returned may have extra elements at the end.
    //  Thus you must ignore get_scores(tasks).size().
    //  The score for (*tasks)[i] will be stored at get_scores(tasks)[i].
    //
    // The pairs should be (edge_additions, edge_removals) pairs.
    const std::vector<long double>& get_scores(
            std::vector<std::pair<
                     std::unordered_set<Edge, EdgeHash>,
                     std::unordered_set<Edge, EdgeHash>>> *tasks);

    void terminate();

protected:
    const size_t num_threads;

    std::vector<std::thread> pool;
    std::vector<NTSparseGraph> graphs;  // One copy of the graph per thread.
    // One editable coloring per thread.
    std::vector<Coloring<Edge, EdgeHash>> edge_colorings;

    std::vector<std::pair<std::unordered_set<Edge, EdgeHash>,
                          std::unordered_set<Edge, EdgeHash>>> *task_vec;
    std::vector<long double> score_vec;

    bool terminate_pool;

    size_t tasks_begun, num_tasks;
    size_t threads_working, threads_launched;

    const CombinatoricUtility& comb_util;
    const Coloring<int>& node_orbit_coloring;
    const Coloring<Edge,EdgeHash>& edge_orbit_coloring;
    const long double log2_p_plus;
    const long double log2_p_minus;
    const long double log2_1_minus_p_plus;
    const long double log2_1_minus_p_minus;

    // Used by running workers to select the next task
    std::mutex m_worker_queue;
    // Used by workers (and get_scores()) to begin and end working on a set of
    //  tasks
    std::mutex m_worker_meta;
    // Used by workers (and get_scores()) to say when all current tasks are done
    std::mutex m_tasks_done;
    // Used in conjunction with m_worker_meta
    std::condition_variable worker_wait_signal;
    // Used in conjunction with m_tasks_done
    std::condition_variable done_wait_signal;

    void run();
};

#endif
