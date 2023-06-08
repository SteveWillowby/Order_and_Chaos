#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"
#include "thread_pool_scorer.h"

#ifdef SYM__THREAD_POOL_SCORER_USE_HEURISTIC
#include "Jonker_Volgenant/src/assignAlgs2D.h"
#endif

#include<condition_variable>
#include<memory>
#include<mutex>
#include<thread>
#include<unordered_set>
#include<vector>

ThreadPoolScorer::ThreadPoolScorer(size_t nt, const Graph& base_graph,
                                   const CombinatoricUtility& comb_util,
                                   const Coloring<int>& node_orbit_coloring,
                                   const Coloring<Edge,EdgeHash>&
                                            edge_orbit_coloring,
                                   long double log2_p_plus,
                                   long double log2_p_minus,
                                   long double log2_1_minus_p_plus,
                                   long double log2_1_minus_p_minus,
                                   size_t max_change_size) : 
                                      num_threads(
                         (nt == 0 ? std::thread::hardware_concurrency() : nt)),
                                      comb_util(comb_util),
                                      node_orbit_coloring(node_orbit_coloring),
                                      edge_orbit_coloring(edge_orbit_coloring),
                                      log2_p_plus(log2_p_plus),
                                      log2_p_minus(log2_p_minus),
                                      log2_1_minus_p_plus(log2_1_minus_p_plus),
                                      log2_1_minus_p_minus(log2_1_minus_p_minus),
                                      max_change_size(max_change_size)
{
    terminate_pool = false;
    tasks_begun = 0;
    num_tasks = 0;
    threads_working = threads_launched = 0;
    score_vec = std::vector<std::pair<long double, long double>>();
    graphs = std::vector<NTSparseGraph>();
    edge_colorings = std::vector<Coloring<Edge, EdgeHash>>();
    pool = std::vector<std::thread>();

#ifdef SYM__THREAD_POOL_SCORER_USE_HEURISTIC
    size_t n = base_graph.num_nodes();

    start_indices = new size_t[n];
    size_t start_idx = 0;
    for (size_t i = 0; i < n; i++) {
        start_indices[i] = start_idx;
        start_idx += n - i;
    }
    difference_matrices_1 = std::vector<double*>();
    difference_matrices_2 = std::vector<double*>();
    u_vec = std::vector<double*>();
    v_vec = std::vector<double*>();
    cost_matrices = std::vector<double*>();
    col_for_row_vec = std::vector<ptrdiff_t*>();
    row_for_col_vec = std::vector<ptrdiff_t*>();
    workspaces = std::vector<void*>();
#endif

    for (size_t i = 0; i < num_threads; i++) {
        graphs.push_back(NTSparseGraph(base_graph));
        edge_colorings.push_back(Coloring<Edge, EdgeHash>());

#ifdef SYM__THREAD_POOL_SCORER_USE_HEURISTIC
        u_vec.push_back(new double[n]);
        v_vec.push_back(new double[n]);
        difference_matrices_1.push_back(new double[n * n]);
        difference_matrices_2.push_back(new double[n * n]);
        cost_matrices.push_back(new double[n * n]);
        col_for_row_vec.push_back(new ptrdiff_t[n]);
        row_for_col_vec.push_back(new ptrdiff_t[n]);
        workspaces.push_back(malloc(assign2DCBufferSize(n, n)));
#endif

        for (auto c = edge_orbit_coloring.colors().begin();
                  c != edge_orbit_coloring.colors().end(); c++) {
            const auto &cell = edge_orbit_coloring.cell(*c);
            for (auto edge = cell.begin(); edge != cell.end(); edge++) {
                edge_colorings[i].set(*edge, *c);
            }
        }
    }

    std::unique_lock<std::mutex> launch_lock(m_launch);
    for (size_t i = 0; i < num_threads; i++) {
        pool.push_back(std::thread(&ThreadPoolScorer::run, this));
        launch_wait_signal.wait(launch_lock);
    }
    launch_lock.unlock();
}

ThreadPoolScorer::~ThreadPoolScorer() {
    terminate();

#ifdef SYM__THREAD_POOL_SCORER_USE_HEURISTIC
    delete start_indices;
    for (size_t i = 0; i < num_threads; i++) {
        delete u_vec[i];
        delete v_vec[i];
        delete cost_matrices[i];
        delete difference_matrices_1[i];
        delete difference_matrices_2[i];
        delete col_for_row_vec[i];
        delete row_for_col_vec[i];
        free(workspaces[i]);
    }
#endif
}

const std::vector<std::pair<long double, long double>>&
            ThreadPoolScorer::get_scores(
                    std::vector<std::unique_ptr<EdgeSetPair>> *tasks) {

    std::unique_lock<std::mutex> begin_lock(m_worker_meta);
    std::unique_lock<std::mutex> done_lock(m_tasks_done);

    // Prepare jobs.
    task_vec = tasks;
    num_tasks = task_vec->size();
    if (score_vec.size() < num_tasks) {
        score_vec =
            std::vector<std::pair<long double, long double>>(num_tasks, {0,0});
    }
    tasks_begun = 0;
    threads_working = 0;

    // Awaken workers.
    begin_lock.unlock();
    worker_wait_signal.notify_all();

    // Wait for completion.
    done_wait_signal.wait(done_lock);
    done_lock.unlock();

    return score_vec; 
}

void ThreadPoolScorer::terminate() {
    if (terminate_pool) {
        return;
    }

    terminate_pool = true;
    worker_wait_signal.notify_all();

    for (size_t i = 0; i < num_threads; i++) {
        pool[i].join();
    }
}

void ThreadPoolScorer::run() {
    std::unique_lock<std::mutex> launch_lock(m_launch);
    size_t thread_id = threads_launched;
    threads_launched++;
    launch_lock.unlock();
    launch_wait_signal.notify_one();

    size_t task_id;

    std::unique_lock<std::mutex> done_lock(m_tasks_done);
    done_lock.unlock();

    std::unique_lock<std::mutex> meta_lock(m_worker_meta);
    while (!terminate_pool) {
        // Wait for the run signal to begin work.
        if (num_tasks == 0 && !terminate_pool) {
            // Releases lock. Re-acquires it when wakened.
            worker_wait_signal.wait(meta_lock);
        }
        if (terminate_pool) {
            break;
        }
        threads_working++;
        meta_lock.unlock();

        // Do scoring.
        std::unique_lock<std::mutex> l(m_worker_queue);  // Locks m_worker_queue
        while (tasks_begun < num_tasks) {
            if (tasks_begun < num_tasks) {
                task_id = tasks_begun;
                tasks_begun++;
                l.unlock();
            } else {
                break;
            }

            std::pair<std::unordered_set<Edge, EdgeHash>,
                      std::unordered_set<Edge, EdgeHash>> e_ne =
                            (*task_vec)[task_id]->edges_and_non_edges();

            // Do task # task_id
#ifdef SYM__THREAD_POOL_SCORER_USE_HEURISTIC
            score_vec[task_id] = score(graphs[thread_id], comb_util,
                                       node_orbit_coloring,
                                       edge_orbit_coloring,
                                       edge_colorings[thread_id],
                                       e_ne.second,
                                       e_ne.first,
                                       log2_p_plus, log2_p_minus,
                                       log2_1_minus_p_plus,
                                       log2_1_minus_p_minus,
                                       max_change_size,
                                       cost_matrices[thread_id],
                                       col_for_row_vec[thread_id],
                                       row_for_col_vec[thread_id],
                                       u_vec[thread_id], v_vec[thread_id],
                                       workspaces[thread_id],
                                       difference_matrices_1[thread_id],
                                       difference_matrices_2[thread_id],
                                       start_indices);
#endif
#ifndef SYM__THREAD_POOL_SCORER_USE_HEURISTIC

            // Right now, no heuristic is used at all.
            score_vec[task_id] =
                std::pair<long double, long double>(
                    score(graphs[thread_id], comb_util,
                                       node_orbit_coloring,
                                       edge_orbit_coloring,
                                       edge_colorings[thread_id],
                                       e_ne.second,
                                       e_ne.first,
                                       log2_p_plus, log2_p_minus,
                                       log2_1_minus_p_plus,
                                       log2_1_minus_p_minus,
                                       max_change_size),
                    0);

                    //  This makes the program get stuck -- too many edges on
                    //    the most unique nodes.
                    // (*task_vec)[task_id]->heuristic_score());

#endif

            l.lock();
        }
        l.unlock();

        meta_lock.lock();
        threads_working--;
        if (threads_working == 0) {
            // Let get_scores() and all other workers know we're done.
            done_lock.lock();
            num_tasks = 0;
            done_lock.unlock();
            worker_wait_signal.notify_all();
            done_wait_signal.notify_all();
        } else {
            // Pause for everyone else to finish.
            worker_wait_signal.wait(meta_lock);
        }
    }
    meta_lock.unlock();
}
