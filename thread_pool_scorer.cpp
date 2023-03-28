#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"
#include "scoring_function.h"
#include "thread_pool_scorer.h"

#include<thread>
#include<unordered_set>
#include<vector>

ThreadPoolScorer::ThreadPoolScorer(size_t num_threads, const Graph& base_graph,
                                   const CombinatoricUtility& comb_util,
                                   const Coloring<int>& node_orbit_coloring,
                                   const Coloring<Edge,EdgeHash>&
                                            edge_orbit_coloring,
                                   double log2_p_plus,
                                   double log2_p_minus,
                                   double log2_1_minus_p_plus,
                                   double log2_1_minus_p_minus) : 
                                      num_threads(
       (num_threads == 0 ? std::thread::hardware_concurrency() : num_threads)),
                                      comb_util(comb_util),
                                      node_orbit_coloring(node_orbit_coloring),
                                      edge_orbit_coloring(edge_orbit_coloring),
                                      log2_p_plus(log2_p_plus),
                                      log2_p_minus(log2_p_minus),
                                      log2_1_minus_p_plus(log2_1_minus_p_plus),
                                      log2_1_minus_p_minus(log2_1_minus_p_minus)
{
    terminate_pool = false;
    tasks_begun = 0;
    tasks_finished = 0;
    num_tasks = 0;

    for (size_t i = 0; i < num_threads; i++) {
        graphs.push_back(NTSparseGraph(base_graph));
        colorings.push_back(Coloring<Edge, EdgeHash>());
        for (auto &c = edge_orbit_coloring.begin();
                  c != edge_orbit_coloring.end(); c++) {
            const auto &cell = edge_orbit_coloring.cell(*c);
            for (auto &edge = cell.begin(); edge != cell.end(); edge++) {
                colorings[i].set(*edge, *c);
            }
        }
    }

    threads_launched = 0;
    for (size_t i = 0; i < num_threads; i++) {
        pool.push_back(std::thread(&ThreadPoolScorer::run, this));
        // the run() method increments threads_launched
        while (threads_launched == i) {}
    }
}

ThreadPoolScorer::~ThreadPoolScorer() {
    terminate();
}

const std::vector<double>& ThreadPoolScorer::get_scores(
                                std::vector<std::pair<
                                  std::unordered_set<Edge, EdgeHash>,
                                  std::unordered_set<Edge, EdgeHash>>> *tasks) {
    task_vec = tasks;
    num_tasks = task_vec->size();
    tasks_finished = tasks_begun = 0;

}

void ThreadPoolScorer::terminate() {
    terminate_pool = true;

    for (size_t i = 0; i < num_threads; i++) {
        pool[i].join();
    }
}

void ThreadPoolScorer::run() {
    size_t thread_id = threads_launched;
    threads_launched++;

    size_t task_id;

    while (!terminate_pool) {
        // Wait for the run signal to begin work.
        std::unique_lock<std::mutex> lock(m_meta);
        if (num_tasks == 0) {
            wait_signal.wait(lock, [this](){
                return num_tasks > 0 || terminate_pool;
            });
        } else {
            lock.unlock();
        }

        if (terminate_pool) {
            break;
        }

        // Do scoring.
        while (tasks_begun < num_tasks && !terminate_pool) {
            std::unique_lock<std::mutex> l(m_queue);
            if (tasks_begun < num_tasks) {
                task_id = tasks_begun;
                tasks_begun++;
                l.unlock();
            } else {
                l.unlock();
                break;
            }

            // Do task # task_id
        }
    }

}
