#include "coloring.h"
#include "edge.h"
#include "nt_sparse_graph.h"
#include "Jonker_Volgenant/src/assignAlgs2D.h"
#include "wl_similarity.h"
#include "thread_pool_wl_sim.h"

#include<condition_variable>
#include<memory>
#include<mutex>
#include<thread>
#include<unordered_set>
#include<vector>

ThreadPoolWLSim::ThreadPoolWLSim(size_t nt, size_t n) : num_threads(
                         (nt == 0 ? std::thread::hardware_concurrency() : nt)),
                         num_nodes(n)
{
    terminate_pool = false;
    tasks_begun = 0;
    num_tasks = 0;
    threads_working = threads_launched = 0;
    uniqueness_vec = std::vector<std::vector<double>>();
    pool = std::vector<std::thread>();

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

    for (size_t i = 0; i < num_threads; i++) {
        u_vec.push_back(new double[n]);
        v_vec.push_back(new double[n]);
        difference_matrices_1.push_back(new double[n * n]);
        difference_matrices_2.push_back(new double[n * n]);
        cost_matrices.push_back(new double[n * n]);
        col_for_row_vec.push_back(new ptrdiff_t[n]);
        row_for_col_vec.push_back(new ptrdiff_t[n]);
        workspaces.push_back(malloc(assign2DCBufferSize(n, n)));
    }

    std::unique_lock<std::mutex> launch_lock(m_launch);
    for (size_t i = 0; i < num_threads; i++) {
        pool.push_back(std::thread(&ThreadPoolWLSim::run, this));
        launch_wait_signal.wait(launch_lock);
    }
    launch_lock.unlock();
}

ThreadPoolWLSim::~ThreadPoolWLSim() {
    terminate();

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
}

const std::vector<std::vector<double>>& ThreadPoolWLSim::get_uniqueness(
                const std::vector<std::pair<const Graph*, size_t>>* tasks) {

    std::unique_lock<std::mutex> begin_lock(m_worker_meta);
    std::unique_lock<std::mutex> done_lock(m_tasks_done);

    // Prepare jobs.
    task_vec = tasks;
    num_tasks = task_vec->size();
    while (uniqueness_vec.size() < num_tasks) {
        uniqueness_vec.push_back(std::vector<double>(num_nodes, 0.0));
    }
    tasks_begun = 0;
    threads_working = 0;

    // Awaken workers.
    begin_lock.unlock();
    worker_wait_signal.notify_all();

    // Wait for completion.
    done_wait_signal.wait(done_lock);
    done_lock.unlock();

    return uniqueness_vec; 
}

void ThreadPoolWLSim::terminate() {
    if (terminate_pool) {
        return;
    }

    terminate_pool = true;
    worker_wait_signal.notify_all();

    for (size_t i = 0; i < num_threads; i++) {
        pool[i].join();
    }
}

void ThreadPoolWLSim::run() {
    std::unique_lock<std::mutex> launch_lock(m_launch);
    size_t thread_id = threads_launched;
    threads_launched++;
    launch_lock.unlock();
    launch_wait_signal.notify_one();

    size_t task_id, a, b;

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

            std::pair<const Graph*, size_t> task = (*task_vec)[task_id];

            // Do task # task_id

            double* sim =wl_similarity_measure(false, NULL,
                                               *(task.first), task.second,
                                               cost_matrices[thread_id],
                                               col_for_row_vec[thread_id],
                                               row_for_col_vec[thread_id],
                                               u_vec[thread_id],
                                               v_vec[thread_id],
                                               workspaces[thread_id],
                                               difference_matrices_1[thread_id],
                                               difference_matrices_2[thread_id],
                                               start_indices);

            // Get the relevant sums                                   
            for (size_t i = 0; i < num_nodes; i++) {
                double v = 0.0;  // Sum of diff between i and all other nodes
                for (size_t j = 0; j < num_nodes; j++) {
                    if (i <= j) {
                        a = i;
                        b = j;
                    } else {
                        a = j;
                        b = i;
                    }
                    v += sim[start_indices[a] + (b - a)];
                }
                uniqueness_vec[task_id][i] = v / ((double) (num_nodes - 1));
            }

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
