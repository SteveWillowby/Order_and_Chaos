#include "nt_wrappers/nauty_traces.h"

#include<condition_variable>
#include<memory>
#include<mutex>
#include<thread>
#include<unordered_set>
#include<vector>

#ifndef SCHENO__THREAD_POOL_WL_SIM_H
#define SCHENO__THREAD_POOL_WL_SIM_H

class ThreadPoolWLSim {
public:
    // `nt` is the number of threads to use. If 0 is passed, the code will
    //      use the default of std::thread::hardware_concurrency()
    //
    // `n` is the number of nodes
    ThreadPoolWLSim(size_t nt, size_t n);

    ~ThreadPoolWLSim();

    // IMPORTANT: The vector returned may have extra elements at the end.
    //  Thus you must ignore get_fuzzy_orbit_sizes().size().
    //
    // A task consists of a graph and a node ID (for node-centric computation).
    //  If node-centric is not desired, pass -1 as the node.
    //
    // DO NOT MODIFY THE RETURNED VECTOR
    const std::vector<std::vector<double>> *get_fuzzy_orbit_sizes(
            const std::vector<std::pair<const Graph*, size_t>>* tasks);

    void terminate();

protected:
    const size_t num_threads;
    const size_t num_nodes;

    const std::vector<std::pair<const Graph*, size_t>>* task_vec;
    std::vector<std::thread> pool;

    const Graph* g_;

    // Editable data for heuristic
    size_t* start_indices;
    std::vector<double*> u_vec;
    std::vector<double*> v_vec;
    std::vector<double*> cost_matrices;
    std::vector<double*> difference_matrices_1;
    std::vector<double*> difference_matrices_2;
    std::vector<ptrdiff_t*> col_for_row_vec;
    std::vector<ptrdiff_t*> row_for_col_vec;
    std::vector<void*> workspaces;
    std::vector<Coloring<int>> node_colorings;

    std::vector<std::vector<double>> orbits_vec;

    bool terminate_pool;

    size_t tasks_begun, num_tasks;
    size_t threads_working, threads_launched;

    // Used by running workers to select the next task
    std::mutex m_worker_queue;
    // Used by workers (and get_scores()) to begin and end working on a set of
    //  tasks
    std::mutex m_worker_meta;
    // Used by workers (and get_scores()) to say when all current tasks are done
    std::mutex m_tasks_done;
    // Used to coordinate thread launches so that IDs are correct.
    std::mutex m_launch;
    // Used in conjunction with m_worker_meta
    std::condition_variable worker_wait_signal;
    // Used in conjunction with m_tasks_done
    std::condition_variable done_wait_signal;
    // Used in conjunction with m_launch
    std::condition_variable launch_wait_signal;

    void run();
};

#endif
