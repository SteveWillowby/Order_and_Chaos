#!/bin/bash

# If using gcc instead of g++, add a link to -lstdc++
#   Removed -Werror due to a warning within the Nauty/Traces code itself
g++ -Wall -Wextra -o executables/nt_graph_test -std=c++11 -O4 testing.cpp nt_partition.cpp graph.cpp sparse_graph.cpp nt_sparse_graph.cpp debugging.cpp nauty_traces.cpp file_utils.cpp nauty27r4_modified/nauty.a

g++ -Wall -Wextra -o executables/mthread_scorer -std=c++11 -O4 test_threaded_scorer.cpp nt_partition.cpp graph.cpp sparse_graph.cpp nt_sparse_graph.cpp nauty_traces.cpp file_utils.cpp edge_sampler.cpp scoring_function.cpp scoring_heuristic.cpp thread_pool_scorer.cpp nauty27r4_modified/nauty.a -lpthread

g++ -Wall -Wextra -o executables/main -std=c++11 -O4 main.cpp nt_partition.cpp graph.cpp sparse_graph.cpp nt_sparse_graph.cpp nauty_traces.cpp file_utils.cpp edge_sampler.cpp scoring_function.cpp scoring_heuristic.cpp thread_pool_scorer.cpp noise_prob_choice.cpp genetic_alg_search.cpp simulated_annealing_search.cpp nauty27r4_modified/nauty.a -lpthread
