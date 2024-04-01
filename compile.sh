#!/bin/bash

# If using gcc instead of g++, add a link to -lstdc++
#   Removed -Werror due to a warning within the Nauty/Traces code itself
g++ -Wall -Wextra -o test_nt_code -std=c++11 -O4 test_nt_code.cpp nt_wrappers/nt_wrappers.a

g++ -Wall -Wextra -o minimal_nt_example -std=c++11 minimal_nt_example.cpp nt_wrappers/nt_wrappers.a

g++ -Wall -Wextra -o executables/graph_enumeration -std=c++11 -O4 graph_enumeration.cpp nt_wrappers/nt_wrappers.a

g++ -Wall -Wextra -o executables/mthread_scorer -std=c++11 -O4 test_threaded_scorer.cpp nt_wrappers/nt_wrappers.a edge_sampler.cpp scoring_function.cpp wl_measures.cpp Jonker_Volgenant/src/assign2DCBasic.c int_edge_sampler.cpp thread_pool_wl_sim.cpp thread_pool_scorer.cpp -lpthread

g++ -Wall -Wextra -o executables/SCHENO_score -std=c++11 -O4 SCHENO_score.cpp nt_wrappers/nt_wrappers.a edge_sampler.cpp scoring_function.cpp wl_measures.cpp Jonker_Volgenant/src/assign2DCBasic.c int_edge_sampler.cpp thread_pool_wl_sim.cpp thread_pool_scorer.cpp noise_prob_choice.cpp -lpthread

g++ -Wall -Wextra -o executables/SCHENO_ga -std=c++11 -O4 SCHENO_ga.cpp nt_wrappers/nt_wrappers.a edge_sampler.cpp scoring_function.cpp wl_measures.cpp Jonker_Volgenant/src/assign2DCBasic.c int_edge_sampler.cpp thread_pool_scorer.cpp thread_pool_wl_sim.cpp noise_prob_choice.cpp genetic_alg_search.cpp -lpthread
