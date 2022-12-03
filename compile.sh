#!/bin/bash

# If using gcc instead of g++, add a link to -lstdc++
#   Removed -Werror due to a warning within the Nauty/Traces code itself
g++ -Wall -Wextra -o executables/nt_graph_test -std=c++11 -O4 testing.cpp nt_partition.cpp graph.cpp sparse_graph.cpp nt_sparse_graph.cpp debugging.cpp nauty_traces.cpp file_utils.cpp nauty27r4_modified/nauty.a
