#!/bin/bash

# If using gcc instead of g++, add a link to -lstdc++
g++ -Wall -Werror -Wextra -o executables/nt_graph_test -std=c++11 nt_graph_test.cpp nt_partition.cpp graph.cpp sparse_graph.cpp nt_sparse_graph.cpp traces.cpp nauty27r4/nauty.a
