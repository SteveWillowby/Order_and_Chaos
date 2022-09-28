#!/bin/bash

# If using gcc instead of g++, add a link to -lstdc++
g++ -Wall -Werror -Wextra -o executables/nt_graph_test -std=c++11 nt_graph_test.cpp graph.cpp augmented_multimap.cpp sparse_graph.cpp nt_sparse_graph.cpp nauty27r4/nauty.a
