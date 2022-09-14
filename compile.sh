#!/bin/bash

# For some reason gcc has issues linking to the stdexcept library on this machine.
g++ -o executables/nt_graph_test -std=c++11 nt_graph_test.cpp graph.cpp sparse_graph.cpp nt_sparse_graph.cpp nauty27r4/nauty.a
