#!/bin/bash

echo ""
echo "######################################################"
echo "###     Compiling C Nauty/Traces and Utilities     ###"
echo "######################################################"

cd nauty27r4
./configure --enable-tls
make

cd ..

echo ""
echo ""
echo "######################################################"
echo "###      Compiling CPP Nauty/Traces Libraries      ###"
echo "######################################################"

cd nauty27r4_modified
# Compile for C++ use. This changes the type used for thread-local storage
#   from _Thread_Local to thread_local .
./configure CC=g++ CFLAGS="-std=c++11 -O4" --enable-tls
make

cd ..

echo ""
echo ""
echo "######################################################"
echo "###      Compiling CPP Wrappers & Utilities        ###"
echo "######################################################"

# Create the .o files
g++ -c nauty_traces.cpp file_utils.cpp debugging.cpp graph.cpp nt_partition.cpp sparse_graph.cpp nt_sparse_graph.cpp

# Create a library
ar rvs nt_wrappers_only.a nauty_traces.o file_utils.o debugging.o graph.o nt_partition.o sparse_graph.o nt_sparse_graph.o

# Combine nauty.a and nt_wrappers_only.a into a single library
ar -M <nt_wrappers.mri
