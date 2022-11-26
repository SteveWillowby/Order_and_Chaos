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
