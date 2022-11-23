#!/bin/bash

cd nauty27r4_modified
# Compile for C++ use. This changes the type used for thread-local storage
#   from _Thread_Local to thread_local .
./configure CC=g++ CFLAGS="-std=c++11" --enable-tls
make

echo "############################################################################"
echo "# If the Compilation Failed at the Compiling of dreadnaut.c, you are fine. #"
echo "############################################################################"
