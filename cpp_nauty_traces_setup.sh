#!/bin/bash

cd nauty27r4_modified
# Compile for C++ use. This changes the type used for thread-local storage
#   from _Thread_Local to thread_local .
./configure CC=g++ CFLAGS="-std=c++11" --enable-tls
make
