#!/bin/bash

python3 prepare_graph.py \
    ../../experiments/real_world_graphs/karate.txt false

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 \
                      --overlap none --iterations 0 \
    graph_file.json >> testing/karate_iterated_output.txt
