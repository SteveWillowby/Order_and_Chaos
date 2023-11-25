#!/bin/bash

python3 prep_graph.py $1
matlab -r run_structureDiscovery

python2.7 MDL/greedySearch_nStop.py DATA/input_graph.txt DATA/input_graph_orderedALL.model >/dev/null 2>&1
mv heuristic* DATA/

python3 post_process.py DATA/input_graph_orderedALL.model DATA/input_graph.txt \
        DATA/heuristicSelection_nStop_ALL_input_graph_orderedALL.model $2

# unweighted_graph='DATA/input_graph.txt'
# model='DATA/input_graph_orderedALL.model'
# python2.7 MDL/greedySearch_nStop.py $unweighted_graph $model >/dev/null 2>&1
# mv heuristic* DATA/
