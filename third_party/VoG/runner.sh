#!/bin/bash

./runner_helper.sh \
    ../experiments/real_world_graphs/karate.txt \
    results/karate

./runner_helper.sh \
    ../experiments/real_world_graphs/season_4_undirected_edges.txt \
    results/season_4

./runner_helper.sh \
    ../experiments/real_world_graphs/cora.txt \
    results/cora

./runner_helper.sh \
    ../experiments/real_world_graphs/eucore.txt \
    results/eucore

./runner_helper.sh \
    ../experiments/real_world_graphs/maayan-foodweb.txt \
    results/foodweb

./runner_helper.sh \
    ../experiments/real_world_graphs/pol_blogs.txt \
    results/pol_blogs
