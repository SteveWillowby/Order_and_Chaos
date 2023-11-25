#!/bin/bash

rm results/*

python3 file_prep.py \
        ../experiments/real_world_graphs/season_4_undirected_edges.txt u
python3 crossspot.py prepared.csv >> results/cfb_season_4.txt
python3 crossspot.py line_graph.csv >> results/cfb_season_4_lg.txt

python3 file_prep.py ../experiments/real_world_graphs/cora.txt
python3 crossspot.py prepared.csv >> results/cora.txt
python3 crossspot.py line_graph.csv >> results/cora_lg.txt

python3 file_prep.py ../experiments/real_world_graphs/eucore.txt
python3 crossspot.py prepared.csv >> results/eucore.txt
python3 crossspot.py line_graph.csv >> results/eucore_lg.txt

python3 file_prep.py ../experiments/real_world_graphs/maayan-foodweb.txt
python3 crossspot.py prepared.csv >> results/foodweb.txt
python3 crossspot.py line_graph.csv >> results/foodweb_lg.txt

python3 file_prep.py ../experiments/real_world_graphs/karate.txt u
python3 crossspot.py prepared.csv >> results/karate.txt
python3 crossspot.py line_graph.csv >> results/karate_lg.txt

python3 file_prep.py ../experiments/real_world_graphs/pol_blogs.txt
python3 crossspot.py prepared.csv >> results/pol_blogs.txt
python3 crossspot.py line_graph.csv >> results/pol_blogs_lg.txt
