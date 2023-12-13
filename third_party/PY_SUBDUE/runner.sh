#!/bin/bash

echo "Karate"
rm testing/karate_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/karate.txt false

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 \
                      --overlap none --iterations 0 \
    graph_file.json >> testing/karate_output.txt

python3 post_process.py testing/karate_output.txt \
    ../experiments/real_world_graphs/karate.txt \
    /tmp/subdue_temp_file.txt testing/karate_noise.txt undirected



echo "Season 4"
rm testing/season_4_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/season_4_undirected_edges.txt false

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/season_4_output.txt

python3 post_process.py testing/season_4_output.txt \
   ../experiments/real_world_graphs/season_4_undirected_edges.txt \
    /tmp/subdue_temp_file.txt testing/season_4_noise.txt undirected



echo "Foodweb"
rm testing/foodweb_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/maayan-foodweb.txt true

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/foodweb_output.txt

python3 post_process.py testing/foodweb_output.txt \
   ../experiments/real_world_graphs/maayan-foodweb.txt \
    /tmp/subdue_temp_file.txt testing/foodweb_noise.txt directed



echo "Foodweb Undirected"
rm testing/foodweb_u_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/maayan-foodweb.txt false

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/foodweb_u_output.txt

python3 post_process.py testing/foodweb_u_output.txt \
   ../experiments/real_world_graphs/maayan-foodweb.txt \
    /tmp/subdue_temp_file.txt testing/foodweb_u_noise.txt undirected



echo "Pol Blogs"
rm testing/pol_blogs_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/pol_blogs.txt true

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/pol_blogs_output.txt

python3 post_process.py testing/pol_blogs_output.txt \
   ../experiments/real_world_graphs/pol_blogs.txt \
    /tmp/subdue_temp_file.txt testing/pol_blogs_noise.txt directed



echo "Pol Blogs Undirected"
rm testing/pol_blogs_u_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/pol_blogs.txt false

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/pol_blogs_u_output.txt

python3 post_process.py testing/pol_blogs_u_output.txt \
   ../experiments/real_world_graphs/pol_blogs.txt \
    /tmp/subdue_temp_file.txt testing/pol_blogs_u_noise.txt undirected



echo "EuCore"
rm testing/eucore_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/eucore.txt true

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/eucore_output.txt

python3 post_process.py testing/eucore_output.txt \
   ../experiments/real_world_graphs/eucore.txt \
    /tmp/subdue_temp_file.txt testing/eucore_noise.txt directed



echo "EuCore Undirected"
rm testing/eucore_u_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/eucore.txt false

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/eucore_u_output.txt

python3 post_process.py testing/eucore_u_output.txt \
   ../experiments/real_world_graphs/eucore.txt \
    /tmp/subdue_temp_file.txt testing/eucore_u_noise.txt undirected



echo "Cora"
rm testing/cora_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/cora.txt true

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/cora_output.txt

python3 post_process.py testing/cora_output.txt \
   ../experiments/real_world_graphs/cora.txt \
    /tmp/subdue_temp_file.txt testing/cora_noise.txt directed



echo "Cora Undirected"
rm testing/cora_u_output.txt

python3 prepare_graph.py \
    ../experiments/real_world_graphs/cora.txt false

python3 src/Subdue.py --minsize 2 --maxsize 10 --numbest 1 --overlap none \
    graph_file.json >> testing/cora_u_output.txt

python3 post_process.py testing/cora_u_output.txt \
   ../experiments/real_world_graphs/cora.txt \
    /tmp/subdue_temp_file.txt testing/cora_u_noise.txt undirected
