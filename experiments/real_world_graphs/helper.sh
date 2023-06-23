#!/bin/bash

# Karate

echo "Karate"

time nice -2 ../../executables/main -graph karate.txt \
                                    -o results/karate \
                                    -n_itr 250 \
                                    >> results/karate.txt
# Maayan Foodweb -- Directed

echo "Foodweb"

time nice -2 ../../executables/main -graph maayan-foodweb.txt \
                                    -nodes maayan-foodweb_nodes.txt
                                    -o results/foodweb \
                                    -n_itr 250 -d \
                                    >> results/foodweb.txt

# College Football Season 4 -- Undirected

echo "College Football Season 4"

time nice -2 ../../executables/main -graph season_4_undirected_edges.txt \
                                    -nodes season_4_undirected_nodes.txt \
                                    -o results/college_football_s4 \
                                    -n_itr 250 \
                                    >> results/college_football_s4.txt
# Political Blogs

echo "Political Blogs"

time nice -2 ../../executables/main -graph pol_blogs.txt \
                                    -nodes pol_blogs_nodes.txt
                                    -o results/pol_blogs \
                                    -n_itr 250 -d \
                                    >> results/pol_blogs.txt
# EU-Core Emails

echo "EU-Core Emails"

time nice -2 ../../executables/main -graph eucore.txt \
                                    -nodes eucore_nodes.txt
                                    -o results/eucore \
                                    -n_itr 250 -d \
                                    >> results/eucore.txt
# Cora Citations

echo "Cora Citations"

time nice -2 ../../executables/main -graph cora.txt \
                                    -nodes cora_nodes.txt
                                    -o results/cora \
                                    -n_itr 250 -d \
                                    >> results/cora.txt
