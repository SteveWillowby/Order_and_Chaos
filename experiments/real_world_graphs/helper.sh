#!/bin/bash

# Karate

# echo "Karate"
# rm results/karate*
# time nice -2 ../../executables/SCHENO_ga -graph karate.txt \
#                                     -o results/karate \
#                                     -n_itr 250 \
#                                     >> results/karate.txt
# Maayan Foodweb -- Directed

# echo "Foodweb"
# rm results/foodweb*
# time nice -2 ../../executables/SCHENO_ga -graph maayan-foodweb.txt \
#                                     -nodes maayan-foodweb_nodes.txt \
#                                     -o results/foodweb \
#                                     -n_itr 250 -d \
#                                     >> results/foodweb.txt

# College Football Season 4 -- Undirected

echo "College Football Season 4"
rm results/college_football_s4*
time nice -2 ../../executables/SCHENO_ga -graph season_4_undirected_edges.txt \
                                    -nodes season_4_undirected_nodes.txt \
                                    -o results/college_football_s4 \
                                    -n_itr 250 \
                                    >> results/college_football_s4.txt
# Political Blogs

echo "Political Blogs"
rm results/pol_blogs*
time nice -2 ../../executables/SCHENO_ga -graph pol_blogs.txt \
                                    -nodes pol_blogs_nodes.txt \
                                    -o results/pol_blogs \
                                    -n_itr 250 -d \
                                    >> results/pol_blogs.txt
# EU-Core Emails

echo "EU-Core Emails"
rm results/eucore*
time nice -2 ../../executables/SCHENO_ga -graph eucore.txt \
                                    -nodes eucore_nodes.txt \
                                    -o results/eucore \
                                    -n_itr 250 -d \
                                    >> results/eucore.txt
# Cora Citations

echo "Cora Citations"
rm results/cora*
time nice -2 ../../executables/SCHENO_ga -graph cora.txt \
                                    -nodes cora_nodes.txt \
                                    -o results/cora \
                                    -n_itr 250 -d \
                                    >> results/cora.txt
