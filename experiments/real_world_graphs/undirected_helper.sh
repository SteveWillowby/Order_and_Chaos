#!/bin/bash

# Karate

# echo "Karate"
# rm results/karate*
# time nice -2 ../../executables/main -graph karate.txt \
#                                     -o results/karate \
#                                     -n_itr 250 \
#                                     >> results/karate.txt

# Maayan Foodweb

echo "Undirected Foodweb"
rm results/undir_foodweb*
time nice -2 ../../executables/main -graph maayan-foodweb.txt \
                                    -nodes maayan-foodweb_nodes.txt \
                                    -o results/undir_foodweb \
                                    -n_itr 250 -u \
                                    >> results/undir_foodweb.txt

# College Football Season 4 -- Undirected

# echo "College Football Season 4"
# rm results/college_football_s4*
# time nice -2 ../../executables/main -graph season_4_undirected_edges.txt \
#                                     -nodes season_4_undirected_nodes.txt \
#                                     -o results/college_football_s4 \
#                                     -n_itr 250 \
#                                     >> results/college_football_s4.txt

# Political Blogs

echo "Undirected Political Blogs"
rm results/undir_pol_blogs*
time nice -2 ../../executables/main -graph pol_blogs.txt \
                                    -nodes pol_blogs_nodes.txt \
                                    -o results/undir_pol_blogs \
                                    -n_itr 250 -u \
                                    >> results/undir_pol_blogs.txt
# EU-Core Emails

echo "Undirected EU-Core Emails"
rm results/undir_eucore*
time nice -2 ../../executables/main -graph eucore.txt \
                                    -nodes eucore_nodes.txt \
                                    -o results/undir_eucore \
                                    -n_itr 250 -u \
                                    >> results/undir_eucore.txt
# Cora Citations

echo "Undirected Cora Citations"
rm results/undir_cora*
time nice -2 ../../executables/main -graph cora.txt \
                                    -nodes cora_nodes.txt \
                                    -o results/undir_cora \
                                    -n_itr 250 -u \
                                    >> results/undir_cora.txt
