#!/bin/bash

# Karate

echo "Karate"
rm results/nsh_karate*
time nice -2 ../../executables/SCHENO_ga -graph karate.txt \
                                    -o results/nsh_karate \
                                    -n_itr 250 \
                                    -no_sample_heuristic \
                                    >> results/nsh_karate.txt

# College Football Season 4 -- Undirected

echo "College Football Season 4"
rm results/nsh_college_football_s4*
time nice -2 ../../executables/SCHENO_ga -graph season_4_undirected_edges.txt \
                                    -nodes season_4_undirected_nodes.txt \
                                    -o results/nsh_college_football_s4 \
                                    -n_itr 250 \
                                    -no_sample_heuristic \
                                    >> results/nsh_college_football_s4.txt
# Political Blogs

echo "Political Blogs"
rm results/nsh_pol_blogs*
time nice -2 ../../executables/SCHENO_ga -graph pol_blogs.txt \
                                    -nodes pol_blogs_nodes.txt \
                                    -o results/nsh_pol_blogs \
                                    -n_itr 250 -d \
                                    -no_sample_heuristic \
                                    >> results/nsh_pol_blogs.txt
# EU-Core Emails

echo "EU-Core Emails"
rm results/nsh_eucore*
time nice -2 ../../executables/SCHENO_ga -graph eucore.txt \
                                    -nodes eucore_nodes.txt \
                                    -o results/nsh_eucore \
                                    -n_itr 250 -d \
                                    -no_sample_heuristic \
                                    >> results/nsh_eucore.txt
# Cora Citations

echo "Cora Citations"
rm results/nsh_cora*
time nice -2 ../../executables/SCHENO_ga -graph cora.txt \
                                    -nodes cora_nodes.txt \
                                    -o results/nsh_cora \
                                    -n_itr 250 -d \
                                    -no_sample_heuristic \
                                    >> results/nsh_cora.txt

# Maayan Foodweb -- Directed

echo "Foodweb"
rm results/nsh_foodweb*
time nice -2 ../../executables/SCHENO_ga -graph maayan-foodweb.txt \
                                    -nodes maayan-foodweb_nodes.txt \
                                    -o results/nsh_foodweb \
                                    -n_itr 250 -d \
                                    -no_sample_heuristic \
                                    >> results/nsh_foodweb.txt
