#!/bin/bash


# Karate

echo "Karate"
rm results/karate_ga*
time nice -2 ../executables/SCHENO_ga -graph results/karate_og_graph.txt \
                                 -nodes results/karate_og_nodes.txt \
                                 -seed  results/karate_og_noise.txt \
                                 -o results/karate_ga \
                                 -n_itr 250 -u \
                                 >> results/karate_ga.txt

# Maayan Foodweb

echo "Undirected Foodweb"
rm results/foodweb_ga*
time nice -2 ../executables/SCHENO_ga -graph results/foodweb_og_graph.txt \
                                 -nodes results/foodweb_og_nodes.txt \
                                 -seed  results/foodweb_og_noise.txt \
                                 -o results/foodweb_ga \
                                 -n_itr 250 -u \
                                 >> results/foodweb_ga.txt

# College Football Season 4 -- Undirected

echo "College Football Season 4"
rm results/season_4_ga*
time nice -2 ../executables/SCHENO_ga -graph results/season_4_og_graph.txt \
                                 -nodes results/season_4_og_nodes.txt \
                                 -seed  results/season_4_og_noise.txt \
                                 -o results/season_4_ga \
                                 -n_itr 250 -u \
                                 >> results/season_4_ga.txt

# Political Blogs

echo "Undirected Political Blogs"
rm results/pol_blogs_ga*
time nice -2 ../executables/SCHENO_ga -graph results/pol_blogs_og_graph.txt \
                                 -nodes results/pol_blogs_og_nodes.txt \
                                 -seed  results/pol_blogs_og_noise.txt \
                                 -o results/pol_blogs_ga \
                                 -n_itr 250 -u \
                                 >> results/pol_blogs_ga.txt
# EU-Core Emails

echo "Undirected EU-Core Emails"
rm results/eucore_ga*
time nice -2 ../executables/SCHENO_ga -graph results/eucore_og_graph.txt \
                                 -nodes results/eucore_og_nodes.txt \
                                 -seed  results/eucore_og_noise.txt \
                                 -o results/eucore_ga \
                                 -n_itr 250 -u \
                                 >> results/eucore_ga.txt
# Cora Citations

echo "Undirected Cora Citations"
rm results/cora_ga*
time nice -2 ../executables/SCHENO_ga -graph results/cora_og_graph.txt \
                                 -nodes results/cora_og_nodes.txt \
                                 -seed  results/cora_og_noise.txt \
                                 -o results/cora_ga \
                                 -n_itr 250 -u \
                                 >> results/cora_ga.txt
