#!/bin/bash

echo "4 Clique"
rm results/4_clique*
time nice -2 ../../executables/main -graph 4_clique.txt \
                                    -o results/4_clique \
                                    -n_itr 2 \
                                    >> results/4_clique.txt

echo "6 Cycle"
rm results/6_cycle*
time nice -2 ../../executables/main -graph 6_cycle.txt \
                                    -o results/6_cycle \
                                    -n_itr 5 \
                                    >> results/6_cycle.txt
echo "7 Cycle"
rm results/7_cycle*
time nice -2 ../../executables/main -graph 7_cycle.txt \
                                    -o results/7_cycle \
                                    -n_itr 15 \
                                    >> results/7_cycle.txt

echo "7 Star"
rm results/7_star*
time nice -2 ../../executables/main -graph 7_star.txt \
                                    -o results/7_star \
                                    -n_itr 15 \
                                    >> results/7_star.txt

echo "3 Triangles"
rm results/3_triangles*
time nice -2 ../../executables/main -graph 3_triangles.txt \
                                    -o results/3_triangles \
                                    -n_itr 20 \
                                    >> results/3_triangles.txt

echo "'Human'"
rm results/human*
time nice -2 ../../executables/main -graph human.txt \
                                    -o results/human \
                                    -n_itr 20 \
                                    >> results/human.txt
