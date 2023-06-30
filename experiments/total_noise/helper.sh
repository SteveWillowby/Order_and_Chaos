#!/bin/bash

# 50 nodes

echo "50 Nodes t1"

time nice -2 ../../executables/main -graph clique_50.txt \
                                    -noise- 0.5 -o results/run_50_01 \
                                    -n_itr 150 \
                                    >> results/run_50_01.txt

echo "50 Nodes t2"

time nice -2 ../../executables/main -graph clique_50.txt \
                                    -noise- 0.5 -o results/run_50_02 \
                                    -n_itr 150 \
                                    >> results/run_50_02.txt

echo "50 Nodes t3"

time nice -2 ../../executables/main -graph clique_50.txt \
                                    -noise- 0.5 -o results/run_50_03 \
                                    -n_itr 150 \
                                    >> results/run_50_03.txt


# 100 nodes

echo "100 Nodes t1"

time nice -2 ../../executables/main -graph clique_100.txt \
                                    -noise- 0.5 -o results/run_100_01 \
                                    -n_itr 150 \
                                    >> results/run_100_01.txt

echo "100 Nodes t2"

time nice -2 ../../executables/main -graph clique_100.txt \
                                    -noise- 0.5 -o results/run_100_02 \
                                    -n_itr 150 \
                                    >> results/run_100_02.txt

echo "100 Nodes t3"

time nice -2 ../../executables/main -graph clique_100.txt \
                                    -noise- 0.5 -o results/run_100_03 \
                                    -n_itr 150 \
                                    >> results/run_100_03.txt

# 150 nodes

echo "150 Nodes t1"

time nice -2 ../../executables/main -graph clique_150.txt \
                                    -noise- 0.5 -o results/run_150_01 \
                                    -n_itr 150 \
                                    >> results/run_150_01.txt

echo "150 Nodes t2"

time nice -2 ../../executables/main -graph clique_150.txt \
                                    -noise- 0.5 -o results/run_150_02 \
                                    -n_itr 150 \
                                    >> results/run_150_02.txt

echo "150 Nodes t3"

time nice -2 ../../executables/main -graph clique_150.txt \
                                    -noise- 0.5 -o results/run_150_03 \
                                    -n_itr 150 \
                                    >> results/run_150_03.txt

# 200 nodes

echo "200 Nodes t1"

time nice -2 ../../executables/main -graph clique_200.txt \
                                    -noise- 0.5 -o results/run_200_01 \
                                    -n_itr 150 \
                                    >> results/run_200_01.txt

echo "200 Nodes t2"

time nice -2 ../../executables/main -graph clique_200.txt \
                                    -noise- 0.5 -o results/run_200_02 \
                                    -n_itr 150 \
                                    >> results/run_200_02.txt

echo "200 Nodes t3"

time nice -2 ../../executables/main -graph clique_200.txt \
                                    -noise- 0.5 -o results/run_200_03 \
                                    -n_itr 150 \
                                    >> results/run_200_03.txt
