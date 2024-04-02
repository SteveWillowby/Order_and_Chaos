#!/bin/bash

# 50 nodes

echo "50 Nodes t1"

# n at 50 --> p = 0.36284
#   1 - p = 0.63716
time nice -2 ../../executables/SCHENO_ga -graph clique_50.txt \
                                    -noise- 0.63716 -o results/run_50_01 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_50_01.txt

echo "50 Nodes t2"

time nice -2 ../../executables/SCHENO_ga -graph clique_50.txt \
                                    -noise- 0.63716 -o results/run_50_02 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_50_02.txt

echo "50 Nodes t3"

time nice -2 ../../executables/SCHENO_ga -graph clique_50.txt \
                                    -noise- 0.63716 -o results/run_50_03 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_50_03.txt


# 100 nodes

echo "100 Nodes t1"

# n at 100 --> p = 0.420843
#   1 - p = 0.579157
time nice -2 ../../executables/SCHENO_ga -graph clique_100.txt \
                                    -noise- 0.579157 -o results/run_100_01 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_100_01.txt

echo "100 Nodes t2"

time nice -2 ../../executables/SCHENO_ga -graph clique_100.txt \
                                    -noise- 0.579157 -o results/run_100_02 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_100_02.txt

echo "100 Nodes t3"

time nice -2 ../../executables/SCHENO_ga -graph clique_100.txt \
                                    -noise- 0.579157 -o results/run_100_03 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_100_03.txt

# 150 nodes

echo "150 Nodes t1"

# n at 100 --> p = 0.44282
#   1 - p = 0.55718
time nice -2 ../../executables/SCHENO_ga -graph clique_150.txt \
                                    -noise- 0.55718 -o results/run_150_01 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_150_01.txt

echo "150 Nodes t2"

time nice -2 ../../executables/SCHENO_ga -graph clique_150.txt \
                                    -noise- 0.55718 -o results/run_150_02 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_150_02.txt

echo "150 Nodes t3"

time nice -2 ../../executables/SCHENO_ga -graph clique_150.txt \
                                    -noise- 0.55718 -o results/run_150_03 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_150_03.txt

# 200 nodes

echo "200 Nodes t1"

# n at 200 --> p = 0.454684
#   1 - p = 0.545316
time nice -2 ../../executables/SCHENO_ga -graph clique_200.txt \
                                    -noise- 0.545316 -o results/run_200_01 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_200_01.txt

echo "200 Nodes t2"

time nice -2 ../../executables/SCHENO_ga -graph clique_200.txt \
                                    -noise- 0.545316 -o results/run_200_02 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_200_02.txt

echo "200 Nodes t3"

time nice -2 ../../executables/SCHENO_ga -graph clique_200.txt \
                                    -noise- 0.545316 -o results/run_200_03 \
                                    -n_itr 150 \
                                    -use_g_as_ses \
                                    >> results/run_200_03.txt
