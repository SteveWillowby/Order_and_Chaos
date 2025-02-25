#!/bin/bash

# Johnson Graph

rm results/nsh_johnson_*

echo "Johnson 0.001"

time nice -2 ../../executables/SCHENO_ga -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.001 -noise- 0.001 -o results/nsh_johnson_001 \
                                    -no_sample_heuristic \
                                    >> results/nsh_johnson_001.txt

echo "Johnson 0.002 t1"

time nice -2 ../../executables/SCHENO_ga -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.002 -noise- 0.002 -o results/nsh_johnson_002_1 \
                                    -no_sample_heuristic \
                                    >> results/nsh_johnson_002_1.txt

echo "Johnson 0.002 t2"

time nice -2 ../../executables/SCHENO_ga -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.002 -noise- 0.002 -o results/nsh_johnson_002_2 \
                                    -no_sample_heuristic \
                                    >> results/nsh_johnson_002_2.txt

echo "Johnson 0.002 t3"

time nice -2 ../../executables/SCHENO_ga -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.002 -noise- 0.002 -o results/nsh_johnson_002_3 \
                                    -no_sample_heuristic \
                                    >> results/nsh_johnson_002_3.txt

echo "Johnson 0.003 t1"

time nice -2 ../../executables/SCHENO_ga -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.003 -noise- 0.003 -o results/nsh_johnson_003_1 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_johnson_003_1.txt

echo "Johnson 0.003 t2"

time nice -2 ../../executables/SCHENO_ga -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.003 -noise- 0.003 -o results/nsh_johnson_003_2 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_johnson_003_2.txt

echo "Johnson 0.003 t3"

time nice -2 ../../executables/SCHENO_ga -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.003 -noise- 0.003 -o results/nsh_johnson_003_3 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_johnson_003_3.txt

# Bin Tree

rm results/nsh_bin_tree_*

echo "Bin Tree 0.01"

time nice -2 ../../executables/SCHENO_ga -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.01 -noise- 0.01 -o results/nsh_bin_tree_01 \
                                    -no_sample_heuristic \
                                    >> results/nsh_bin_tree_01.txt

echo "Bin Tree 0.02"

time nice -2 ../../executables/SCHENO_ga -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.02 -noise- 0.02 -o results/nsh_bin_tree_02 \
                                    -no_sample_heuristic \
                                    >> results/nsh_bin_tree_02.txt

echo "Bin Tree 0.03"

time nice -2 ../../executables/SCHENO_ga -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.03 -noise- 0.03 -o results/nsh_bin_tree_03 \
                                    -no_sample_heuristic \
                                    >> results/nsh_bin_tree_03.txt

echo "Bin Tree 0.04"

time nice -2 ../../executables/SCHENO_ga -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.04 -noise- 0.04 -o results/nsh_bin_tree_04 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_bin_tree_04.txt

echo "Bin Tree 0.05"

time nice -2 ../../executables/SCHENO_ga -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.05 -noise- 0.05 -o results/nsh_bin_tree_05 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_bin_tree_05.txt

# Ring

rm results/nsh_ring_*

echo "Ring 0.01"

time nice -2 ../../executables/SCHENO_ga -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.01 -noise- 0.01 -o results/nsh_ring_01 \
                                    -no_sample_heuristic \
                                    >> results/nsh_ring_01.txt

echo "Ring 0.02"

time nice -2 ../../executables/SCHENO_ga -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.02 -noise- 0.02 -o results/nsh_ring_02 \
                                    -no_sample_heuristic \
                                    >> results/nsh_ring_02.txt

echo "Ring 0.03"

time nice -2 ../../executables/SCHENO_ga -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.03 -noise- 0.03 -o results/nsh_ring_03 \
                                    -no_sample_heuristic \
                                    >> results/nsh_ring_03.txt

echo "Ring 0.04"

time nice -2 ../../executables/SCHENO_ga -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.04 -noise- 0.04 -o results/nsh_ring_04 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_ring_04.txt

echo "Ring 0.05"

time nice -2 ../../executables/SCHENO_ga -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.05 -noise- 0.05 -o results/nsh_ring_05 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_ring_05.txt

# Wreath

rm results/nsh_wreath_*

echo "Wreath 0.0025"

time nice -2 ../../executables/SCHENO_ga -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.0025 -noise- 0.0025 -o results/nsh_wreath_0025 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_wreath_0025.txt

echo "Wreath 0.005"

time nice -2 ../../executables/SCHENO_ga -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.005 -noise- 0.005 -o results/nsh_wreath_005 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_wreath_005.txt

echo "Wreath 0.01"

time nice -2 ../../executables/SCHENO_ga -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.01 -noise- 0.01 -o results/nsh_wreath_01 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_wreath_01.txt

echo "Wreath 0.02"

time nice -2 ../../executables/SCHENO_ga -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.02 -noise- 0.02 -o results/nsh_wreath_02 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_wreath_02.txt

echo "Wreath 0.04"

time nice -2 ../../executables/SCHENO_ga -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.04 -noise- 0.04 -o results/nsh_wreath_04 \
                                    -n_itr 120 \
                                    -no_sample_heuristic \
                                    >> results/nsh_wreath_04.txt
