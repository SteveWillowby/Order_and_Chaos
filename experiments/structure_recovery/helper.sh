#!/bin/bash

# Johnson Graph

echo "Johnson 0.001"

time nice -2 ../../executables/main -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.001 -noise- 0.001 -o results/johnson_001 \
                                    >> results/johnson_001.txt

echo "Johnson 0.002 t1"

time nice -2 ../../executables/main -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.002 -noise- 0.002 -o results/johnson_002_1 \
                                    >> results/johnson_002_1.txt

echo "Johnson 0.002 t2"

time nice -2 ../../executables/main -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.002 -noise- 0.002 -o results/johnson_002_2 \
                                    >> results/johnson_002_2.txt

echo "Johnson 0.002 t3"

time nice -2 ../../executables/main -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.002 -noise- 0.002 -o results/johnson_002_3 \
                                    >> results/johnson_002_3.txt

echo "Johnson 0.003 t1"

time nice -2 ../../executables/main -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.003 -noise- 0.003 -o results/johnson_003_1 \
                                    -n_itr 120 \
                                    >> results/johnson_003_1.txt

echo "Johnson 0.003 t2"

time nice -2 ../../executables/main -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.003 -noise- 0.003 -o results/johnson_003_2 \
                                    -n_itr 120 \
                                    >> results/johnson_003_2.txt

echo "Johnson 0.003 t3"

time nice -2 ../../executables/main -graph johnson_10_3_120_edges.txt \
                                    -nodes johnson_10_3_120_nodes.txt \
                                    -noise+ 0.003 -noise- 0.003 -o results/johnson_003_3 \
                                    -n_itr 120 \
                                    >> results/johnson_003_3.txt

# Bin Tree

echo "Bin Tree 0.01"

time nice -2 ../../executables/main -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.01 -noise- 0.01 -o results/bin_tree_01 \
                                    >> results/bin_tree_01.txt

echo "Bin Tree 0.02"

time nice -2 ../../executables/main -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.02 -noise- 0.02 -o results/bin_tree_02 \
                                    >> results/bin_tree_02.txt

echo "Bin Tree 0.03"

time nice -2 ../../executables/main -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.03 -noise- 0.03 -o results/bin_tree_03 \
                                    >> results/bin_tree_03.txt

echo "Bin Tree 0.04"

time nice -2 ../../executables/main -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.04 -noise- 0.04 -o results/bin_tree_04 \
                                    -n_itr 120 \
                                    >> results/bin_tree_04.txt

echo "Bin Tree 0.05"

time nice -2 ../../executables/main -graph binary_tree_127_edges.txt \
                                    -nodes binary_tree_127_nodes.txt \
                                    -noise+ 0.05 -noise- 0.05 -o results/bin_tree_05 \
                                    -n_itr 120 \
                                    >> results/bin_tree_05.txt

# Ring

echo "Ring 0.01"

time nice -2 ../../executables/main -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.01 -noise- 0.01 -o results/ring_01 \
                                    >> results/ring_01.txt

echo "Ring 0.02"

time nice -2 ../../executables/main -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.02 -noise- 0.02 -o results/ring_02 \
                                    >> results/ring_02.txt

echo "Ring 0.03"

time nice -2 ../../executables/main -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.03 -noise- 0.03 -o results/ring_03 \
                                    >> results/ring_03.txt

echo "Ring 0.04"

time nice -2 ../../executables/main -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.04 -noise- 0.04 -o results/ring_04 \
                                    -n_itr 120 \
                                    >> results/ring_04.txt

echo "Ring 0.05"

time nice -2 ../../executables/main -graph ring_128_edges.txt \
                                    -nodes ring_128_nodes.txt \
                                    -noise+ 0.05 -noise- 0.05 -o results/ring_05 \
                                    -n_itr 120 \
                                    >> results/ring_05.txt

# Wreath

echo "Wreath 0.0025"

time nice -2 ../../executables/main -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.0025 -noise- 0.0025 -o results/wreath_0025 \
                                    -n_itr 120 \
                                    >> results/wreath_0025.txt

echo "Wreath 0.005"

time nice -2 ../../executables/main -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.005 -noise- 0.005 -o results/wreath_005 \
                                    -n_itr 120 \
                                    >> results/wreath_005.txt

echo "Wreath 0.01"

time nice -2 ../../executables/main -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.01 -noise- 0.01 -o results/wreath_01 \
                                    -n_itr 120 \
                                    >> results/wreath_01.txt

echo "Wreath 0.02"

time nice -2 ../../executables/main -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.02 -noise- 0.02 -o results/wreath_02 \
                                    -n_itr 120 \
                                    >> results/wreath_02.txt

echo "Wreath 0.04"

time nice -2 ../../executables/main -graph wreath_d7_128_edges.txt \
                                    -nodes wreath_d7_128_nodes.txt \
                                    -noise+ 0.04 -noise- 0.04 -o results/wreath_04 \
                                    -n_itr 120 \
                                    >> results/wreath_04.txt
