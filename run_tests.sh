#!/bin/bash

rm test_results/bin_tree_127.txt
nice -2 executables/main_bin_tree_127 >> test_results/bin_tree_127.txt
rm test_results/ring_128.txt
nice -2 executables/main_ring_128 >> test_results/ring_128.txt
rm test_results/wreath_128.txt
nice -2 executables/main_wreath_128 >> test_results/wreath_128.txt
rm test_results/johnson_120.txt
nice -2 executables/main_johnson_120 >> test_results/johnson_120.txt

rm test_results/season_4_undir_p.5.txt
nice -2 executables/main_season_4_undir >> test_results/season_4_undir_p.5.txt
rm test_results/season_4_dir_p.5.txt
nice -2 executables/main_season_4_dir >> test_results/season_4_dir_p.5.txt
rm test_results/season_3_undir_p.5.txt
nice -2 executables/main_season_3_undir >> test_results/season_3_undir_p.5.txt
rm test_results/season_3_dir_p.5.txt
nice -2 executables/main_season_3_dir >> test_results/season_3_dir_p.5.txt
