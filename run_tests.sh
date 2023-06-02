#!/bin/bash

# rm test_results/bin_tree_127.txt
# time nice -2 executables/main_bin_tree_127 >> test_results/bin_tree_127.txt
# rm test_results/ring_128.txt
# time nice -2 executables/main_ring_128 >> test_results/ring_128.txt
# rm test_results/wreath_128.txt
# time nice -2 executables/main_wreath_128 >> test_results/wreath_128.txt
# rm test_results/johnson_120.txt
# time nice -2 executables/main_johnson_120 >> test_results/johnson_120.txt

rm test_results/season_4_undir_p.5.txt
time nice -2 executables/main_season_4_undir >> test_results/season_4_undir_p.5.txt
rm test_results/season_4_dir_p.5.txt
time nice -2 executables/main_season_4_dir >> test_results/season_4_dir_p.5.txt
# rm test_results/season_3_undir_p.5.txt
# time nice -2 executables/main_season_3_undir >> test_results/season_3_undir_p.5.txt
# rm test_results/season_3_dir_p.5.txt
# time nice -2 executables/main_season_3_dir >> test_results/season_3_dir_p.5.txt

rm test_results/jazz_colab_undir.txt
time nice -2 executables/main_jazz_undir >> test_results/jazz_undir.txt
rm test_results/c_elegans_dir.txt
time nice -2 executables/main_c_elegans_dir >> test_results/c_elegans_dir.txt
rm test_results/species_brain_dir.txt
time nice -2 executables/main_species_brain_dir >> test_results/species_brain_dir.txt
