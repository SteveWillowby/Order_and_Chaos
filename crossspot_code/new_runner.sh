#!/bin/bash


python3 file_prep.py \
    ../experiments/real_world_graphs/season_4_undirected_edges.txt u
python3 crossspot.py prepared.csv > results/cfb_season_4.txt
python3 post_process.py \
    ../experiments/real_world_graphs/season_4_undirected_edges.txt \
    results/cfb_season_4.txt u > results/cfb_season_4_score.txt
echo "No Edges" >> results/cfb_season_4_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/season_4_undirected_edges.txt \
    -nodes nodes.txt -u \
    -edges no_edges.txt >> results/cfb_season_4_score.txt
echo "CrossSpot Edges" >> results/cfb_season_4_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/season_4_undirected_edges.txt \
    -nodes nodes.txt -u \
    -edges chosen_edges.txt >> results/cfb_season_4_score.txt
cat chosen_edges.txt >> results/cfb_season_4_score.txt


python3 file_prep.py ../experiments/real_world_graphs/karate.txt u
python3 crossspot.py prepared.csv > results/karate.txt
python3 post_process.py \
    ../experiments/real_world_graphs/karate.txt \
    results/karate.txt u > results/karate_score.txt
echo "No Edges" >> results/karate_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/karate.txt \
    -nodes nodes.txt -u \
    -edges no_edges.txt >> results/karate_score.txt
echo "CrossSpot Edges" >> results/karate_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/karate.txt \
    -nodes nodes.txt -u \
    -edges chosen_edges.txt >> results/karate_score.txt
cat chosen_edges.txt >> results/karate_score.txt


python3 file_prep.py ../experiments/real_world_graphs/cora.txt d
python3 crossspot.py prepared.csv > results/cora.txt
python3 post_process.py \
    ../experiments/real_world_graphs/cora.txt \
    results/cora.txt d > results/cora_score.txt
echo "No Edges" >> results/cora_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/cora.txt \
    -nodes nodes.txt -d \
    -edges no_edges.txt >> results/cora_score.txt
echo "CrossSpot Edges" >> results/cora_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/cora.txt \
    -nodes nodes.txt -d \
    -edges chosen_edges.txt >> results/cora_score.txt
cat chosen_edges.txt >> results/cora_score.txt


python3 file_prep.py ../experiments/real_world_graphs/eucore.txt d
python3 crossspot.py prepared.csv > results/eucore.txt
python3 post_process.py \
    ../experiments/real_world_graphs/eucore.txt \
    results/eucore.txt d > results/eucore_score.txt
echo "No Edges" >> results/eucore_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/eucore.txt \
    -nodes nodes.txt -d \
    -edges no_edges.txt >> results/eucore_score.txt
echo "CrossSpot Edges" >> results/eucore_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/eucore.txt \
    -nodes nodes.txt -d \
    -edges chosen_edges.txt >> results/eucore_score.txt
cat chosen_edges.txt >> results/eucore_score.txt


python3 file_prep.py ../experiments/real_world_graphs/maayan-foodweb.txt d
python3 crossspot.py prepared.csv > results/foodweb.txt
python3 post_process.py \
    ../experiments/real_world_graphs/maayan-foodweb.txt \
    results/foodweb.txt d > results/foodweb_score.txt
echo "No Edges" >> results/foodweb_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/maayan-foodweb.txt \
    -nodes nodes.txt -d \
    -edges no_edges.txt >> results/foodweb_score.txt
echo "CrossSpot Edges" >> results/foodweb_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/maayan-foodweb.txt \
    -nodes nodes.txt -d \
    -edges chosen_edges.txt >> results/foodweb_score.txt
cat chosen_edges.txt >> results/foodweb_score.txt


python3 file_prep.py ../experiments/real_world_graphs/pol_blogs.txt d
python3 crossspot.py prepared.csv > results/pol_blogs.txt
python3 post_process.py \
    ../experiments/real_world_graphs/pol_blogs.txt \
    results/pol_blogs.txt d > results/pol_blogs_score.txt
echo "No Edges" >> results/pol_blogs_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/pol_blogs.txt \
    -nodes nodes.txt -d \
    -edges no_edges.txt >> results/pol_blogs_score.txt
echo "CrossSpot Edges" >> results/pol_blogs_score.txt
../executables/score_only \
    -graph ../experiments/real_world_graphs/pol_blogs.txt \
    -nodes nodes.txt -d \
    -edges chosen_edges.txt >> results/pol_blogs_score.txt
cat chosen_edges.txt >> results/pol_blogs_score.txt
