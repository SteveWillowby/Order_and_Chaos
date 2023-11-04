#!/bin/bash

echo ""
echo "Karate Undirected - All Noise"
executables/score_only -graph experiments/real_world_graphs/karate.txt \
                       -edges experiments/real_world_graphs/karate.txt \
                       -u

echo ""
echo "Karate Undirected - No Noise"
executables/score_only -graph experiments/real_world_graphs/karate.txt \
                       -edges experiments/no_edges.txt \
                       -u

echo ""
echo "Foodweb Undirected - All Noise"
executables/score_only -graph experiments/real_world_graphs/maayan-foodweb.txt \
                       -edges experiments/real_world_graphs/maayan-foodweb.txt \
                       -u

echo ""
echo "Foodweb Undirected - No Noise"
executables/score_only -graph experiments/real_world_graphs/maayan-foodweb.txt \
                       -edges experiments/no_edges.txt \
                       -u

echo ""
echo "Foodweb Directed - All Noise"
executables/score_only -graph experiments/real_world_graphs/maayan-foodweb.txt \
                       -edges experiments/real_world_graphs/maayan-foodweb.txt \
                       -d

echo ""
echo "Foodweb Directed - No Noise"
executables/score_only -graph experiments/real_world_graphs/maayan-foodweb.txt \
                       -edges experiments/no_edges.txt \
                       -d

echo ""
echo "Cora Undirected - All Noise"
executables/score_only -graph experiments/real_world_graphs/cora.txt \
                       -edges experiments/real_world_graphs/cora.txt \
                       -u

echo ""
echo "Cora Undirected - No Noise"
executables/score_only -graph experiments/real_world_graphs/cora.txt \
                       -edges experiments/no_edges.txt \
                       -u

echo ""
echo "Cora Directed - All Noise"
executables/score_only -graph experiments/real_world_graphs/cora.txt \
                       -edges experiments/real_world_graphs/cora.txt \
                       -d

echo ""
echo "Cora Directed - No Noise"
executables/score_only -graph experiments/real_world_graphs/cora.txt \
                       -edges experiments/no_edges.txt \
                       -d

echo ""
echo "Eucore Undirected - All Noise"
executables/score_only -graph experiments/real_world_graphs/eucore.txt \
                       -edges experiments/real_world_graphs/eucore.txt \
                       -u

echo ""
echo "Eucore Undirected - No Noise"
executables/score_only -graph experiments/real_world_graphs/eucore.txt \
                       -edges experiments/no_edges.txt \
                       -u

echo ""
echo "Eucore Directed - All Noise"
executables/score_only -graph experiments/real_world_graphs/eucore.txt \
                       -edges experiments/real_world_graphs/eucore.txt \
                       -d

echo ""
echo "Eucore Directed - No Noise"
executables/score_only -graph experiments/real_world_graphs/eucore.txt \
                       -edges experiments/no_edges.txt \
                       -d

echo ""
echo "Pol Blogs Undirected - All Noise"
executables/score_only -graph experiments/real_world_graphs/pol_blogs.txt \
                       -edges experiments/real_world_graphs/pol_blogs.txt \
                       -u

echo ""
echo "Pol Blogs Undirected - No Noise"
executables/score_only -graph experiments/real_world_graphs/pol_blogs.txt \
                       -edges experiments/no_edges.txt \
                       -u

echo ""
echo "Pol Blogs Directed - All Noise"
executables/score_only -graph experiments/real_world_graphs/pol_blogs.txt \
                       -edges experiments/real_world_graphs/pol_blogs.txt \
                       -d

echo ""
echo "Pol Blogs Directed - No Noise"
executables/score_only -graph experiments/real_world_graphs/pol_blogs.txt \
                       -edges experiments/no_edges.txt \
                       -d

echo ""
echo "Season 4 Undirected - All Noise"
executables/score_only \
    -graph experiments/real_world_graphs/season_4_undirected_edges.txt \
    -edges experiments/real_world_graphs/season_4_undirected_edges.txt \
    -u

echo ""
echo "Season 4 Undirected - No Noise"
executables/score_only \
    -graph experiments/real_world_graphs/season_4_undirected_edges.txt \
    -edges experiments/no_edges.txt \
    -u

echo ""
echo "Season 4 Directed - All Noise"
executables/score_only \
    -graph experiments/real_world_graphs/season_4_directed_edges.txt \
    -edges experiments/real_world_graphs/season_4_directed_edges.txt \
    -d

echo ""
echo "Season 4 Directed - No Noise"
executables/score_only \
    -graph experiments/real_world_graphs/season_4_directed_edges.txt \
    -edges experiments/no_edges.txt \
    -d
