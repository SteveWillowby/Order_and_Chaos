#!/bin/bash

python3 score_decompositions.py Ktruss > results/Ktruss_score_file.txt
python3 score_decompositions.py Ktruss core_only > results/core_only/Ktruss_score_file.txt

python3 score_decompositions.py SUBDUE > results/Subdue_score_file.txt
python3 score_decompositions.py SUBDUE core_only > results/core_only/Subdue_score_file.txt

python3 score_decompositions.py VoG > results/VoG_score_file.txt
python3 score_decompositions.py VoG core_only > results/core_only/VoG_score_file.txt
