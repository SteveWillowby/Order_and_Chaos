#!/bin/bash

rm results/karate_vog.txt
nice -2 python3 alternations.py karate >> results/karate_vog.txt
rm greedySelection*

rm results/season_4_vog.txt
nice -2 python3 alternations.py season_4 >> results/season_4_vog.txt
rm greedySelection*

rm results/foodweb_vog.txt
nice -2 python3 alternations.py foodweb >> results/foodweb_vog.txt
rm greedySelection*
