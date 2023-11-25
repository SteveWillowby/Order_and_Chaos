#!/bin/bash

rm karate_deets.txt
nice -2 python3 alternations.py karate >> karate_deets.txt
rm greedySelection*

rm season_4_deets.txt
nice -2 python3 alternations.py season_4 >> season_4_deets.txt
rm greedySelection*

rm foodweb_deets.txt
nice -2 python3 alternations.py foodweb >> foodweb_deets.txt
rm greedySelection*
