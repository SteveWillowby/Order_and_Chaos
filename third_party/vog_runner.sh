#!/bin/bash

rm karate_deets.txt
python3 alternations.py karate >> karate_deets.txt

rm season_4_deets.txt
python3 alternations.py season_4 >> season_4_deets.txt

rm foodweb_deets.txt
python3 alternations.py foodweb >> foodweb_deets.txt
