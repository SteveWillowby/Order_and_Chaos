#!/bin/bash

rm results/karate_subdue.txt
nice -2 python3 alternations.py karate subdue >> results/karate_subdue.txt

rm results/season_4_subdue.txt
nice -2 python3 alternations.py season_4 subdue >> results/season_4_subdue.txt

rm results/foodweb_subdue.txt
nice -2 python3 alternations.py foodweb subdue >> results/foodweb_subdue.txt

rm results/pol_blogs_subdue.txt
nice -2 python3 alternations.py pol_blogs subdue >> results/pol_blogs_subdue.txt

rm results/eucore_subdue.txt
nice -2 python3 alternations.py eucore subdue >> results/eucore_subdue.txt

rm results/cora_subdue.txt
nice -2 python3 alternations.py cora subdue >> results/cora_subdue.txt
