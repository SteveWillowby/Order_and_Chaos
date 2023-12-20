#!/bin/bash

rm results/karate_kcore.txt
nice -2 python3 alternations.py karate kcore >> results/karate_kcore.txt

rm results/season_4_kcore.txt
nice -2 python3 alternations.py season_4 kcore >> results/season_4_kcore.txt

rm results/foodweb_kcore.txt
nice -2 python3 alternations.py foodweb kcore >> results/foodweb_kcore.txt

rm results/pol_blogs_kcore.txt
nice -2 python3 alternations.py pol_blogs kcore >> results/pol_blogs_kcore.txt

rm results/eucore_kcore.txt
nice -2 python3 alternations.py eucore kcore >> results/eucore_kcore.txt

rm results/cora_kcore.txt
nice -2 python3 alternations.py cora kcore >> results/cora_kcore.txt
