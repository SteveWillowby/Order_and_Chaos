#!/bin/bash

# If using gcc instead of g++, add a link to -lstdc++
#   Removed -Werror due to a warning within the Nauty/Traces code itself
g++ -Wall -Wextra -o test_nt_code -std=c++11 -O4 test_nt_code.cpp scheno/scheno.a

g++ -Wall -Wextra -o minimal_nt_example -std=c++11 minimal_nt_example.cpp scheno/scheno.a

g++ -Wall -Wextra -o executables/graph_enumeration -std=c++11 -O4 graph_enumeration.cpp scheno/scheno.a

g++ -Wall -Wextra -o executables/mthread_scorer -std=c++11 -O4 test_threaded_scorer.cpp scheno/scheno.a -lpthread

g++ -Wall -Wextra -o executables/SCHENO_score -std=c++11 -O4 SCHENO_score.cpp scheno/scheno.a -lpthread

g++ -Wall -Wextra -o executables/SCHENO_ga -std=c++11 -O4 SCHENO_ga.cpp scheno/scheno.a -lpthread
