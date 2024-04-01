#!/bin/bash

# If using gcc instead of g++, add a link to -lstdc++
#   Removed -Werror due to a warning within the Nauty/Traces code itself
g++ -Wall -Wextra -o test_nt_code -std=c++11 -O4 test_nt_code.cpp nt_wrappers/nt_wrappers.a

g++ -Wall -Wextra -o minimal_nt_example -std=c++11 minimal_nt_example.cpp nt_wrappers/nt_wrappers.a
