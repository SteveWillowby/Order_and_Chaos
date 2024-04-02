#!/bin/bash

echo "######################################################"
echo "###         Compiling the SCHENO Utilities         ###"
echo "######################################################"
echo ""

# Create the .o files
g++ -c edge_sampler.cpp genetic_alg_search.cpp int_edge_sampler.cpp noise_prob_choice.cpp scoring_function.cpp thread_pool_scorer.cpp thread_pool_wl_sim.cpp wl_measures.cpp Jonker_Volgenant/src/assign2DCBasic.c

# Create a library
ar rvs scheno_only.a edge_sampler.o genetic_alg_search.o int_edge_sampler.o noise_prob_choice.o scoring_function.o thread_pool_scorer.o thread_pool_wl_sim.o wl_measures.o assign2DCBasic.o

# Combine the nauty & traces library with this library
ar -M <scheno.mri
