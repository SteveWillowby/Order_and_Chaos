#!/bin/bash

# Delete ALL former results

rm results/*

./helper.sh &>> results/runtimes.txt
