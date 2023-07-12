#!/bin/bash

# Delete ALL former results

rm results/run*

./half_helper.sh &>> results/runtimes.txt
