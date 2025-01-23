#!/bin/bash

source cpu_virtual_env/bin/activate
python link_pred.py $1 $2
deactivate
