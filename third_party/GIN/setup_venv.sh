#!/bin/bash

python3 -m venv cpu_virtual_env
source cpu_virtual_env/bin/activate
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install torch_geometric
pip install scikit-learn
pip install pyg-lib -f https://data.pyg.org/whl/torch-2.5.0+cpu.html
deactivate
