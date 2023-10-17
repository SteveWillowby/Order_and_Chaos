# CrossSpot

This repository contains the code package for the ICDM'15 / TKDE'16 paper:

**[A General Suspiciousness Metric for Dense Blocks in Multimodal Data](http://www.meng-jiang.com/pubs/crossspot-icdm15/crossspot-icdm15-paper.pdf).** 

**[Spotting Suspicious Behaviors in Multimodal Data: A General Metric and Algorithms](http://www.meng-jiang.com/pubs/crossspot-tkde16/crossspot-tkde16-paper.pdf).** 

[Meng Jiang](http://meng-jiang.com/), Alex Beutel, Peng Cui, Bryan Hooi, Shiqiang Yang, Christos Faloutsos.

## Usage

**crossspot.py**
1. set [c_local]: count of injected, foreground block, 1000 as default
2. generate_data(): generate random data and inject the block. See global variables for data information:
- [MAX_NUM_SEED]: number of seeds in the algorithm
- [k_data]: number of modes
- [vec_n_local]: size vector of block
- [vec_n_global]: size vector of data
- [c_global]: capital C for count of the data
3. load_data(): load from file (data.csv) to
- data: entry list + value
- item2lineno: [k_data] maps, each map is {item:no. entry in [data] (lineno)}
4. CrossSpot Algorithm
- Output: screen output with best accuracy performance (maximum F1 score) with precision, recall, and F1 (average F1 score).

**crossspot-less-dense.py**
- change [c_local] from 3000 down to 400, from denser block to less dense block: generate data, run CrossSpot algorithm
- Output: in report.csv, best accuracy performance, avarage F1 score

## Citation
If you find this repository useful in your research, please cite our paper:

```bibtex
@inproceedings{jiang2015general,
  title={A general suspiciousness metric for dense blocks in multimodal data},
  author={Jiang, Meng and Beutel, Alex and Cui, Peng and Hooi, Bryan and Yang, Shiqiang and Faloutsos, Christos},
  booktitle={2015 IEEE International Conference on Data Mining},
  pages={781--786},
  year={2015},
  organization={IEEE}
}

@article{jiang2016spotting,
  title={Spotting suspicious behaviors in multimodal data: A general metric and algorithms},
  author={Jiang, Meng and Beutel, Alex and Cui, Peng and Hooi, Bryan and Yang, Shiqiang and Faloutsos, Christos},
  journal={IEEE Transactions on Knowledge and Data Engineering},
  volume={28},
  number={8},
  pages={2187--2200},
  year={2016},
  publisher={IEEE}
}
