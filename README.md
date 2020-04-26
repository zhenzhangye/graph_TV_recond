# graph_TV_recond
This code implements the following [paper](https://arxiv.org/abs/2002.12236):

> **Optimization of Graph Total Variation via Active-Set-based Combinatorial Reconditioning**
> *Ye, Z., Möllenhoff, T., Wu, T., Cremers, D.; In Proceedings of the International Conference on Artificial Intelligence and Statistics (AISTATS), 2020*

We show that nested-forest decomposition of the inactive edges yields a guaranteed local linear convergence rate. Further, we propose a practical greedy heuristic which realizes such nested decompositions and show in several numerical experiments that our reconditioning strategy, when applied to proximal gradient or primal-dual hybrid gradient algorithm, achieves competitive performances.

## 1. Requirements

This code has following dependencies:

0) MATLAB (code was tested on R2019a)

1) [Boost C++ Library](https://www.boost.org/)

## 2. Installation

- Clone this repository and run `setup_mex` in MATLAB.
- If you want to run the code on graph cut problems, download the datasets in the following paper:
```
@inproceedings{goldberg2011maximum,
  title={Maximum flows by incremental breadth-first search},
  author={Goldberg, Andrew V and Hed, Sagi and Kaplan, Haim and Tarjan, Robert E and Werneck, Renato F},
  booktitle={European Symposium on Algorithms},
  pages={457--468},
  year={2011},
  organization={Springer}
}
```
and change the `files` variable in `E_graphcut/main.m` to the corresponding directories.

## 3. Reproduce the results in paper

To reproduce the results in paper, run the main file in each folder to get the figures respectively.

## 4. Publication

If you make use of the library in any form in a scientific publication, please refer to `https://github.com/zhenzhangye/graph_TV_recond` and cite the paper

```
@article{ye2020optimization,
  title={Optimization of Graph Total Variation via Active-Set-based Combinatorial Reconditioning},
  author={Ye, Zhenzhang and M{\"o}llenhoff, Thomas and Wu, Tao and Cremers, Daniel},
  journal={arXiv preprint arXiv:2002.12236},
  year={2020}
}
```
