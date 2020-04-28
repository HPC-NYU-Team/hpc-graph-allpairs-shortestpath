# Parallel All Pairs Shortest Path in Large Graphs
## High Performance Computing


### Setup:
- ```pip install networkx```
- ```conda env create -f environment.yml```
- ```conda activate hpc-project```
- ```python createGraph.py <order> <degree>```


### Running adj-matrix 
- ```srun --nodes=10 --tasks-per-node=1 --cpus-per-task=8 --pty /bin/bash```
- ```mpic++ adj-matrix-parallel-imp.cpp -fopenmp -o mpi ```
- ``` mpirun --mca btl_tcp_if_include eth0 --mca btl '^openib' --np 10 mpi 1000 5 100 ``` 
- ``` mpirun --mca btl_tcp_if_include eth0 --mca btl '^openib' --np 10 mpi <order> <degree> <chunk-size> ```
