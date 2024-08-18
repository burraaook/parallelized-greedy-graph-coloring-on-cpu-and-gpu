# PBitCo: Parallelization of BitColor Algorithm via Multithreading and GPU
Acceleration of Greedy Graph Coloring with the Application of Bitcolor on applying it to CPU and GPU.

## Dependencies
- CUDA (Tested on 12.4)
- g++ (Tested on 11.4.0)

## Parameters
### Algorithms
- Basic Greedy: Regular greedy coloring algorithm that does nothing special.
- Bitwise Greedy: Uses bit operations instead of color traversal at each node.
- PBitCo: Parallelized version of Bitwise Greedy.

### Number of Processing Units
- CPU: Number of threads
- GPU: Number of blocks, and size of each block.

### Sorted Graph
- In the preprocessing step graph can be sorted which results in higher degree vertices will have lower vertex id. Algorithm can work on both unsorted and sorted graphs.

## Usage
### Compilation
```
make
```

### Run CPU
```
./bitcolor_cpu <algorithm_name> <dataset_path>
```

- Algorithms: basic_greedy, bitwise_greedy, p_bitcolor

- options:   
[NUM_THREADS] [sort_option]  
sort_option: sort_yes, sort_no  
#### Examples
```
./bitcolor_cpu p_bitcolor graph.txt 8 sort_no
```
```
./bitcolor_cpu bitwise_greedy graph.txt sort_yes
```
### Run GPU
```
./bitcolor_gpu <dataset_path> <number_of_blocks> <block_size> <sort_option>
```
sort_option: sort_yes or sort_no
#### Example
```
./bitcolor_gpu graph.txt 1 1024 sort_no
```
## Input File Format
```
<vertex_id> <neighbor_id>
```

### Example
```
0   1
0   4
0   2
1   2
1   4
4   5
1   3
2   6
3   5
3   7
3   6
1   6
5   7
6   7
0   8
2   8
4   8
```
