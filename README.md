# accelerating-greedy-graph-coloring
Acceleration of Greedy Graph Coloring with the Application of Bitcolor on applying it to CPU and GPU.

## compilation
```
make
```

## run cpu
```
./bitcolor_cpu <algorithm_name> <dataset_path>
```

- Algorithms: basic_greedy, bitwise_greedy, p_bitcolor

- p_bitcolor options:   
[NUM_THREADS] [sort_option]  
sort_option: sort_yes, sort_no  



## run gpu
```
./bitcolor_gpu <dataset_path> <number_of_blocks> <block_size> <sort_option>
```
sort_option: sort_yes or sort_no
