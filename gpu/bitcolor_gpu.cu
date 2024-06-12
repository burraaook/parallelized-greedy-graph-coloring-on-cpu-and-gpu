#include <iostream>
#include <cuda_runtime.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

// atomicCAS
#include <device_atomic_functions.h>

#define TEST_NUM 20

// Define the maximum number of bits for the BitSet
#define BITS_PER_WORD 32
#define MAX_BITS 512
#define INVALID_MAX_SIZE 1000

#define IDLE 0
#define RUNNING 1

// cuda error checking
#define cudaCheckError() {                                          \
    cudaError_t e=cudaGetLastError();                                 \
    if(e!=cudaSuccess) {                                              \
        printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
        exit(0); \
    }                                                                 \
}

// BitSet structure definition
typedef struct {
    uint32_t bits[MAX_BITS / BITS_PER_WORD];
} BitSet;

// adj list graph structure, csr format
struct Graph
{
    int num_nodes;
    int num_edges;
    int* offsets;
    int* adj_list;
};


// Initialize the BitSet to zero
__device__ void bitset_init(BitSet* bs) {
    for (size_t i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        bs->bits[i] = 0;
    }
}

// Set a bit
__device__ void bitset_set(BitSet* bs, size_t index) {
    if (index >= MAX_BITS) return;
    bs->bits[index / BITS_PER_WORD] |= (1U << (index % BITS_PER_WORD));
}

// Clear a bit
__device__ void bitset_clear(BitSet* bs, size_t index) {
    if (index >= MAX_BITS) return;
    bs->bits[index / BITS_PER_WORD] &= ~(1U << (index % BITS_PER_WORD));
}

// Get a bit
__device__ int bitset_get(const BitSet* bs, size_t index) {
    if (index >= MAX_BITS) return 0;
    return (bs->bits[index / BITS_PER_WORD] & (1U << (index % BITS_PER_WORD))) != 0;
}

// Reset all bits
__device__ void bitset_reset(BitSet* bs) {
    for (size_t i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        bs->bits[i] = 0;
    }
}

// Bitwise OR
__device__ void bitset_or(BitSet* dest, const BitSet* src) {
    for (size_t i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        dest->bits[i] |= src->bits[i];
    }
}

// Bitwise NOT
__device__ void bitset_not(BitSet* dest, const BitSet* src) {
    for (size_t i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        dest->bits[i] = ~src->bits[i];
    }
}

// Increment BitSet
__device__ void bitset_increment(BitSet* bs) {
    int carry = 1;
    for (size_t i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        uint64_t temp = (uint64_t)bs->bits[i] + carry;
        bs->bits[i] = temp & 0xFFFFFFFF;
        carry = temp >> 32;
        if (carry == 0) break;
    }
}

__device__ void bitset_assign_one(BitSet* bs) {
    // Reset all bits to 0
    bitset_reset(bs);
    // Set the rightmost bit to 1
    bitset_set(bs, 0);
}

// Helper function to print the BitSet (for debugging)
__device__ void bitset_print(const BitSet* bs) {
    for (int i = MAX_BITS - 1; i >= 0; --i) {
        printf("%d", bitset_get(bs, i));
    }
    printf("\n");
}

// is 0
__device__ bool bitset_is_zero(const BitSet* bs) {
    for (int i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        if (bs->bits[i] != 0) {
            return false;
        }
    }
    return true;
}

// Bitwise AND
__device__ void bitset_and(BitSet* dest, const BitSet* src1, const BitSet* src2) {
    for (size_t i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        dest->bits[i] = src1->bits[i] & src2->bits[i];
    }
}


// kernel function which operates depending on warp
__global__ void kernel(const int* offsets, const int* adj_list, BitSet* result, int* queue, 
    bool* conflict_table, int* warp_state, int* num_vertices, int* num_warps, bool* terminate,
    int* warp_vertex_table)
{
    // get the thread id
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    // get the warp id
    int wid = tid / 32;
    // get the lane id
    int lane = tid % 32;
    // first warp is the master warp
    if (wid == 0 && lane == 0) {
        
        // cursor for the queue
        int cursor = 0;

        while (true) {
            // if queue is empty, wait for all warps to finish, and then terminate
            if (cursor >= *num_vertices)
            {

                for (int i = 1; i < *num_warps; ++i) {
                    while (warp_state[i] != IDLE)
                    {
                        __threadfence();
                    }
                }

                // set terminate flag to true
                *terminate = true;
                __threadfence();

                break;
            }

            // iterate through the warp table
            for (int i = 1; i < *num_warps; ++i) {
                // get warp id and status
                int state = warp_state[i];
                // if warp is IDLE, assign tasks to the warp
                if (state == IDLE) {

                    // assign vertex
                    warp_vertex_table[i] = queue[cursor];

                    // set dct false
                    conflict_table[queue[cursor]] = false;

                    // increment cursor
                    cursor++;

                    atomicExch(&warp_state[i], RUNNING);
                    __threadfence();
                }

                if (cursor >= *num_vertices) {
                    break;
                }
            }
        }
    }

    // first thread of other warps
    else if (lane == 0) {
        while(true) {
            // wait for the master warp to assign tasks
            while (warp_state[wid] != RUNNING && !(*terminate))
            {
                __threadfence();
            }
            if (*terminate) {
                return;
            }

            // get the vertex id
            int vertex = warp_vertex_table[wid];

            // color_state
            BitSet color_state;
            bitset_init(&color_state);
            BitSet color;
            bitset_init(&color);

            // local vertex buffer
            int local_vertex[INVALID_MAX_SIZE];
            int invalid_size = 0;

            // traverse the neighbors
            for (int i = offsets[vertex]; i < offsets[vertex + 1]; ++i) {
                int neighbor = adj_list[i];
                if (neighbor > vertex)
                    continue;

                // check if data is valid
                if (conflict_table[neighbor] == true)
                {
                    // get the color of the neighbor
                    BitSet neighbor_color = result[neighbor];

                    // bitwise OR
                    bitset_or(&color_state, &neighbor_color);
                }
                else
                {
                    local_vertex[invalid_size] = neighbor;
                    invalid_size++;
                }
            }
            while (true) {
                bool all_valid = true;
                for (int i = 0; i < invalid_size; ++i) {
                    if (i >= INVALID_MAX_SIZE)
                    {
                        printf("Invalid size exceeded\n");
                        break;
                    }
                    if (local_vertex[i] == -1)
                        continue;
                    if (conflict_table[local_vertex[i]] == true)
                    {
                        // get the color of the neighbor
                        BitSet neighbor_color = result[local_vertex[i]];

                        // bitwise OR
                        bitset_or(&color_state, &neighbor_color);
                        local_vertex[i] = -1;
                    }
                    else
                        all_valid = false;
                }
                if (all_valid)
                    break;
                __threadfence();
            }
            // get the color
            if (bitset_is_zero(&color_state)) {
                bitset_assign_one(&color);
            } else {
                BitSet temp;
                bitset_not(&temp, &color_state);
                bitset_increment(&color_state);
                bitset_and(&color, &temp, &color_state);
            }

            // assign the color to the vertex
            result[vertex] = color;

            // get yourself idle
            atomicExch(&warp_state[wid], IDLE);

            // set dct true
            conflict_table[vertex] = true;
            __threadfence();
            if (*terminate) {
                return;
            }        
        }
    }

    __syncwarp();
}



void free_graph(Graph& graph) {
    delete[] graph.offsets;
    delete[] graph.adj_list;
}

// sort based on the degree of the nodes
// smaller vertex id means higher degree
void sort_adj_list(std::unordered_map<int, std::unordered_set<int>>& adj_list, bool sortOption=true)
{
    std::vector<std::pair<int, int>> degree;
    for (const auto& p : adj_list) {
        degree.push_back({p.first, p.second.size()});
    }

    if (!sortOption)
        return;

    std::cerr << "\nSorting the graph" << std::endl;
    // sort based on the degree, smaller vertex id means higher degree
    // stable sort
    std::sort(degree.begin(), degree.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second > b.second || (a.second == b.second && a.first < b.first);
    });

    // mapping from old vertex id to new vertex id
    std::unordered_map<int, int> mapping;
    for (int i = 0; i < degree.size(); ++i) {
        mapping[degree[i].first] = i;
    }
    
    // remap the vertices
    std::unordered_map<int, std::unordered_set<int>> new_adj_list;
    for (const auto& p : adj_list) {
        int u = mapping[p.first];
        for (int v : p.second) {
            // map the vertex id
            int new_v = mapping[v];
            new_adj_list[u].insert(new_v);

        }
    }

    adj_list = new_adj_list;

    std::cerr << "Sorting done\n" << std::endl;
}

// example format
// 0 1
// 0 2
// 1 2
// 2 3
void read_graph(const std::string& filename, Graph& graph, bool sortOption=true)
{
    std::cerr << "\nReading graph from file: " << filename << std::endl;
    std::ifstream file(filename);
    std::string line;
    std::unordered_map<int, std::unordered_set<int>> adj_list;
    int num_edges = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) {
            break;
        }
        adj_list[u].insert(v);
        adj_list[v].insert(u);
        num_edges += 2;
    }
    file.close();

    // sort the graph
    sort_adj_list(adj_list, sortOption);
        // check if it is starting from 0, if not make it start from 0
    if (!sortOption && adj_list.find(0) == adj_list.end())
    {
        std::cerr << "Vertex 0 not found, remapping the vertices" << std::endl;
        std::unordered_map<int, int> mapping;
        int new_vertex_id = 0;
        for (const auto& p : adj_list) {
            mapping[p.first] = new_vertex_id;
            new_vertex_id++;
        }

        // remap the vertices
        std::unordered_map<int, std::unordered_set<int>> new_adj_list;
        for (const auto& p : adj_list) {
            int u = mapping[p.first];
            for (int v : p.second) {
                // map the vertex id
                int new_v = mapping[v];
                new_adj_list[u].insert(new_v);

            }
        }

        adj_list = new_adj_list;
    }

    graph.num_nodes = adj_list.size();
    graph.num_edges = num_edges;
    graph.offsets = new int[graph.num_nodes + 1];
    graph.adj_list = new int[num_edges];
    int offset = 0;
    graph.offsets[0] = 0;
    for (int i = 0; i < graph.num_nodes; ++i) {
        graph.offsets[i + 1] = graph.offsets[i] + adj_list[i].size();
        for (int v : adj_list[i]) {
            graph.adj_list[offset++] = v;
        }
    }

    std::cerr << "Reading done\n" << std::endl;
}

void print_graph(const Graph& graph) {

    std::cerr << "Number of nodes: " << graph.num_nodes << std::endl;
    std::cerr << "Number of edges: " << graph.num_edges << std::endl;

    // vertex and its neighbors
    for (int i = 0; i < graph.num_nodes; ++i) {
        std::cerr << i << ": ";
        for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; ++j) {
            std::cerr << graph.adj_list[j] << " ";
        }
        std::cerr << std::endl;
    }
}

int gpu_info()
{
    // print total number of cores
    int num_cores = 0;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    // maximum number of streaming multiprocessors
    num_cores = prop.multiProcessorCount;
    std::cerr << "Number of cores: " << num_cores << std::endl;

    // maximum number of threads per block
    int max_threads = prop.maxThreadsPerBlock;
    std::cerr << "Maximum number of threads per block: " << max_threads << std::endl;

    // number of warps, threads per block / warp size * number of cores
    int num_warps = max_threads / 32 * num_cores;
    std::cerr << "Number of warps: " << num_warps << std::endl;

    return num_warps;
}

bool compare_bitset(const BitSet& a, const BitSet& b) 
{
    for (int i = 0; i < MAX_BITS / BITS_PER_WORD; ++i) {
        if (a.bits[i] != b.bits[i]) {
            return false;
        }
    }
    return true;
}

// bitset get host
int bitset_get_host(const BitSet& bs, int index) 
{
    if (index >= MAX_BITS) return 0;
    return (bs.bits[index / BITS_PER_WORD] & (1U << (index % BITS_PER_WORD))) != 0;
}

std::string bitset_get_str(const BitSet& bs) 
{
    std::string result;
    // iterate through the bits
    for (int i = MAX_BITS - 1; i >= 0; --i) {
        result += std::to_string(bitset_get_host(bs, i));
    }

    return result;

}

int check_conflict(const Graph& graph, const BitSet* result) 
{
    int num_conflicts = 0;
    for (int i = 0; i < graph.num_nodes; ++i) {
        for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; ++j) {
            if (compare_bitset(result[i], result[graph.adj_list[j]])) {
                // std::cerr << "Conflict: " << i << " " << graph.adj_list[j] << std::endl;
                num_conflicts++;
            }
        }
    }
    std::cerr << "Number of conflicts: " << num_conflicts << std::endl;
    return num_conflicts;
}

void count_colors(const Graph& graph, const BitSet* result) 
{
    std::unordered_set<std::string> colors;
    for (int i = 0; i < graph.num_nodes; ++i) {
        colors.insert(bitset_get_str(result[i]));
    }
    std::cerr << "Number of colors: " << colors.size() << std::endl;
}
void gpu_process(const Graph& graph, int num_blocks = 30, int block_size = 1024, std::string dataset_name = "dataset")
{
    // allocate the result on the device
    // vertex, color
    BitSet* result;

    // store graph on the device
    int* offsets;
    int* adj_list;

    // store queue on the device
    // it is a queue of vertices, from 0 to num_nodes - 1
    int* queue = new int[graph.num_nodes];
    for (int i = 0; i < graph.num_nodes; ++i) {
        queue[i] = i;
    }

    // copy the result back to the host
    BitSet* result_host = new BitSet[graph.num_nodes];

    // measure time
    cudaEvent_t start, stop;
    // launch the kernel
    std::cerr << "Launching kernel" << std::endl;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    cudaMalloc(&result, graph.num_nodes * sizeof(BitSet));
    cudaCheckError();


    cudaMalloc(&offsets, (graph.num_nodes + 1) * sizeof(int));
    cudaMalloc(&adj_list, graph.num_edges * sizeof(int));
    cudaMemcpy(offsets, graph.offsets, (graph.num_nodes + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(adj_list, graph.adj_list, graph.num_edges * sizeof(int), cudaMemcpyHostToDevice);
    cudaCheckError();

    int* d_queue;
    cudaMalloc(&d_queue, graph.num_nodes * sizeof(int));
    cudaMemcpy(d_queue, queue, graph.num_nodes * sizeof(int), cudaMemcpyHostToDevice);
    cudaCheckError();

    // store number of vertices on the device
    int* num_vertices;
    cudaMalloc(&num_vertices, sizeof(int));
    cudaMemcpy(num_vertices, &graph.num_nodes, sizeof(int), cudaMemcpyHostToDevice);
    cudaCheckError();

    // create data conflict table
    // <vertex, valid> map
    bool* conflict_table;
    // initialize to true
    cudaMalloc(&conflict_table, graph.num_nodes * sizeof(bool));
    cudaMemset(conflict_table, 1, graph.num_nodes * sizeof(bool));
    cudaCheckError();

    // get the number of warps
    // int num_warps = gpu_info();
    // gpu_info();
    // int num_blocks = 30;
    // int block_size = 1024;
    int num_warps = num_blocks * block_size / 32;

    // create warp state table
    int* warp_state;
    cudaMalloc(&warp_state, num_warps * sizeof(int));
    // initialize to IDLE
    cudaMemset(warp_state, IDLE, num_warps * sizeof(int));
    cudaCheckError();

    // store number of warps on the device
    int* num_warps_device;
    cudaMalloc(&num_warps_device, sizeof(int));
    cudaMemcpy(num_warps_device, &num_warps, sizeof(int), cudaMemcpyHostToDevice);
    cudaCheckError();

    // store terminate flag on the device
    bool* terminate;
    cudaMalloc(&terminate, sizeof(bool));
    // initialize to false
    cudaMemset(terminate, 0, sizeof(bool));
    cudaCheckError();

    // create table for <warp, vertex> mapping
    int* warp_vertex_table;
    cudaMalloc(&warp_vertex_table, num_warps * sizeof(int));
    // initialize to -1
    cudaMemset(warp_vertex_table, -1, num_warps * sizeof(int));
    cudaCheckError();

    kernel<<<num_blocks, block_size>>>(offsets, adj_list, result, d_queue, conflict_table, warp_state, num_vertices, num_warps_device, terminate, warp_vertex_table);
    cudaCheckError();
    cudaDeviceSynchronize();
    cudaMemcpy(result_host, result, graph.num_nodes * sizeof(BitSet), cudaMemcpyDeviceToHost);
    cudaCheckError();
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cerr << "Time: " << milliseconds << " ms" << std::endl;


    // check if there is a conflict
    std::cerr << "\nChecking conflicts" << std::endl;
    int conflict = check_conflict(graph, result_host);

    // count the number of colors
    std::cerr << "\nCounting colors" << std::endl;
    count_colors(graph, result_host);

    // free the memory
    delete[] result_host;
    delete[] queue;
    cudaFree(result);
    cudaFree(offsets);
    cudaFree(adj_list);
    cudaFree(d_queue);
    cudaFree(conflict_table);
    cudaFree(warp_state);
    cudaFree(num_vertices);
    cudaFree(num_warps_device);
    cudaFree(terminate);
    cudaFree(warp_vertex_table);


    // write result to file as: dataset_name num_blocks block_size time conflict
    std::ofstream file("results.txt", std::ios_base::app);
    file << dataset_name << " " << num_blocks << " " << block_size << " " << milliseconds << " " << conflict << std::endl;
    file.close();



}
void run_test() {
    std::vector<std::string> filenames = {"../datasets/EF.txt", "../datasets/CD.txt", "../datasets/RC.txt", 
                                         "../datasets/CA.txt", "../datasets/RP.txt", "../datasets/RT.txt",
                                         "../datasets/CL.txt"};

    std::vector<int> num_blocks = {1, 8, 16, 30, 1, 1, 1};
    std::vector<int> block_size = {1024, 1024, 1024, 1024, 512, 256, 128};

    for (int i = 0; i < filenames.size(); ++i) {
        for (int k = 0; k < num_blocks.size(); ++k) {
            for (int j = 0; j < TEST_NUM; ++j) {
                std::cerr << "\n\nDataset: " << filenames[i] << " Test: " << j << " Blocks: " << num_blocks[k] << " Block size: " << block_size[k] << std::endl;
                Graph graph;
                read_graph(filenames[i], graph);
                gpu_process(graph, num_blocks[k], block_size[k], filenames[i]);
                free_graph(graph);
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc == 2 && std::string(argv[1]) == "run_tests") {
        run_test();
        return 0;
    }

    if (argc == 2 && std::string(argv[1]) == "gpu_info") {
        gpu_info();
        return 0;
    }

    if (argc == 2 && std::string(argv[1]) == "-h") {
        std::cerr << "Format: " << argv[0] << " <dataset_path> <number_of_blocks> <block_size> <sort_option>" << std::endl;
        std::cerr << "sort_option: sort_yes or sort_no" << std::endl;
        return 0;
    }

    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <dataset_path> <number_of_blocks> <block_size> <sort_option>" << std::endl;
        std::cerr << "Usage 2: " << argv[0] << " run_tests" << std::endl;
        std::cerr << "Usage 3: " << argv[0] << " gpu_info" << std::endl;
        std::cerr << "Usage 4: " << argv[0] << " -h" << std::endl;
        return 1;
    }

    std::string dataset_path = argv[1];
    int num_blocks = std::stoi(argv[2]);
    int block_size = std::stoi(argv[3]);
    
    // sort_option -> sort_yes or sort_no arguments
    bool sort_option = std::string(argv[4]) == "sort_yes";

    if (block_size % 32 != 0 || block_size > 1024 || block_size < 64) {
        std::cerr << "Error: block_size must be a multiple of 32, between 64 and 1024" << std::endl;
        return 1;
    }

    Graph graph;
    read_graph(dataset_path, graph, sort_option);
    gpu_process(graph, num_blocks, block_size, dataset_path);
    free_graph(graph);

    return 0;
}