
#ifndef GRAPH_HPP
#define GRAPH_HPP
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <set>
#include <bitset>

#include <thread>
#include <queue>
#include <functional>

#include "dct.hpp"

#define BITSET_SIZE 512

// graph with adjacency list
class Graph 
{
protected:

    // adjacency list
    std::unordered_map<int, std::list<int>> adj;
    std::vector<int> vertices;
    int num_color = 0;
public:
    Graph(std::string filename) { read_graph(filename); }
    void addEdge(int v, int w) { adj[v].push_back(w); adj[w].push_back(v); }
    int read_graph(std::string filename);
    int getNumVertices() { return adj.size(); }
    std::list<int> getNeighbors(int v) { return adj[v]; }
    std::unordered_map<int, std::list<int>> getAdj() { return adj; }
    bool checkConflict(std::unordered_map<int, std::bitset<BITSET_SIZE>> result);

    std::unordered_map<int, int> basicGreedyColoring();
    int countColorsBitset(std::unordered_map<int, std::bitset<BITSET_SIZE>> bitset_map);
    std::unordered_map<int, std::bitset<BITSET_SIZE>> bitWiseColoring();
    std::bitset<BITSET_SIZE> increment_bitset(std::bitset<BITSET_SIZE> bitset);
    int getNumColor() { return num_color; }

    // count colors
    int countColors(std::unordered_map<int, int> color_array);

    ~Graph() { adj.clear(); }
};

// thread status enum
enum ThreadStatus
{
    IDLE,
    RUNNING,
    WAITING,
    TERMINATED
};

// inherit threaded graph coloring
class ParallelGraph : public Graph
{
private:
    // result
    std::unordered_map<int, std::bitset<BITSET_SIZE>> result;

    // number of threads
    int num_threads;

    // threads
    std::vector<std::thread> threads;

    // FIFO, tasks
    // std::vector<int> queue;

    // data conflict table
    DataConflictTable dct;

    // thread table: <thread_id, (v_id, run)>, run: thread status
    std::map<int, std::pair<int, ThreadStatus>> thread_table;

    // vertex thread table: <vertex_id, thread_id>
    std::unordered_map<int, int> vertex_thread_table;

    // flag for termination
    bool terminate = false;

    // sort option
    bool sortOption = false;

    // initialize tables and queue
    // void initialize_tables();

    // master thread
    void master_thread();

    // each processing thread
    void processing_thread(int thread_id);

    // worker function: thread_id, vertex_id, neighbors
    void worker(int thread_id, int vertex, std::list<int> *neighbors);

public:

    ParallelGraph(std::string filename, int num_threads, bool sortOption = false);
    // overloaded bit wise coloring
    std::unordered_map<int, std::bitset<BITSET_SIZE>> bitWiseColoring();
    // overloaded sort vertices by degree
    void sortVerticesByDegree();

    bool checkConflict();
};

#endif // GRAPH_HPP