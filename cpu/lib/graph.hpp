
#ifndef GRAPH_HPP
#define GRAPH_HPP
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <set>
#include <bitset>

#include <thread>
#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>
// lock
#include <shared_mutex>

#include "dct.hpp"

// graph with adjacency list
class Graph 
{
protected:

    // adjacency list
    std::map<int, std::list<int>> adj;
    int num_color = 0;
public:
    Graph(std::string filename) { read_graph(filename); }
    void addEdge(int v, int w) { adj[v].push_back(w); adj[w].push_back(v); }
    int read_graph(std::string filename);
    int getNumVertices() { return adj.size(); }
    std::list<int> getNeighbors(int v) { return adj[v]; }
    std::map<int, std::list<int>> getAdj() { return adj; }
    bool checkConflict(std::map<int, std::bitset<512>> result);

    std::map<int, int> basicGreedyColoring();
    int countColorsBitset(std::map<int, std::bitset<512>> bitset_map);
    std::map<int, std::bitset<512>> bitWiseColoring();
    std::bitset<512> increment_bitset(std::bitset<512> bitset);
    int getNumColor() { return num_color; }
    // void sortVerticesByDegree();

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
    std::map<int, std::bitset<512>> result;

    // number of threads
    int num_threads;

    // threads
    std::vector<std::thread> threads;

    // FIFO, tasks
    std::vector<int> queue;

    // data conflict table
    DataConflictTable dct;

    // thread table: <thread_id, (v_id, run)>, run: thread status
    std::map<int, std::pair<int, ThreadStatus>> thread_table;

    // vertex thread table: <vertex_id, thread_id>
    std::map<int, int> vertex_thread_table;


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

    // mutex
    std::mutex mtx;

public:

    ParallelGraph(std::string filename, int num_threads, bool sortOption = false);
    // overloaded bit wise coloring
    std::map<int, std::bitset<512>> bitWiseColoring();
    // overloaded sort vertices by degree
    void sortVerticesByDegree();

    bool checkConflict();
};

#endif // GRAPH_HPP