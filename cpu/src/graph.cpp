#include "../lib/graph.hpp"


int Graph::read_graph(std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(file, line)) {
        int v, w;
        std::istringstream iss(line);

        // if line is %, #, or empty, skip
        if (line[0] == '%' || line[0] == '#' || line.empty()) {
            continue;
        }
        if (!(iss >> v >> w)) {
            std::cerr << "Error: could not read line " << line << std::endl;
            return 1;
        }
        this->addEdge(v, w);
    }

    file.close();
    return 0;
}

std::map<int, int> Graph::basicGreedyColoring() {

    // initalize color_array, and color_flag
    std::map<int, int> color_array; // <vertex, color>
    std::map<int, int> color_flag; // <color, flag> flag = 1 if color is used, flag = 0 if color is not used

    // iterate through all vertices
    for (auto const& x : adj) {
        // stage 0: neighbor vertices traversal

        // get the vertex
        int vertex = x.first;

        // get the neighbors of the vertex
        std::list<int> *neighbors = &adj[vertex];

        // iterate through the neighbors
        for (auto const& y : *neighbors) {

            // get the neighbor
            int neighbor = y;

            // if the neighbor has been colored, set the color_flag to 1
            if (color_array.find(neighbor) != color_array.end()) {
                color_flag[color_array[neighbor]] = 1;
            }
        }

        // stage 1: color traversal
        int color = 0;
        while (color_flag[color] == 1) {
            color++;
        }

        // set the color_flag array back to 0
        for (auto const& y : *neighbors) {
            if (color_array.find(y) != color_array.end()) {
                color_flag[color_array[y]] = 0;
            }
        }

        // stage 2: color assignment
        color_array[vertex] = color;
    }

    // set the number of colors used
    num_color = color_flag.size();

    // return the color_array
    return color_array;
}

std::map<int, std::bitset<512>> Graph::bitWiseColoring() {
    std::map<int, std::bitset<512>> color_array; // <vertex, color>

    // color_state is used for checking if a color is used, represented with a bitset
    std::bitset<512> color_state;
    color_state.reset();

    // iterate through all vertices
    for (auto const& x : adj) {
        // stage 0: neighbor vertices traversal

        // get the vertex
        int vertex = x.first;

        // get the neighbors of the vertex
        std::list<int> *neighbors = &adj[vertex];

        // iterate through the neighbors
        for (auto const& y : *neighbors) {

            // get the neighbor
            // int neighbor = y;

            // get color of neighbor
            // std::bitset<512> neighbor_color = color_array[neighbor];
            // color_state |= neighbor_color;
            color_state |= color_array[y];
        }
        // stage 1: color traversal
        // if color state = 0, then color = 1
        // else, color = (~color_state) & (color_state + 1'b1)
        std::bitset<512> color;
        // std::cout << "color_state: " << color_state << "\n";
        if (color_state == 0) {
            color = 1;
        } 
        else {
            color = (~color_state) & (increment_bitset(color_state));
        }

        // stage 2: color assignment
        color_array[vertex] = color;

        // reset color_state
        color_state.reset();
    }
    return color_array;
}

std::bitset<512> Graph::increment_bitset(std::bitset<512> bitset) {
    
    for (std::size_t i = 0; i < bitset.size(); i++) {
        if (bitset[i] == 0) {
            bitset[i] = 1;
            break;
        } else {
            bitset[i] = 0;
        }
    }
    return bitset;
}

int Graph::countColorsBitset(std::map<int, std::bitset<512>> bitset_map) {
    std::set<std::string> colors;
    for (auto const& x : bitset_map) {
        colors.insert(x.second.to_string());
    }

    this->num_color = colors.size();
    return this->num_color;
}

// parallel graph coloring 
ParallelGraph::ParallelGraph(std::string filename, int num_threads) : Graph(filename), num_threads(num_threads) {
    // sort vertices by degree
    sortVerticesByDegree();

    // initialize threads
    threads.resize(num_threads);

    // fill result with empty bitset
    for (auto const& x : adj) {
        result[x.first].reset();
        dct.insert(x.first);
        result[x.first] = 0;
    }

    // initialize dct
    for (auto const& x : adj) {
        dct.insert(x.first);

        // set all vertices to valid, and uncolored
        dct.setValid(x.first);
        // dct.setUncolored(x.first);
    }

    // std::cerr << "ParallelGraph initialized" << std::endl;
}

// sort vertices by degree  
void ParallelGraph::sortVerticesByDegree() {
    // std::cerr << "Sorting vertices by degree" << std::endl;

    // Step 1: Sort vertices by degree
    for (auto const& x : adj) {
        queue.push_back(x.first);
    }
    // return;

    std::sort(queue.begin(), queue.end(), [this](int a, int b) {
        return adj[a].size() > adj[b].size();
    });
    // std::cerr << "Vertices sorted by degree" << std::endl;

    // Step 2: Precompute vertex positions in the sorted queue
    std::map<int, int> vertex_positions;
    for (size_t i = 0; i < queue.size(); ++i) {
        vertex_positions[queue[i]] = i;
    }

    // Step 3: Update adjacency list indexes
    std::map<int, std::list<int>> updated_adj;
    for (auto const& vertex : queue) {
        // Step 4: Update neighbor lists using precomputed positions
        std::list<int> updated_neighbors;
        for (auto const& neighbor : adj[vertex]) {
            auto neighbor_it = vertex_positions.find(neighbor);
            if (neighbor_it != vertex_positions.end()) {
                updated_neighbors.push_back(neighbor_it->second);
            }
        }
        updated_adj[vertex_positions[vertex]] = updated_neighbors;
    }
    // update queue indexes
    for (size_t i = 0; i < queue.size(); ++i) {
        queue[i] = vertex_positions[queue[i]];
    }
    
    // std::cerr << "Adjacency list updated" << std::endl;

    // Step 5: Reverse the queue, because algorithm uses pop_back
    std::reverse(queue.begin(), queue.end());

    // Step 6: Update the adjacency list
    adj = updated_adj;
}

// parallel bit wise coloring
std::map<int, std::bitset<512>> ParallelGraph::bitWiseColoring() {

    // start master thread
    master_thread();

    // join threads
    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
    }

    return result;
}

// master thread
void ParallelGraph::master_thread() {

    // start threads
    for (int i = 0; i < num_threads; i++) {
        threads[i] = std::thread(&ParallelGraph::processing_thread, this, i);

        // assign to thread table
        thread_table[i] = {queue[i], IDLE};
    }

    // assign tasks to threads till all vertices are colored
    while (true) {
        // check if all vertices are colored. check fifo queue
        if (queue.empty()) {
            // set terminate flag
            
            // wait for all threads to become idle
            for (int i = 0; i < num_threads; i++) {
                while (thread_table[i].second != IDLE);
            }

            // set terminate flag
            terminate = true;
            // break the loop
            break;
        }

        // iterate through the thread table
        for (auto const& x : thread_table) {
            int thread_id = x.first;
            ThreadStatus status = x.second.second;
            int vertex = -1;            
            // if the thread is not running, assign a task
            if (status == IDLE) {

                if (queue.empty()) {
                    break;
                }
                vertex = queue.back();

                // set dct to invalid
                dct.setInvalid(vertex);
                
                // add to vertex thread table
                // std::cerr << "coloring vertex= " << vertex << std::endl; 
                vertex_thread_table[vertex] = thread_id;

                // update thread table, get largest vertex from queue
                thread_table[thread_id] = {vertex, RUNNING};

                // remove the vertex from the queue
                queue.pop_back();
            }
        }        
    }
}

// processing thread
void ParallelGraph::processing_thread(int thread_id) {

    int vertex_id = thread_table[thread_id].first;

    // run until terminate flag is set
    while (true)
    {
        while (thread_table[thread_id].second == IDLE && !terminate);
        if (terminate) {
            break;
        }
        vertex_id = thread_table[thread_id].first;

        std::list<int> *neighbors = &adj[vertex_id];

        // worker function
        worker(thread_id, vertex_id, neighbors);
        // set thread status to idle
        thread_table[thread_id].second = IDLE;

        if (terminate) {
            break;
        }
    }

    // std::cerr << "Thread " << thread_id << " terminated" << std::endl;
}

int bit2num (std::bitset<512> bitset) {
    int num = 0;
    // start from left first 1 is the number for that bitset
    for (int i = 512; i >= 0; i--) {
        if (bitset[i] == 1) {
            num = i;
            return num;
        }
    }

    return num;
}


// worker function, address of neighbors
void ParallelGraph::worker(int thread_id, int vertex, std::list<int> *neighbors) { 

    // color_state is used for checking if a color is used, represented with a bitset
    std::bitset<512> color_state = 0;
    color_state.reset();

    // invalid vertices linked list
    std::list<int> invalid_vertices;

    // iterate through the neighbors, prune valid and uncolored neighbors, wait invalid neighbors
    for (auto const& neighbor : *neighbors) {

        // prune
        if (neighbor > vertex)
            continue;
            
        // check if the neighbor is valid. if thread id is smaller don't wait
        // if (dct.isValid(neighbor) && dct.isColored(neighbor))
        if (dct.isValid(neighbor))
        {
            color_state |= result[neighbor];
        }
        // add to invalid list
        else 
        {
            // std::cerr << "thread waiting: vertex " << vertex << " neighbor " << neighbor << std::endl;
            // if (vertex > neighbor)
            invalid_vertices.push_back(neighbor);
        }  

    }

    // check invalid list, wait till all invalid neighbors are colored
    // set state as WAITING
    // thread_table[thread_id].second = WAITING;
    while (true) {
        // traverse invalid list with iterator
        for (auto it = invalid_vertices.begin(); it != invalid_vertices.end();) {
            int neighbor = *it;
            // std::cerr << "thread waiting: vertex " << vertex << " neighbor " << neighbor << std::endl;
            // check if the neighbor is valid
            if (dct.isValid(neighbor)) {

                // check if the neighbor is colored
                color_state |= result[neighbor];

                // remove from invalid list
                it = invalid_vertices.erase(it);
            }
            else {
                // increment iterator
                it++;
            }
        }
        if (invalid_vertices.empty())
            break;
    }
    // set state as RUNNING
    // thread_table[thread_id].second = RUNNING;

    std::bitset<512> color;
    if (color_state == 0) {
        color = 1;
    } 
    else {
        color = (~color_state) & (increment_bitset(color_state));
    }

    result[vertex] = color;
    // dct.setColored(vertex);
    dct.setValid(vertex);
    // set color on result
    // result[vertex] = color;
}

// check conflict
// true: conflict exists
bool ParallelGraph::checkConflict() {
    int num = 0;
    int total_edges = 0;
    for (auto const& x : result) {
        int vertex = x.first;
        std::list<int> neighbors = adj[vertex];
        std::bitset<512> color = result[vertex];
        for (auto const& y : neighbors) {
            total_edges++;
            if (result[y] == color) {
                num++;
            }
        }
    }
    std::cerr << "Number of conflicts: " << num / 2 << std::endl;
    num = num / 2;
    total_edges = total_edges / 2;
    // ratio to total edges
    double ratio = (double) num / total_edges;
    ratio = ratio * 100;
    std::cerr << "Total edges: " << total_edges << std::endl;
    std::cerr << "Conflict ratio: " << ratio << "%" << std::endl;
    if (num > 0) {
        return true;
    }
    return false;
}

// check conflict
bool Graph::checkConflict(std::map<int, std::bitset<512>> result) {
    int num = 0;
    int total_edges = 0;
    for (auto const& x : result) {
        int vertex = x.first;
        std::list<int> neighbors = adj[vertex];
        std::bitset<512> color = result[vertex];
        for (auto const& y : neighbors) {
            total_edges++;
            if (result[y] == color) {
                num++;
            }
        }
    }
    std::cerr << "Number of conflicts: " << num / 2 << std::endl;
    num = num / 2;
    total_edges = total_edges / 2;
    // ratio to total edges
    double ratio = (double) num / total_edges;
    ratio = ratio * 100;
    std::cerr << "Total edges: " << total_edges << std::endl;
    std::cerr << "Conflict ratio: " << ratio << "%" << std::endl;
    if (num > 0) {
        return true;
    }
    return false;
}