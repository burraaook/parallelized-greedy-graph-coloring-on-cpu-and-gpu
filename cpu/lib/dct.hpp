
#ifndef DCT_HPP
#define DCT_HPP

#include <iostream>
#include <vector>
#include <map>
#include <shared_mutex>
#include <atomic>

// struct to hold valid, colored tuple
struct ValidColoredTuple
{
    // atomic
    std::atomic<bool> valid;
    std::atomic<bool> colored;
};

// data conflict table
// vertex id, valid: valid means is the vertex is colored or not
// <vertex_id, validColoredTuple>
// valid: true, invalid: false
class DataConflictTable
{
private:
    // data conflict table
    // vertex id -> valid, colored
    std::map<int, ValidColoredTuple> table;

    // read write mutex
    // mutable std::shared_mutex rw_mutex;
public:
    // constructor
    DataConflictTable() {
        table.clear();
    }

    // insert vertex into data conflict table
    void insert(int vertex) {

        // write lock
        // std::unique_lock<std::shared_mutex> lock(rw_mutex);
        table[vertex].valid.store(true);
        table[vertex].colored.store(false);
    }

    // remove vertex from data conflict table
    void remove(int vertex) {

        // write lock
        // std::unique_lock<std::shared_mutex> lock(rw_mutex);
        table.erase(vertex);
    }

    // is vertex valid
    bool isValid(int vertex) {

        // read lock
        // std::shared_lock<std::shared_mutex> lock(rw_mutex);
        
        // if it is colored, and valid
        return table[vertex].valid.load();
    }

    // set vertex to valid
    void setValid(int vertex) {

        // write lock
        // std::unique_lock<std::shared_mutex> lock(rw_mutex);
        table[vertex].valid.store(true);
    }

    // set vertex to invalid
    void setInvalid(int vertex) {

        // write lock
        // std::unique_lock<std::shared_mutex> lock(rw_mutex);
        table[vertex].valid.store(false);
    }    

    // is vertex colored
    bool isColored(int vertex) {

        // read lock
        // std::shared_lock<std::shared_mutex> lock(rw_mutex);
        
        // if it is colored, and valid
        return table[vertex].colored.load();
    }

    // set vertex to colored
    void setColored(int vertex) {

        // write lock
        // std::unique_lock<std::shared_mutex> lock(rw_mutex);
        table[vertex].colored.store(true);
    }

    // set vertex to uncolored
    void setUncolored(int vertex) {

        // write lock
        // std::unique_lock<std::shared_mutex> lock(rw_mutex);
        table[vertex].colored.store(false);
    }

    // get the size of the data conflict table
    int size() {

        // read lock
        // std::shared_lock<std::shared_mutex> lock(rw_mutex);
        return table.size();
    }
};

#endif