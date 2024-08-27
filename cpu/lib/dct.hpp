
#ifndef DCT_HPP
#define DCT_HPP

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
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
    std::unordered_map<int, ValidColoredTuple> table;

public:
    // constructor
    DataConflictTable() {
        table.clear();
    }

    // insert vertex into data conflict table
    void insert(int vertex) {

        table[vertex].valid.store(true);
        table[vertex].colored.store(false);
    }

    // remove vertex from data conflict table
    void remove(int vertex) {

        table.erase(vertex);
    }

    // is vertex valid
    bool isValid(int vertex) {

        // if it is colored, and valid
        return table[vertex].valid.load();
    }

    // set vertex to valid
    void setValid(int vertex) {

        table[vertex].valid.store(true);
    }

    // set vertex to invalid
    void setInvalid(int vertex) {

        table[vertex].valid.store(false);
    }    

    // is vertex colored
    bool isColored(int vertex) {

        // if it is colored, and valid
        return table[vertex].colored.load();
    }

    // set vertex to colored
    void setColored(int vertex) {
        table[vertex].colored.store(true);
    }

    // set vertex to uncolored
    void setUncolored(int vertex) {
        table[vertex].colored.store(false);
    }

    // get the size of the data conflict table
    int size() {
        return table.size();
    }

    // reserve
    void reserve(size_t size) {
        table.reserve(size);
    }
};

#endif