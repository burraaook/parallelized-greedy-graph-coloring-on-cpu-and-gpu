#include "../lib/graph.hpp"

// measure time
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <bitset>
#include <string>
#include <vector>

// define the number of tests
#define NUM_TESTS 50

// filenames
std::string bitcolor_test_filename = "output/bitcolor_test.txt";
std::string bitcolor_parallel_results_filename = "output/bitcolor_parallel_results.txt";

int basic_greedy(std::string filename, std::ofstream &file);
int bit_wise_greedy_normal(std::string filename, std::ofstream &file);
void write_performance_result(std::ofstream &file, std::string algorithm, std::string dataset, int time, int colors, int conflict=0);
int parallel_graph_coloring(std::string filename, int num_threads, std::ofstream &file, bool sortOption=true);
void write_parallel_result(std::ofstream &file, std::string algorithm, int threads, std::string dataset, int time, int colors, int conflict);

void run_tests();

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        std::cerr << "Usage: ./bitcolor_cpu <algorithm_name> <dataset_path>" << std::endl;
        return 1;
    }


    std::string algorithm = argv[1];

    if (algorithm == "run_tests") {
        run_tests();
    } 

    // -h
    else if (algorithm == "-h") {
        std::cout << "Usage: ./bitcolor_cpu <algorithm_name> <dataset_path>" << std::endl;
        std::cout << "Algorithms: basic_greedy, bitwise_greedy, p_bitcolor" << std::endl;
        std::cout << "To run tests: ./bitcolor_cpu run_tests" << std::endl;
        std::cout << "p_bitcolor options: " << std::endl;
        std::cout << "[NUM_THREADS] [sort_option]" << std::endl;
        std::cout << "sort_option: sort_yes, sort_no" << std::endl;
    }

    else if (algorithm == "basic_greedy") {
        if (argc != 3) {
            std::cerr << "Usage: ./bitcolor_cpu basic_greedy <dataset_path>" << std::endl;
            return 1;
        }

        std::string dataset = argv[2];
        std::ofstream bitcolor_test_file(bitcolor_test_filename, std::ios::app);
        if (!bitcolor_test_file.is_open()) {
            std::cerr << "Error: could not open file " << bitcolor_test_filename << std::endl;
            return 1;
        }

        basic_greedy(dataset, bitcolor_test_file);
        bitcolor_test_file.close();
    } else if (algorithm == "bitwise_greedy") {
        if (argc != 3) {
            std::cerr << "Usage: ./bitcolor_cpu bitwise_greedy <dataset_path>" << std::endl;
            return 1;
        }

        std::string dataset = argv[2];
        std::ofstream bitcolor_test_file(bitcolor_test_filename, std::ios::app);
        if (!bitcolor_test_file.is_open()) {
            std::cerr << "Error: could not open file " << bitcolor_test_filename << std::endl;
            return 1;
        }

        bit_wise_greedy_normal(dataset, bitcolor_test_file);
        bitcolor_test_file.close();
    } else if (algorithm == "p_bitcolor") {
        if (argc != 5) {
            std::cerr << "Usage: ./bitcolor_cpu p_bitcolor <dataset_path> <num_threads> <sort_option>" << std::endl;
            return 1;
        }

        std::string dataset = argv[2];
        int num_threads = std::stoi(argv[3]);

        // check if the number of threads is greater than 1
        if (num_threads < 1) {
            std::cerr << "Error: number of threads must be greater than 1" << std::endl;
            return 1;
        }
        else if (num_threads > 64) {
            std::cerr << "Error: number of threads must be less than 64" << std::endl;
            return 1;
        }

        num_threads -= 1;
        std::ofstream bitcolor_parallel_results_file(bitcolor_parallel_results_filename, std::ios::app);
        if (!bitcolor_parallel_results_file.is_open()) {
            std::cerr << "Error: could not open file " << bitcolor_parallel_results_filename << std::endl;
            return 1;
        }
        
        // sort_yes, sort_no are the possibilities
        bool sortOption = false;
        if (std::string(argv[4]) == "sort_yes") {
            sortOption = true;
        }


        parallel_graph_coloring(dataset, num_threads, bitcolor_parallel_results_file, sortOption);
        bitcolor_parallel_results_file.close();
    } else {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    return 0;
}

void run_tests()
{
    std::vector<std::string> filenames = {
        "datasets/EF.txt",
        "datasets/CA.txt",
        "datasets/CL.txt",
        "datasets/CD.txt",
        "datasets/RC.txt",
        "datasets/RT.txt",
        "datasets/RP.txt"
    };

    std::ofstream bitcolor_test_file(bitcolor_test_filename, std::ios::app);
    if (!bitcolor_test_file.is_open()) {
        std::cerr << "Error: could not open file " << bitcolor_test_filename << std::endl;
        return;
    }

    std::ofstream bitcolor_parallel_results_file(bitcolor_parallel_results_filename, std::ios::app);
    if (!bitcolor_parallel_results_file.is_open()) {
        std::cerr << "Error: could not open file " << bitcolor_parallel_results_filename << std::endl;
        return;
    }

    for (auto const& filename : filenames) {
        for (int i = 0; i < NUM_TESTS; i++) {
            basic_greedy(filename, bitcolor_test_file);
            bit_wise_greedy_normal(filename, bitcolor_test_file);
            parallel_graph_coloring(filename, 11, bitcolor_test_file);
            bitcolor_test_file << std::endl; // Empty line between each iteration
        }

        // Parallel tests with different thread numbers
        std::vector<int> thread_numbers = {3, 7, 11, 15/*, 31*/};
        for (int threads : thread_numbers) {
            for (int i = 0; i < NUM_TESTS-40; i++) {
                parallel_graph_coloring(filename, threads, bitcolor_parallel_results_file);
            }
        }
    }

    bitcolor_test_file.close();
    bitcolor_parallel_results_file.close();
}

int basic_greedy(std::string filename, std::ofstream &file)
{
    std::cout << "\nalgorithm: basic_greedy" << std::endl;
    // std::cout << "dataset: " << filename << std::endl;

    Graph g(filename);

    // measure time
    std::cout << "Starting coloring" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::map<int, int> result = g.basicGreedyColoring();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << " milliseconds\n" << std::endl;
    std::cout << "Number of colors used: " << g.getNumColor() << std::endl << std::endl;

    write_performance_result(file, "basic_greedy", filename, duration.count(), g.getNumColor());

    return 0;
}

int bit_wise_greedy_normal(std::string filename, std::ofstream &file)
{
    std::cout << "\nalgorithm: bit_wise_greedy" << std::endl;
    // std::cout << "dataset: " << filename << std::endl;

    Graph g(filename);

    // measure time
    std::cout << "Starting coloring" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::map<int, std::bitset<512>> result = g.bitWiseColoring();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << " milliseconds" << std::endl << std::endl;
    int max_color = g.countColorsBitset(result);
    std::cout << "Number of colors used: " << max_color << std::endl;

    bool conflict = g.checkConflict(result);
    int conflict_int = conflict ? 1 : 0;

    // check conflict
    // std::cout << "Checking conflict: " << conflict << std::endl;

    write_performance_result(file, "bit_wise_greedy_normal", filename, duration.count(), max_color, conflict_int);

    return 0;
}

int parallel_graph_coloring(std::string filename, int num_threads, std::ofstream &file, bool sortOption)
{
    std::cout << "\nalgorithm: parallel_graph_coloring" << std::endl;
    // std::cout << "dataset: " << filename << std::endl;
    ParallelGraph g(filename, num_threads, sortOption);

    // measure time
    std::cout << "Starting coloring" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::map<int, std::bitset<512>> result = g.bitWiseColoring();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << " milliseconds" << std::endl << std::endl;
    int max_color = g.countColorsBitset(result);
    std::cout << "Number of colors used: " << max_color << std::endl;

    // conflict result
    bool conflict = g.checkConflict();

    // check conflict
    // std::cout << "Checking conflict: " << conflict << std::endl;

    int conflict_int = conflict ? 1 : 0;

    write_parallel_result(file, "parallel_graph_coloring", num_threads, filename, duration.count(), max_color, conflict_int);

    return 0;
}

void write_performance_result(std::ofstream &file, std::string algorithm, std::string dataset, int time, int colors, int conflict)
{
    file << algorithm << " " << dataset << " " << time << " " << colors << " " << conflict << std::endl;
}

void write_parallel_result(std::ofstream &file, std::string algorithm, int threads, std::string dataset, int time, int colors, int conflict)
{
    file << algorithm << " " << threads << " " << dataset << " " << time << " " << colors << " " << conflict << std::endl;
}
