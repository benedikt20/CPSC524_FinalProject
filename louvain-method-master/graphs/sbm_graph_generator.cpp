// File: sbm_graph_generator.cpp
// -- Graph generator from the stochastic block model
// -----------------------------------------------------------------------------
// Simulation of a graph consisting of edge from node i to j with associated weight w
// with the line format: i j w

// This script generates four graphs with 10,000, 50,000, 100,000 and 200,000 nodes with:
// and 50 communities, high intra-community edge probability of pi=0.01, and sparse inter-community edges of po=0.001.
// These parameters are fixed in the script but can be modified.

// The resulting graph is saved as SBM_n<NODES>_k50.txt 

// This program must not be distributed without agreement of the above mentionned authors.
// -----------------------------------------------------------------------------
// Author   : Benedikt Farag
// Email    : benedikt.farag@yale.edu
// Location : New Haven CT
// Time	    : April 2025
// -----------------------------------------------------------------------------
// To compile  : g++-12 -O2 -o sbm_gen sbm_graph_generator.cpp
// To execute  : ./sbm_gen 


#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <filesystem>
#include <cstdint>

namespace fs = std::filesystem;

// Generate SBM graph with edge weights
void generate_sbm(int N, int K, double p_in, double p_out, const std::string& base_name) {
    std::vector<int> community(N);
    int nodes_per_block = N / K;

    // Assign communities to nodes
    for (int i = 0; i < N; ++i)
        community[i] = i / nodes_per_block;

    // Output files: one for the graph and another for weights
    std::ofstream edge_file(base_name + ".txt");

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> edge_decider(0.0, 1.0);
    std::uniform_real_distribution<float> weight_gen(0.1f, 1.0f);  // weights between 0.1 and 1.0

    // Generate edges and their weights
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double prob = (community[i] == community[j]) ? p_in : p_out;
            if (edge_decider(gen) < prob) {
                float w = weight_gen(gen); // Generate a random weight for the edge
                edge_file << i << " " << j << " " << w << "\n";  // Write edge and weight in a single line
            }
        }
    }

    edge_file.close();
    std::cout << "Saved: " << base_name << ".txt\n";
}

int main() {
    fs::create_directory("graphs");

    // Define graph sizes and other parameters
    std::vector<int> sizes = {10000, 50000, 100000, 200000};
    int K = 50;
    double p_in = 0.01;
    double p_out = 0.001;

    // Generate graphs for each size
    for (int N : sizes) {
        // to save in a subfolder 'graphs', use: "graphs/SBM_n" here below instead of "SBM_n"
        std::string base = "SBM_n" + std::to_string(N) + "_k" + std::to_string(K);
        generate_sbm(N, K, p_in, p_out, base);
    }

    return 0;
}
