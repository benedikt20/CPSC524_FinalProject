#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <string>
#include <limits>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <atomic>

// Include ParlayLib for parallelism
#include "parlaylib/include/parlay/primitives.h"
#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/sequence.h"
#include "parlaylib/include/parlay/utilities.h"

// Graph representation using adjacency list
struct Edge {
    int to;
    int weight;
};

using Graph = std::vector<std::vector<Edge>>;

// Sequential Dijkstra's algorithm
std::vector<int> dijkstra_sequential(const Graph& graph, int source) {
    int n = graph.size();
    std::vector<int> dist(n, std::numeric_limits<int>::max());
    dist[source] = 0;
    
    // Priority queue with min-heap property
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    pq.push({0, source});
    
    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();
        
        // Skip outdated distance information
        if (d > dist[u]) continue;
        
        for (const auto& edge : graph[u]) {
            int v = edge.to;
            int weight = edge.weight;
            
            if (dist[u] != std::numeric_limits<int>::max() && 
                dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pq.push({dist[v], v});
            }
        }
    }
    
    return dist;
}

// Parallel Dijkstra's algorithm using ParlayLib
std::vector<int> dijkstra_parallel(const Graph& graph, int source) {
    int n = graph.size();
    std::vector<std::atomic<int>> dist(n);
    std::vector<bool> visited(n, false);
    
    // Initialize distances
    for (int i = 0; i < n; i++) {
        dist[i].store(std::numeric_limits<int>::max());
    }
    dist[source].store(0);

    for (int i = 0; i < n; i++) {
        // Find the vertex with minimum distance value
        int u = -1;
        int min_dist = std::numeric_limits<int>::max();
        
        for (int j = 0; j < n; j++) {
            if (!visited[j] && dist[j].load() < min_dist) {
                min_dist = dist[j].load();
                u = j;
            }
        }
        
        if (min_dist == std::numeric_limits<int>::max()) break;
        
        visited[u] = true;
        
        // Relax all the adjacent vertices in parallel
        if (!graph[u].empty()) {
            parlay::parallel_for(0, graph[u].size(), [&](size_t j) {
                int v = graph[u][j].to;
                int weight = graph[u][j].weight;
                
                if (!visited[v]) {
                    int current_dist = dist[u].load();
                    if (current_dist != std::numeric_limits<int>::max()) {
                        int new_dist = current_dist + weight;
                        int old_dist = dist[v].load();
                        
                        // Use atomic compare-and-swap for thread safety
                        if (new_dist < old_dist) {
                            int expected = old_dist;
                            dist[v].compare_exchange_strong(expected, new_dist);
                        }
                    }
                }
            });
        }
    }
    
    // Convert atomic vector to regular vector for the return
    std::vector<int> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = dist[i].load();
    }
    
    return result;
}

// Function to load a graph from a file
Graph load_graph(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }
    
    int n, m;
    file >> n >> m;
    
    Graph graph(n);
    
    for (int i = 0; i < m; i++) {
        int u, v;
        file >> u >> v;
        // For simplicity, we use a weight of 1 for all edges
        graph[u].push_back({v, 1});
        graph[v].push_back({u, 1}); // For undirected graph
    }
    
    return graph;
}

// Function to run benchmark and return execution time in milliseconds
double benchmark_algorithm(const Graph& graph, int source, bool use_parallel) {
    auto start = std::chrono::high_resolution_clock::now();
    
    if (use_parallel) {
        dijkstra_parallel(graph, source);
    } else {
        dijkstra_sequential(graph, source);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    
    return duration.count();
}

// Function to extract graph parameters from filename
std::pair<int, double> extract_graph_params(const std::string& filename) {
    // Assuming filename format like "graph_X_n100_p005.txt"
    int n = 0;
    double p = 0.0;
    
    size_t n_pos = filename.find("_n");
    size_t p_pos = filename.find("_p");
    
    if (n_pos != std::string::npos && p_pos != std::string::npos) {
        std::string n_str = filename.substr(n_pos + 2, p_pos - n_pos - 2);
        std::string p_str = filename.substr(p_pos + 2, filename.find(".txt") - p_pos - 2);
        
        n = std::stoi(n_str);
        
        // Handle p values with leading zeros like "001" for 0.01
        double divisor = std::pow(10, p_str.length());
        p = std::stod(p_str) / divisor;
    }
    
    return {n, p};
}

int main() {
    std::vector<std::string> graph_files = {
        "graphs/graph_A_n100_p005.txt",
        "graphs/graph_B_n100_p02.txt",
        "graphs/graph_C_n1000_p001.txt",
        "graphs/graph_D_n1000_p01.txt",
        "graphs/graph_E_n10000_p0001.txt",
        "graphs/graph_F_n10000_p001.txt",
        "graphs/graph_G_n50000_p00002.txt",
        "graphs/graph_H_n50000_p0002.txt"
    };
    
    // Open CSV file for results
    std::ofstream csv_file("dijkstra.csv");
    if (!csv_file.is_open()) {
        std::cerr << "Error creating output file" << std::endl;
        return 1;
    }
    
    // Write CSV header
    csv_file << "Graph,Nodes,Edge_Probability,Sequential_Time_ms,Parallel_Time_ms,Speedup\n";
    
    // Run benchmark for each graph
    for (const auto& filename : graph_files) {
        std::cout << "Processing " << filename << "..." << std::endl;
        
        try {
            Graph graph = load_graph(filename);
            int source = 0; // Use node 0 as the source
            
            // Extract graph parameters
            auto [n, p] = extract_graph_params(filename);
            char graph_id = filename.substr(filename.find("graph_") + 6, 1)[0]; // Extract the graph ID (A, B, C, etc.)
            
            // Run sequential algorithm (average of 3 runs)
            double sequential_time = 0;
            for (int i = 0; i < 3; i++) {
                sequential_time += benchmark_algorithm(graph, source, false);
            }
            sequential_time /= 3;
            
            // Run parallel algorithm (average of 3 runs)
            double parallel_time = 0;
            for (int i = 0; i < 3; i++) {
                parallel_time += benchmark_algorithm(graph, source, true);
            }
            parallel_time /= 3;
            
            // Calculate speedup
            double speedup = sequential_time / parallel_time;
            
            // Write results to CSV
            csv_file << "Graph_" << graph_id << "," << n << "," << p << ","
                     << std::fixed << std::setprecision(5) << sequential_time << ","
                     << parallel_time << "," << speedup << "\n";
            
            std::cout << "Graph " << graph_id << " (N=" << n << ", p=" << p << "):\n";
            std::cout << "  Sequential: " << sequential_time << " ms\n";
            std::cout << "  Parallel: " << parallel_time << " ms\n";
            std::cout << "  Speedup: " << speedup << "x\n";
            
        } catch (const std::exception& e) {
            std::cerr << "Error processing " << filename << ": " << e.what() << std::endl;
        }
    }
    
    csv_file.close();
    std::cout << "Results saved to dijkstra.csv" << std::endl;
    
    return 0;
}

// What we're currently doing:
// * No data partitioning: The distance array and adjacency list are shared among all threads.
// * Only edge relaxation is parallelized: We're only parallelizing the edge relaxation step (using parlay::parallel_for on the edges of the current vertex), but not the more critical parts of the algorithm.
// * Sequential bottleneck: We're still performing the minimum-distance vertex selection sequentially, which is a major bottleneck.


// compiler:
// g++-12 -std=c++17 -O2 dijkstra.cpp -o djstra 
// g++-12 -std=c++17 -O2 dijkstra.cpp -I /Users/BFF/parlaylib -o djstra 

// run:
// ./djstra graph_100_01.txt