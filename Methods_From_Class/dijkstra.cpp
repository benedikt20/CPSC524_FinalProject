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
#include <thread>
#include <mutex>

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

// Improved parallel Dijkstra's algorithm with partitioning
std::vector<int> dijkstra_parallel_partitioned(const Graph& graph, int source) {
    int n = graph.size();
    std::vector<int> dist(n, std::numeric_limits<int>::max());
    std::vector<bool> visited(n, false);
    dist[source] = 0;
    
    // Determine number of threads/partitions
    unsigned int num_threads = std::thread::hardware_concurrency();
    // Ensure we have at least 2 threads, but no more than the graph size
    num_threads = std::max(2u, std::min(num_threads, static_cast<unsigned int>(n / 10 + 1)));
    
    // Create thread-local distance arrays for each partition
    std::vector<std::vector<int>> local_dist(num_threads, std::vector<int>(n, std::numeric_limits<int>::max()));
    for (unsigned int t = 0; t < num_threads; t++) {
        local_dist[t][source] = 0;
    }
    
    // Mutex for thread synchronization when updating the global visited array
    std::mutex visited_mutex;
    
    // Partition size (number of vertices per thread)
    int partition_size = (n + num_threads - 1) / num_threads;
    
    for (int iter = 0; iter < n; iter++) {
        // Each thread finds min distance vertex in its partition
        std::vector<std::pair<int, int>> local_min_vertices(num_threads, {std::numeric_limits<int>::max(), -1});
        
        parlay::parallel_for(0, num_threads, [&](size_t t) {
            int start_idx = t * partition_size;
            int end_idx = std::min(start_idx + partition_size, n);
            
            for (int i = start_idx; i < end_idx; i++) {
                if (!visited[i] && dist[i] < local_min_vertices[t].first) {
                    local_min_vertices[t] = {dist[i], i};
                }
            }
        });
        
        // Find the global minimum across all partitions
        auto global_min = *std::min_element(local_min_vertices.begin(), local_min_vertices.end());
        int u = global_min.second;
        
        if (u == -1 || global_min.first == std::numeric_limits<int>::max()) {
            break;  // All reachable vertices processed
        }
        
        // Mark vertex as visited
        visited[u] = true;
        
        // Partition the edges of the current vertex for parallel processing
        std::vector<Edge> edges = graph[u];
        int edges_per_thread = (edges.size() + num_threads - 1) / num_threads;
        
        parlay::parallel_for(0, num_threads, [&](size_t t) {
            int start_edge = t * edges_per_thread;
            int end_edge = std::min(start_edge + edges_per_thread, static_cast<int>(edges.size()));
            
            // Process assigned edges
            for (int i = start_edge; i < end_edge; i++) {
                int v = edges[i].to;
                int weight = edges[i].weight;
                
                if (!visited[v]) {
                    if (dist[u] != std::numeric_limits<int>::max() && 
                        dist[u] + weight < local_dist[t][v]) {
                        local_dist[t][v] = dist[u] + weight;
                    }
                }
            }
        });
        
        // Merge the local distance arrays to get the global minimum distances
        parlay::parallel_for(0, n, [&](size_t v) {
            if (!visited[v]) {
                int min_dist = dist[v];
                for (unsigned int t = 0; t < num_threads; t++) {
                    min_dist = std::min(min_dist, local_dist[t][v]);
                }
                dist[v] = min_dist;
                
                // Reset local distance arrays for next iteration
                for (unsigned int t = 0; t < num_threads; t++) {
                    local_dist[t][v] = min_dist;
                }
            }
        });
    }
    
    return dist;
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
        dijkstra_parallel_partitioned(graph, source);
    } else {
        dijkstra_sequential(graph, source);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    
    return duration.count();
}

// Function to validate that both algorithms produce the same results
bool validate_results(const Graph& graph, int source) {
    auto seq_result = dijkstra_sequential(graph, source);
    auto par_result = dijkstra_parallel_partitioned(graph, source);
    
    if (seq_result.size() != par_result.size()) {
        return false;
    }
    
    for (size_t i = 0; i < seq_result.size(); i++) {
        if (seq_result[i] != par_result[i]) {
            return false;
        }
    }
    
    return true;
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
    csv_file << "Graph,Nodes,Edge_Probability,Sequential_Time_ms,Parallel_Time_ms,Speedup,Valid_Results\n";
    
    // Run benchmark for each graph
    for (const auto& filename : graph_files) {
        std::cout << "Processing " << filename << "..." << std::endl;
        
        try {
            Graph graph = load_graph(filename);
            int source = 0; // Use node 0 as the source
            
            // Extract graph parameters
            auto [n, p] = extract_graph_params(filename);
            char graph_id = filename.substr(filename.find("graph_") + 6, 1)[0]; // Extract the graph ID (A, B, C, etc.)
            
            // Validate results
            bool is_valid = validate_results(graph, source);
            std::cout << "  Validation: " << (is_valid ? "Passed" : "Failed") << std::endl;
            
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
                     << std::fixed << std::setprecision(2) << sequential_time << ","
                     << parallel_time << "," << speedup << "," << (is_valid ? "Yes" : "No") << "\n";
            
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
// compiler:
// g++-12 -std=c++17 -O2 dijkstra.cpp -o djstra 
// g++-12 -std=c++17 -O2 dijkstra.cpp -I /Users/BFF/parlaylib -o djstra 

// run:
// ./djstra graph_100_01.txt