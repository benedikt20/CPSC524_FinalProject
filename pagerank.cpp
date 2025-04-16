#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cmath>
#include <unordered_map>
#include <iomanip>

// Include ParlayLib for parallelism
#include "parlaylib/include/parlay/primitives.h"
#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/sequence.h"
#include "parlaylib/include/parlay/utilities.h"

// Constants for PageRank
const double DAMPING = 0.85;     // Damping factor
const double EPSILON = 1e-6;     // Convergence threshold
const int MAX_ITERATIONS = 100;  // Maximum number of iterations

// Store results
struct BenchmarkResult {
    std::string graph_name;
    int vertices;
    int edges;
    double avg_degree;
    double sequential_time;
    double parallel_time;
    double speedup;
    int seq_iterations;
    int par_iterations;
    double avg_diff; // difference between sequential and parallel results
};


// Graph representation
class Graph {
private:
    int n;                        // Number of vertices
    int m;                        // Number of edges
    std::vector<std::vector<int>> out_edges;  // Outgoing edges
    std::vector<std::vector<int>> in_edges;   // Incoming edges
    std::vector<int> out_degree;  // Outgoing degrees

public:
    Graph() : n(0), m(0) {}

    // Load graph from edge list file
    bool load_from_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        // Read number of vertices and edges
        file >> n >> m;
        
        // Initialize data structures
        out_edges.resize(n);
        in_edges.resize(n);
        out_degree.resize(n, 0);

        // Read edges
        for (int i = 0; i < m; i++) {
            int u, v;
            file >> u >> v;
            
            // Add edge
            out_edges[u].push_back(v);
            in_edges[v].push_back(u);
            out_degree[u]++;
        }

        file.close();
        return true;
    }

    // Get basic graph info
    void print_info() const {
        std::cout << "Graph Info:" << std::endl;
        std::cout << "  Vertices: " << n << std::endl;
        std::cout << "  Edges: " << m << std::endl;
        
        // Calculate average degree
        double avg_degree = static_cast<double>(m * 2) / n;
        std::cout << "  Average degree: " << avg_degree << std::endl;
    }
    
    // Get average degree
    double get_avg_degree() const {
        return static_cast<double>(m * 2) / n;
    }
    
    // Get number of vertices
    int get_vertices() const {
        return n;
    }
    
    // Get number of edges
    int get_edges() const {
        return m;
    }

    // Sequential PageRank implementation
    std::vector<double> pagerank_sequential(int max_iter = MAX_ITERATIONS, double epsilon = EPSILON, int* iterations = nullptr) {
        // Initialize PageRank scores
        std::vector<double> scores(n, 1.0 / n);
        std::vector<double> new_scores(n, 0.0);
        
        // Random teleportation value
        double random_jump = (1.0 - DAMPING) / n;
        
        // Main loop
        int iter;
        for (iter = 0; iter < max_iter; iter++) {
            // Reset new_scores
            std::fill(new_scores.begin(), new_scores.end(), 0.0);
            
            // Update scores
            for (int u = 0; u < n; u++) {
                if (out_degree[u] > 0) {
                    double contribution = DAMPING * scores[u] / out_degree[u];
                    for (int v : out_edges[u]) {
                        new_scores[v] += contribution;
                    }
                } else {
                    // Distribute score evenly for dangling nodes
                    for (int v = 0; v < n; v++) {
                        new_scores[v] += DAMPING * scores[u] / n;
                    }
                }
            }
            
            // Add random jump component
            for (int v = 0; v < n; v++) {
                new_scores[v] += random_jump;
            }
            
            // Check for convergence
            double diff = 0.0;
            for (int v = 0; v < n; v++) {
                diff += std::abs(new_scores[v] - scores[v]);
            }
            
            // Swap scores
            scores.swap(new_scores);
            
            if (diff < epsilon) {
                std::cout << "  Sequential PageRank converged after " << iter + 1 << " iterations." << std::endl;
                break;
            }
            
            if (iter == max_iter - 1) {
                std::cout << "  Sequential PageRank reached maximum iterations." << std::endl;
            }
        }
        
        if (iterations) {
            *iterations = iter + 1;
        }
        
        return scores;
    }

    // Parallel PageRank implementation using ParlayLib
    parlay::sequence<double> pagerank_parallel(int max_iter = MAX_ITERATIONS, double epsilon = EPSILON, int* iterations = nullptr) {
        // Initialize PageRank scores
        parlay::sequence<double> scores(n, 1.0 / n);
        parlay::sequence<double> new_scores(n, 0.0);
        
        // Random teleportation value
        double random_jump = (1.0 - DAMPING) / n;
        
        // Main loop
        int iter;
        for (iter = 0; iter < max_iter; iter++) {
            // Reset new_scores
            parlay::parallel_for(0, n, [&](size_t i) {
                new_scores[i] = 0.0;
            });
            
            // Process contributions from vertices with out-edges
            parlay::parallel_for(0, n, [&](size_t u) {
                if (out_degree[u] > 0) {
                    double contribution = DAMPING * scores[u] / out_degree[u];
                    for (int v : out_edges[u]) {
                        // Atomic add to avoid race conditions
                        std::atomic_ref<double>(new_scores[v]).fetch_add(contribution, std::memory_order_relaxed);
                    }
                } else {
                    // Distribute score evenly for dangling nodes
                    double dangling_contrib = DAMPING * scores[u] / n;
                    parlay::parallel_for(0, n, [&](size_t v) {
                        std::atomic_ref<double>(new_scores[v]).fetch_add(dangling_contrib, std::memory_order_relaxed);
                    });
                }
            });
            
            // Add random jump component
            parlay::parallel_for(0, n, [&](size_t v) {
                new_scores[v] += random_jump;
            });
            
            // Check for convergence
            double diff = parlay::reduce(
                parlay::tabulate(n, [&](size_t i) { return std::abs(new_scores[i] - scores[i]); }),
                0.0, std::plus<double>()
            );
            
            // Swap scores
            scores.swap(new_scores);
            
            if (diff < epsilon) {
                std::cout << "  Parallel PageRank converged after " << iter + 1 << " iterations." << std::endl;
                break;
            }
            
            if (iter == max_iter - 1) {
                std::cout << "  Parallel PageRank reached maximum iterations." << std::endl;
            }
        }
        
        if (iterations) {
            *iterations = iter + 1;
        }
        
        return scores;
    }
    
    // Get number of vertices
    int num_vertices() const {
        return n;
    }
};

// Utility function to print top PageRank scores
template<typename T>
void print_top_scores(const T& scores, int num_to_show = 10) {
    // Create a vector of (score, vertex) pairs
    std::vector<std::pair<double, int>> ranked;
    for (int i = 0; i < static_cast<int>(scores.size()); i++) {
        ranked.push_back({scores[i], i});
    }
    
    // Sort by score in descending order
    std::sort(ranked.begin(), ranked.end(), 
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // Print top scores
    std::cout << "Top " << std::min(num_to_show, static_cast<int>(ranked.size())) 
              << " PageRank scores:" << std::endl;
    for (int i = 0; i < std::min(num_to_show, static_cast<int>(ranked.size())); i++) {
        std::cout << "  Vertex " << ranked[i].second << ": " 
                  << std::fixed << std::setprecision(6) << ranked[i].first << std::endl;
    }
}

// Save top PageRank scores to CSV file
template<typename T>
void save_top_scores_to_csv(const std::string& filename, const std::string& graph_name, 
                           const T& scores, int num_to_show = 10) {
    // Create a vector of (score, vertex) pairs
    std::vector<std::pair<double, int>> ranked;
    for (int i = 0; i < static_cast<int>(scores.size()); i++) {
        ranked.push_back({scores[i], i});
    }
    
    // Sort by score in descending order
    std::sort(ranked.begin(), ranked.end(), 
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // Open file for writing (append mode)
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    
    // Write top scores
    for (int i = 0; i < std::min(num_to_show, static_cast<int>(ranked.size())); i++) {
        file << graph_name << "," << ranked[i].second << "," 
             << std::fixed << std::setprecision(6) << ranked[i].first << std::endl;
    }
    
    file.close();
}

// Benchmark function to run PageRank and measure time
BenchmarkResult benchmark_pagerank(const std::string& graph_file, const std::string& graph_name) {
    std::cout << "\nBenchmarking " << graph_name << " from file " << graph_file << std::endl;
    
    BenchmarkResult result;
    result.graph_name = graph_name;
    
    // Load graph
    Graph g;
    if (!g.load_from_file(graph_file)) {
        std::cerr << "Failed to load graph. Skipping..." << std::endl;
        return result;
    }
    
    g.print_info();
    
    // Store graph info in result
    result.vertices = g.get_vertices();
    result.edges = g.get_edges();
    result.avg_degree = g.get_avg_degree();
    
    // Run sequential PageRank and measure time
    int seq_iters = 0;
    auto start_seq = std::chrono::high_resolution_clock::now();
    auto seq_scores = g.pagerank_sequential(MAX_ITERATIONS, EPSILON, &seq_iters);
    auto end_seq = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> seq_time = end_seq - start_seq;
    std::cout << "Sequential PageRank completed in " << seq_time.count() << " seconds." << std::endl;
    
    // Store sequential results
    result.sequential_time = seq_time.count();
    result.seq_iterations = seq_iters;
    
    // Run parallel PageRank and measure time
    int par_iters = 0;
    auto start_par = std::chrono::high_resolution_clock::now();
    auto par_scores = g.pagerank_parallel(MAX_ITERATIONS, EPSILON, &par_iters);
    auto end_par = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> par_time = end_par - start_par;
    std::cout << "Parallel PageRank completed in " << par_time.count() << " seconds." << std::endl;
    
    // Store parallel results
    result.parallel_time = par_time.count();
    result.par_iterations = par_iters;
    
    // Calculate speedup
    double speedup = seq_time.count() / par_time.count();
    std::cout << "Speedup: " << speedup << "x" << std::endl;
    result.speedup = speedup;
    
    // Verify results match (approximately)
    double total_diff = 0.0;
    for (int i = 0; i < g.num_vertices(); i++) {
        total_diff += std::abs(seq_scores[i] - par_scores[i]);
    }
    double avg_diff = total_diff / g.num_vertices();
    std::cout << "Average difference between sequential and parallel results: " 
              << std::scientific << avg_diff << std::endl;
    result.avg_diff = avg_diff;
    
    // Print top scores from parallel version
    print_top_scores(par_scores, 5);
    
    // Save top scores to CSV
    save_top_scores_to_csv("pagerank_scores.csv", graph_name, par_scores, 10);
    
    std::cout << "-----------------------------------------------------" << std::endl;
    
    return result;
}

// Save benchmark results to CSV file
void save_results_to_csv(const std::string& filename, const std::vector<BenchmarkResult>& results) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    
    // Write header
    file << "Graph,Vertices,Edges,AvgDegree,SeqTime,ParTime,Speedup,SeqIters,ParIters,AvgDiff" << std::endl;
    
    // Write data
    for (const auto& r : results) {
        file << r.graph_name << ","
             << r.vertices << ","
             << r.edges << ","
             << std::fixed << std::setprecision(2) << r.avg_degree << ","
             << std::fixed << std::setprecision(6) << r.sequential_time << ","
             << r.parallel_time << ","
             << r.speedup << ","
             << r.seq_iterations << ","
             << r.par_iterations << ","
             << std::scientific << r.avg_diff << std::endl;
    }
    
    file.close();
    std::cout << "Results saved to " << filename << std::endl;
}

int main(int argc, char** argv) {
    std::vector<std::pair<std::string, std::string>> graphs = {
        {"graphs/graph_A_n100_p005.txt", "Graph A (N=100, p=0.05)"},
        {"graphs/graph_B_n100_p02.txt", "Graph B (N=100, p=0.2)"},
        {"graphs/graph_C_n1000_p001.txt", "Graph C (N=1000, p=0.01)"},
        {"graphs/graph_D_n1000_p01.txt", "Graph D (N=1000, p=0.1)"},
        {"graphs/graph_E_n10000_p0001.txt", "Graph E (N=10000, p=0.001)"},
        {"graphs/graph_F_n10000_p001.txt", "Graph F (N=10000, p=0.01)"},
        {"graphs/graph_G_n50000_p00002.txt", "Graph G (N=50000, p=0.0002)"},
        {"graphs/graph_H_n50000_p0002.txt", "Graph H (N=50000, p=0.002)"}
    };
    
    // Create header for scores CSV
    {
        std::ofstream file("pagerank_scores.csv");
        if (file.is_open()) {
            file << "Graph,VertexID,Score" << std::endl;
            file.close();
        }
    }
    
    // Run benchmarks for all graphs
    std::vector<BenchmarkResult> results;
    for (const auto& [file, name] : graphs) {
        BenchmarkResult result = benchmark_pagerank(file, name);
        results.push_back(result);
    }
    
    // Save results to CSV
    save_results_to_csv("pagerank_results.csv", results);
    
    // Generate summary
    std::cout << "\nPerformance Summary:" << std::endl;
    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << "| Graph Name                      | Seq Time | Par Time | Speedup |" << std::endl;
    std::cout << "---------------------------------------------------------------" << std::endl;
    for (const auto& r : results) {
        std::cout << "| " << std::left << std::setw(30) << r.graph_name
                  << " | " << std::right << std::setw(8) << std::fixed << std::setprecision(4) << r.sequential_time
                  << " | " << std::setw(8) << r.parallel_time
                  << " | " << std::setw(7) << r.speedup << " |" << std::endl;
    }
    std::cout << "---------------------------------------------------------------" << std::endl;
    
    return 0;
}

// compiler:
// g++-12 -std=c++17 -O2 pagerank.cpp -o pgrank 