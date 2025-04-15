#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <chrono>

// Include ParlayLib for parallelism
#include "parlaylib/include/parlay/primitives.h"
#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/sequence.h"
#include "parlaylib/include/parlay/utilities.h"

// Edge structure
struct Edge {
    int u, v;
    Edge(int uu = 0, int vv = 0) : u(uu), v(vv) {}
};

// Graph representation
class Graph {
private:
    int n; // Number of vertices
    std::vector<std::vector<int>> adj; // Adjacency list
    parlay::sequence<Edge> edges; // Edge list
    
public:
    Graph(int num_vertices) : n(num_vertices) {
        adj.resize(n);
    }
    
    // Add an undirected edge
    void add_edge(int u, int v) {
        if (u != v) { // No self-loops
            adj[u].push_back(v);
            adj[v].push_back(u);
            edges.push_back(Edge(u, v));
        }
    }
    
    // Generate Erdős–Rényi graph G(n,p)
    void generate_erdos_renyi(double p) {
        std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<> dist(0.0, 1.0);
        
        for (int u = 0; u < n; u++) {
            for (int v = u + 1; v < n; v++) {
                if (dist(gen) < p) {
                    add_edge(u, v);
                }
            }
        }
    }
    
    // Get number of edges
    int get_edge_count() const {
        return edges.size();
    }
    
    // Calculate average degree
    double calculate_average_degree() const {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += adj[i].size();
        }
        return sum / n;
    }
    
    // Estimate diameter using BFS from a randomly chosen node
    int estimate_diameter() const {
        if (n == 0) return 0;
        
        // Choose a random start vertex
        std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<> dist(0, n - 1);
        int start = dist(gen);
        
        // Find the farthest node from start
        std::vector<int> dist1(n, -1);
        int farthest1 = bfs(start, dist1);
        
        // Find the farthest node from farthest1
        std::vector<int> dist2(n, -1);
        int farthest2 = bfs(farthest1, dist2);
        
        // The distance from farthest1 to farthest2 is an estimate of the diameter
        return dist2[farthest2];
    }
    
    // Helper BFS function for diameter estimation
    int bfs(int start, std::vector<int>& distance) const {
        std::queue<int> q;
        q.push(start);
        distance[start] = 0;
        
        int farthest = start;
        int max_dist = 0;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (int v : adj[u]) {
                if (distance[v] == -1) {
                    distance[v] = distance[u] + 1;
                    q.push(v);
                    
                    if (distance[v] > max_dist) {
                        max_dist = distance[v];
                        farthest = v;
                    }
                }
            }
        }
        
        return farthest;
    }
    
    // Calculate clustering coefficient
    double calculate_clustering_coefficient() const {
        double total_coefficient = 0.0;
        
        parlay::parallel_for(0, n, [&](int u) {
            int possible_triangles = adj[u].size() * (adj[u].size() - 1) / 2;
            if (possible_triangles == 0) return;
            
            int triangles = 0;
            for (size_t i = 0; i < adj[u].size(); i++) {
                for (size_t j = i + 1; j < adj[u].size(); j++) {
                    int v = adj[u][i];
                    int w = adj[u][j];
                    
                    // Check if edge v-w exists
                    if (std::find(adj[v].begin(), adj[v].end(), w) != adj[v].end()) {
                        triangles++;
                    }
                }
            }
            
            double local_coef = static_cast<double>(triangles) / possible_triangles;
            
            #pragma omp atomic
            total_coefficient += local_coef;
        });
        
        return total_coefficient / n;
    }
    
    // Save graph to file in edge list format
    bool save_to_file(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return false;
        }
        
        // Write number of vertices and edges
        file << n << " " << edges.size() << std::endl;
        
        // Write edges
        for (const auto& edge : edges) {
            file << edge.u << " " << edge.v << std::endl;
        }
        
        file.close();
        return true;
    }
};

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <num_vertices> <probability> <output_file>" << std::endl;
        return 1;
    }
    
    int n = std::stoi(argv[1]);
    double p = std::stod(argv[2]);
    std::string output_file = argv[3];
    
    // Check valid parameters
    if (n <= 0 || p < 0.0 || p > 1.0) {
        std::cerr << "Error: Invalid parameters. n should be positive and p should be between 0 and 1." << std::endl;
        return 1;
    }
    
    std::cout << "Generating Erdős–Rényi graph G(" << n << "," << p << ")..." << std::endl;
    
    // Create and generate the graph
    Graph graph(n);
    graph.generate_erdos_renyi(p);
    
    // Compute metrics
    int edge_count = graph.get_edge_count();
    double avg_degree = graph.calculate_average_degree();
    int diameter = graph.estimate_diameter();
    double clustering = graph.calculate_clustering_coefficient();
    
    // Print results
    std::cout << "Graph metrics:" << std::endl;
    std::cout << "  Number of vertices: " << n << std::endl;
    std::cout << "  Number of edges: " << edge_count << std::endl;
    std::cout << "  Average degree: " << avg_degree << std::endl;
    std::cout << "  Estimated diameter: " << diameter << std::endl;
    std::cout << "  Clustering coefficient: " << clustering << std::endl;
    
    // Save graph to file
    if (graph.save_to_file(output_file)) {
        std::cout << "Graph saved to " << output_file << std::endl;
    } else {
        std::cerr << "Failed to save graph." << std::endl;
        return 1;
    }
    
    return 0;
}

// In current location:
// compile
// g++-12 -std=c++17 -O2 GraphSim.cpp -o ergen 