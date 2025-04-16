// Watts–Strogatz Graph Generator (small world graph like in biology)
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>

using namespace std;

struct Edge {
    int u, v;
    Edge(int u, int v) : u(u), v(v) {}
};

class Graph {
public:
    int n;  // number of nodes
    vector<vector<int>> adj;  // adjacency list
    vector<Edge> edges;       // list of edges

    Graph(int n) : n(n) {
        adj.resize(n);
    }
    
    // Add an undirected edge (if not a self-loop).
    void add_edge(int u, int v) {
        if(u == v) return;
        // To keep things simple, we assume no duplicate check here.
        adj[u].push_back(v);
        adj[v].push_back(u);
        edges.push_back(Edge(u, v));
    }
    
    // Generate Watts–Strogatz graph.
    // k: Each node is connected to k nearest neighbors (must be even).
    // beta: rewiring probability.
    void generate_watts_strogatz(int k, double beta) {
        if(k % 2 != 0) {
            cerr << "Error: k must be an even number." << endl;
            return;
        }
        // Step 1: Build the initial ring lattice.
        for (int i = 0; i < n; i++) {
            for (int j = 1; j <= k / 2; j++) {
                int neighbor = (i + j) % n;
                // To avoid duplicate edges, add only if i < neighbor.
                if (i < neighbor) {
                    add_edge(i, neighbor);
                }
            }
        }
        
        // Random number generator.
        mt19937 gen(chrono::system_clock::now().time_since_epoch().count());
        uniform_real_distribution<> dis(0.0, 1.0);
        
        // Step 2: Rewire edges with probability beta.
        // We traverse each node's neighbors where (u,v) with u < v.
        for (int u = 0; u < n; u++) {
            // Create a copy of the neighbor list to iterate over
            vector<int> neighbors = adj[u];
            for (int v : neighbors) {
                if (v > u) {  // process each edge only once.
                    if (dis(gen) < beta) {
                        // Attempt to choose a new neighbor w not equal to u and not already connected.
                        int w = u;
                        for (int attempt = 0; attempt < 100; attempt++) {
                            w = gen() % n;
                            if (w == u) continue;
                            // Check if w is already a neighbor of u.
                            bool alreadyNeighbor = false;
                            for (int nbr : adj[u]) {
                                if (nbr == w) {
                                    alreadyNeighbor = true;
                                    break;
                                }
                            }
                            if (!alreadyNeighbor) break;
                        }
                        // Remove the existing edge (u, v) from both u and v.
                        auto it = remove(adj[u].begin(), adj[u].end(), v);
                        adj[u].erase(it, adj[u].end());
                        auto it2 = remove(adj[v].begin(), adj[v].end(), u);
                        adj[v].erase(it2, adj[v].end());
                        
                        // Remove the edge from the edges list (simple scan; can be optimized).
                        for (auto it_e = edges.begin(); it_e != edges.end(); ++it_e) {
                            if ((it_e->u == u && it_e->v == v) || (it_e->u == v && it_e->v == u)) {
                                edges.erase(it_e);
                                break;
                            }
                        }
                        // Add the new edge (u, w).
                        add_edge(u, w);
                    }
                }
            }
        }
    }
    
    // Save the graph to file in edge-list format.
    bool save_to_file(const string &filename) const {
        ofstream file(filename);
        if (!file.is_open()){
            cerr << "Error: Could not open file " << filename << " for writing." << endl;
            return false;
        }
        file << n << " " << edges.size() << "\n";
        for (const auto &e : edges) {
            file << e.u << " " << e.v << "\n";
        }
        file.close();
        return true;
    }
};

int main(int argc, char** argv) {
    if(argc < 5) {
        cerr << "Usage: " << argv[0] << " <num_nodes> <k> <beta> <output_file>\n";
        return 1;
    }
    
    int n = stoi(argv[1]);  // total nodes
    int k = stoi(argv[2]);  // each node connected to k nearest neighbors
    double beta = stod(argv[3]);  // rewiring probability [0,1]
    string output_file = argv[4];
    
    if(n <= 0 || k <= 0 || beta < 0.0 || beta > 1.0) {
        cerr << "Error: Invalid parameters." << endl;
        return 1;
    }
    
    Graph g(n);
    g.generate_watts_strogatz(k, beta);
    cout << "Generated Watts–Strogatz graph with " << n << " nodes." << endl;
    
    if(g.save_to_file(output_file)) {
        cout << "Graph saved to " << output_file << endl;
    } else {
        cerr << "Failed to save graph." << endl;
        return 1;
    }
    
    return 0;
}

// compile: g++ -std=c++17 -O2 GraphSim_watts.cpp -o watts_graph_generator
// then run:
//          Watts–Strogatz graph with 1000 nodes, where each node is initially 
//          connected to 10 nearest neighbors (make sure k is even), and each edge is rewired with probability 0.1.
// ./watts_graph_generator 1000 10 0.1 watts_graph_1000_10_01.txt
