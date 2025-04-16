// Created by Elisa Sommer on 15.04.25.
// louvain_sequential.cpp
// This program implements a simple one-level Louvain method for community detection.
// It reads a graph from a file, performs local optimization of community assignments,
// computes the modularity of the resulting partition, and prints the results.

// -----------------------------
// Includes and Global Debug Flag
// -----------------------------
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>
#include <cstdlib>
#include <algorithm>

using namespace std;

// Global debug flag: set to true to enable debug printouts.
const bool DEBUG = false;

// -----------------------------
// Graph Structure and Functions
// -----------------------------

// Structure representing an undirected graph using an adjacency list.
struct Graph {
    int n;                      // Number of vertices.
    vector<vector<int>> adj;    // Adjacency list: each vertex's list of neighbors.
};

// Reads a graph from a file.
// Expected file format:
//   First line: <number_of_vertices> <number_of_edges>
//   Each subsequent line: <u> <v> (edge between vertices u and v)
// Assumes vertices are 0-indexed.
Graph readGraph(const string &filename) {
    ifstream infile(filename);
    if (!infile.is_open()){
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }
    
    int n, m;
    infile >> n >> m;
    
    // Debug: Output graph header info.
    if (DEBUG) {
        cout << "[DEBUG] readGraph: n = " << n << ", m = " << m << endl;
    }
    
    Graph g;
    g.n = n;
    g.adj.resize(n);
    
    int u, v;
    for (int i = 0; i < m; i++){
        if (!(infile >> u >> v)) {
            cerr << "Error reading edge " << i << endl;
            break;
        }
        // Debug: Print the first few edges.
        if (DEBUG && i < 5) {
            cout << "[DEBUG] Edge " << i << ": " << u << " - " << v << endl;
        }
        // Add edge to both nodes (undirected graph).
        g.adj[u].push_back(v);
        g.adj[v].push_back(u);
    }
    
    infile.close();
    return g;
}

// Computes the degree (number of neighbors) for each vertex.
vector<int> computeDegrees(const Graph &g) {
    vector<int> degree(g.n, 0);
    for (int i = 0; i < g.n; i++){
        degree[i] = g.adj[i].size();
    }
    return degree;
}

// -----------------------------
// Simple One-level Louvain Method (with Debug Printouts)
// -----------------------------
// In this implementation, each node starts in its own community.
// The algorithm iterates over each node and moves it to one of its neighbor's communities
// if such a move yields a positive (simplified) modularity gain.
// The iteration stops when no single node move improves the measure.
vector<int> simpleLouvain(const Graph &g, const vector<int> &initComm, bool &changed) {
    int n = g.n;
    vector<int> comm = initComm;         // Current community assignment.
    vector<int> degree = computeDegrees(g); // Degrees for all nodes.
    
    // Compute total number of edges: m2 = (sum of degrees) / 2.
    double m2 = 0;
    for (int d : degree)
        m2 += d;
    m2 /= 2.0;
    
    changed = false;
    bool localChange = true;
    int iteration = 0;
    
    // Loop until no node is moved.
    while (localChange) {
        iteration++;
        int moves = 0;
        if (DEBUG) {
            cout << "[DEBUG] simpleLouvain: Iteration " << iteration << " starting." << endl;
        }
        localChange = false;
        
        // Process each node.
        for (int i = 0; i < n; i++){
            int currentComm = comm[i];  // Current community of node i.
            int bestComm = currentComm; // Best candidate community.
            double bestGain = 0.0;      // Best gain found.
            
            // Count number of edges from node i to each community present in its neighbors.
            unordered_map<int, int> neighborCommCount;
            for (int neighbor : g.adj[i]) {
                neighborCommCount[comm[neighbor]]++;
            }
            
            // Evaluate potential moves: try moving node i to each neighboring community.
            for (auto &entry : neighborCommCount) {
                int targetComm = entry.first;
                if (targetComm == currentComm)
                    continue;
                int k_i_in = entry.second;  // Edges from i to nodes in target community.
                int sumTot = 0;
                // Calculate total degree of nodes in target community.
                for (int j = 0; j < n; j++){
                    if (comm[j] == targetComm)
                        sumTot += degree[j];
                }
                
                // Simplified modularity gain calculation:
                // ΔQ ≈ [k_i_in - (degree[i] * sumTot) / (2*m2)] / m2
                double deltaQ = (static_cast<double>(k_i_in) - 
                                 (static_cast<double>(degree[i]) * sumTot) / (2.0 * m2)) / m2;
                if (deltaQ > bestGain) {
                    bestGain = deltaQ;
                    bestComm = targetComm;
                }
            }
            
            // If moving node i yields improvement, update its community.
            if (bestComm != currentComm) {
                if (DEBUG) {
                    cout << "[DEBUG] simpleLouvain: Node " << i << " moves from community " 
                         << currentComm << " to " << bestComm 
                         << " (gain: " << bestGain << ")" << endl;
                }
                comm[i] = bestComm;
                moves++;
                localChange = true;
                changed = true;
            }
        }
        if (DEBUG) {
            cout << "[DEBUG] simpleLouvain: Iteration " << iteration << " completed with " 
                 << moves << " moves." << endl;
        }
        // Safety check to avoid infinite loops.
        if (iteration > 10 * g.n) {
            cout << "[DEBUG] simpleLouvain: Exceeded 10 * N iterations; breaking out." << endl;
            break;
        }
    }
    return comm;
}

// -----------------------------
// Modularity Computation
// -----------------------------
// Computes the modularity Q for the given community assignment on the original graph.
// The formula is: Q = (1/(2m)) * Σ[i,j] [A(i,j) - (k_i*k_j)/(2m)] * δ(c_i,c_j),
// where m is the total number of edges and δ is 1 if nodes i and j are in the same community.
double computeModularity(const Graph &g, const vector<int>& comm) {
    int n = g.n;
    vector<int> degree = computeDegrees(g);
    double m = 0;
    for (int d : degree) {
        m += d;
    }
    m /= 2.0;
    
    double Q = 0.0;
    // Sum over all edges (each edge appears twice in the loops).
    for (int i = 0; i < n; i++){
        for (int j : g.adj[i]){
            if (comm[i] == comm[j])
                Q += 1 - (static_cast<double>(degree[i]) * degree[j]) / (2.0 * m);
        }
    }
    Q /= (2.0 * m);
    return Q;
}

// -----------------------------
// Main Function
// -----------------------------
int main(int argc, char* argv[]){
    // Disable buffering for real-time output.
    setvbuf(stdout, NULL, _IONBF, 0);
    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <graph_file>" << endl;
        return 1;
    }
    
    string filename = argv[1];
    if (DEBUG) {
        cout << "[DEBUG] Main: Reading graph from file: " << filename << endl;
    }
    Graph g = readGraph(filename);
    
    if (DEBUG) {
        cout << "[DEBUG] Main: Graph reading completed. Total nodes: " << g.n << endl;
    }
    
    // Initialize community assignment: each node in its own community.
    vector<int> initComm(g.n);
    for (int i = 0; i < g.n; i++){
        initComm[i] = i;
    }
    
    // Run the simple one-level Louvain algorithm.
    auto start = chrono::high_resolution_clock::now();
    bool changed;
    vector<int> communities = simpleLouvain(g, initComm, changed);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> execTime = end - start;
    
    // Output final community assignments.
    cout << "Final community assignments for each node:" << endl;
    for (int i = 0; i < g.n; i++){
        cout << "Node " << i << " -> Community " << communities[i] << endl;
    }
    
    // Compute and output the modularity.
    double modularity = computeModularity(g, communities);
    cout << "Modularity: " << modularity << endl;
    cout << "Simple one-level Louvain execution time: " << execTime.count() << " seconds." << endl;
    
    return 0;
}

// compile: g++ -std=c++17 -O2 Louvain_C_Files/louvain_sequential.cpp -o Louvain_C_Files/louvain_sequential
// then run automatically from Python file
// note: debug variable in this code can be set to TRUE to get progress printouts