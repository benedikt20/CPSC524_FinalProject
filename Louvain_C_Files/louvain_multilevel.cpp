// Created by Elisa Sommer on 15.04.25.
// louvain_multilevel.cpp
// This program implements the multi-level (hierarchical) Louvain method for community detection.
// It reads a graph from a file, performs iterative local optimizations and graph aggregations,
// and finally computes and prints the modularity of the partition.

// -----------------------------
// Includes and Global Debug Flag
// -----------------------------
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <queue>
#include <cstdlib>

using namespace std;

// Global debug flag: set to true to enable debug printouts.
const bool DEBUG = false;

// --------------------------------------
// Graph Structure and Helper Functions
// --------------------------------------

// Structure representing an undirected graph using an adjacency list.
struct Graph {
    int n;                      // Number of vertices in the graph.
    vector<vector<int>> adj;    // Adjacency list: for each vertex, a list of its neighbors.
};

// Reads a graph from a file.
// The expected format:
//   First line: <number_of_vertices> <number_of_edges>
//   Then each subsequent line: <u> <v>  (edge between vertices u and v)
// Vertices are assumed to be 0-indexed.
Graph readGraph(const string &filename) {
    ifstream infile(filename);
    if (!infile.is_open()){
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }
    
    int n, m;
    infile >> n >> m;
    if (DEBUG) {
        cout << "[DEBUG] readGraph: n = " << n << ", m = " << m << endl;
    }
    
    Graph g;
    g.n = n;
    g.adj.resize(n);
    
    int u, v;
    // Read each edge.
    for (int i = 0; i < m; i++){
        infile >> u >> v;
        // Optionally print the first few edges for verification.
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

// Computes the degree (i.e., number of neighbors) for each vertex.
vector<int> computeDegrees(const Graph &g) {
    vector<int> degree(g.n, 0);
    for (int i = 0; i < g.n; i++){
        degree[i] = g.adj[i].size();
    }
    return degree;
}

// --------------------------------------
// Louvain Method Functions
// --------------------------------------

// Performs one level of local optimization (one pass of the Louvain method).
// For each node, it evaluates whether moving the node to one of its neighbors'
// communities improves modularity (using a simplified gain formula).
// Parameters:
//   g: The input graph.
//   initComm: Initial community assignment (typically each node in its own community).
//   changed: Output flag set to true if at least one node changes community in this pass.
// Returns: Updated community assignment after local optimization.
vector<int> localLouvain(const Graph &g, const vector<int> &initComm, bool &changed) {
    int n = g.n;
    vector<int> comm = initComm;          // Current community assignment.
    vector<int> degree = computeDegrees(g); // Degrees for all nodes.
    
    // Compute total number of edges: m2 = (sum of degrees) / 2.
    double m2 = 0;
    for (int d : degree) {
        m2 += d;
    }
    m2 /= 2.0;
    
    changed = false;
    bool localChange = true;
    int iteration = 0;  // To track number of iterations in local optimization.
    
    // Continue until no move produces improvement.
    while (localChange) {
        iteration++;
        if (DEBUG) {
            cout << "[DEBUG] localLouvain - Iteration " << iteration << " starting." << endl;
        }
        localChange = false;
        
        // Process every node.
        for (int i = 0; i < n; i++){
            int currentComm = comm[i];  // Current community of node i.
            int bestComm = currentComm; // Best candidate community for node i.
            double bestGain = 0.0;      // Best modularity gain found.
            
            // Count links from node i to each community among its neighbors.
            unordered_map<int, int> neighborCommCount;
            for (int j : g.adj[i]) {
                neighborCommCount[comm[j]]++;
            }
            
            // Evaluate potential move: try moving node i to each neighbor's community.
            for (auto &entry : neighborCommCount) {
                int targetComm = entry.first;
                if (targetComm == currentComm)
                    continue;
                
                int k_i_in = entry.second;  // Count of edges from i to target community.
                int sumTot = 0;
                // Sum degrees of nodes in target community.
                for (int j = 0; j < n; j++){
                    if (comm[j] == targetComm)
                        sumTot += degree[j];
                }
                
                // Simplified gain formula:
                // ΔQ ≈ [k_i_in - (degree[i]*sumTot)/(2*m2)] / m2
                double deltaQ = (static_cast<double>(k_i_in) - 
                                 (static_cast<double>(degree[i]) * sumTot) / (2.0 * m2)) / m2;
                if (deltaQ > bestGain) {
                    bestGain = deltaQ;
                    bestComm = targetComm;
                }
            }
            
            // If a beneficial move is found, update the community assignment.
            if (bestComm != currentComm) {
                if (DEBUG) {
                    cout << "[DEBUG] localLouvain: Node " << i << " moves from community " 
                         << currentComm << " to " << bestComm 
                         << " (gain: " << bestGain << ")" << endl;
                }
                comm[i] = bestComm;
                localChange = true;
                changed = true;
            }
        }
        if (DEBUG) {
            cout << "[DEBUG] localLouvain - Iteration " << iteration << " completed." << endl;
        }
        // (Optional) Safety check: break if too many iterations occur.
        if (iteration > 10 * g.n) {
            cout << "[DEBUG] localLouvain: Exceeded 10 * N iterations; breaking to avoid infinite loop." << endl;
            break;
        }
    }
    return comm;
}

// Aggregates the graph according to the current community assignment.
// Every community becomes a single node in the aggregated graph.
// Parameters:
//   g: The original graph.
//   comm: Current community assignment.
//   new_n: (Output) number of communities in the aggregated graph.
//   commMap: (Output) mapping from original nodes to aggregated community IDs.
// Returns: A new aggregated graph.
Graph aggregateGraph(const Graph &g, const vector<int> &comm, int &new_n, vector<int> &commMap) {
    unordered_map<int, int> remap; // Map old community labels to new contiguous labels.
    new_n = 0;
    
    // Build remapping.
    for (int i = 0; i < g.n; i++){
        int c = comm[i];
        if (remap.find(c) == remap.end()){
            remap[c] = new_n;
            new_n++;
        }
    }
    // Build mapping from original nodes to new community ID.
    commMap.resize(g.n, -1);
    for (int i = 0; i < g.n; i++){
        commMap[i] = remap[comm[i]];
    }
    
    // Create new aggregated graph.
    Graph newGraph;
    newGraph.n = new_n;
    newGraph.adj.resize(new_n);
    
    // Use a hash key to avoid duplicate edges.
    unordered_map<long long, bool> added;
    for (int i = 0; i < g.n; i++){
        for (int j : g.adj[i]){
            int ci = remap[comm[i]];
            int cj = remap[comm[j]];
            if (ci == cj) continue;  // Skip self-loop.
            long long key = ((long long)min(ci, cj) << 32) | max(ci, cj);
            if (added.find(key) == added.end()){
                newGraph.adj[ci].push_back(cj);
                newGraph.adj[cj].push_back(ci);
                added[key] = true;
            }
        }
    }
    if (DEBUG) {
        cout << "[DEBUG] aggregateGraph: Aggregated graph has " << new_n << " nodes." << endl;
    }
    return newGraph;
}

// Implements the full multilevel Louvain method.
// The algorithm repeatedly performs local optimization (using localLouvain)
// and then aggregates the graph until no further improvement is possible.
// Returns the final community assignment for the original graph.
vector<int> multiLevelLouvain(const Graph &g) {
    int n = g.n;
    // Initially, assign each node to its own community.
    vector<int> comm(n);
    for (int i = 0; i < n; i++){
        comm[i] = i;
    }
    
    Graph currentGraph = g;                // Graph to be optimized.
    vector<vector<int>> hierarchy;           // Holds community assignments for each level.
    hierarchy.push_back(comm);               // Level 0: initial assignment.
    
    int level = 0;
    bool improvement = true;
    while (improvement) {
        if (DEBUG) {
            cout << "[DEBUG] multiLevelLouvain: Starting optimization at level " << level << endl;
        }
        bool changed = false;
        // Perform local optimization on the current graph.
        vector<int> newComm = localLouvain(currentGraph, hierarchy.back(), changed);
        hierarchy.back() = newComm;  // Update current level.
        
        // If no node moved, we are done.
        if (!changed) break;
        
        int new_n;
        vector<int> commMap;
        // Aggregate the current graph using the new community assignment.
        Graph aggregated = aggregateGraph(currentGraph, newComm, new_n, commMap);
        // If aggregation does not reduce the number of nodes, stop.
        if (new_n == currentGraph.n) break;
        hierarchy.push_back(commMap);  // Save the mapping from current to aggregated level.
        currentGraph = aggregated;       // Work on the aggregated graph next.
        level++;
        if (DEBUG) {
            cout << "[DEBUG] multiLevelLouvain: Completed level " << level << ". New graph size: " << new_n << endl;
        }
        // Safety break: limit number of levels if needed.
        if (level > 10 * g.n) {
            cout << "[DEBUG] multiLevelLouvain: Exceeded 10 * N levels; breaking to avoid infinite loop." << endl;
            break;
        }
    }
    
    // Reconstruct the final community assignment for the original graph.
    vector<int> finalComm(g.n);
    finalComm = hierarchy[0];
    for (size_t lvl = 1; lvl < hierarchy.size(); lvl++){
        vector<int> mapping = hierarchy[lvl];
        for (int i = 0; i < g.n; i++){
            finalComm[i] = mapping[finalComm[i]];
        }
    }
    if (DEBUG) {
        cout << "[DEBUG] multiLevelLouvain: Final community reconstruction completed." << endl;
    }
    return finalComm;
}

// --------------------------------------
// Modularity Computation Function
// --------------------------------------
// Computes the modularity Q for a given community assignment on the original graph.
// Q is defined as:
//   Q = (1/(2*m)) * Σ[i,j] [A(i,j) - (k_i * k_j)/(2*m)] * δ(c_i, c_j)
// where m is the total number of edges and δ(c_i, c_j)=1 if nodes i and j are in the same community.
double compute_modularity(const Graph &g, const vector<int> &comm) {
    int n = g.n;
    vector<int> degree = computeDegrees(g);
    double m = 0;
    for (int d : degree) {
        m += d;
    }
    m /= 2.0;  // Total number of edges.
    
    double Q = 0.0;
    // Sum over all edges (note: each undirected edge is encountered twice).
    for (int i = 0; i < n; i++){
        for (int j : g.adj[i]){
            if (comm[i] == comm[j])
                Q += 1 - (static_cast<double>(degree[i]) * degree[j]) / (2 * m);
        }
    }
    Q /= (2 * m);
    return Q;
}

// --------------------------------------
// Main Function
// --------------------------------------
// Reads the graph from file, applies the multi-level Louvain method,
// computes the modularity, and outputs the final community assignments and performance metrics.
int main(int argc, char* argv[]){
    if (argc < 2){
        cerr << "Usage: " << argv[0] << " <graph_file>" << endl;
        return 1;
    }
    
    string filename = argv[1];
    if (DEBUG) {
        cout << "[DEBUG] Main: Reading graph from file: " << filename << endl;
    }
    Graph g = readGraph(filename);
    if (DEBUG) {
        cout << "[DEBUG] Main: Graph read completed. Total nodes: " << g.n << endl;
    }
    
    // Record start time.
    auto start = chrono::high_resolution_clock::now();
    // Run the multi-level Louvain method to compute community assignments.
    vector<int> communities = multiLevelLouvain(g);
    // Record end time.
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> execTime = end - start;
    
    // Output the final community assignments for each node.
    cout << "Final community assignments for each node:" << endl;
    for (int i = 0; i < static_cast<int>(communities.size()); i++){
        cout << "Node " << i << " -> Community " << communities[i] << endl;
    }
    
    // Compute and output the modularity.
    double mod = compute_modularity(g, communities);
    cout << "Modularity: " << mod << endl;
    cout << "Multi-level Louvain execution time: " << execTime.count() << " seconds." << endl;
    
    return 0;
}

//  compile: g++ -std=c++17 -O2 Louvain_C_Files/louvain_multilevel.cpp -o Louvain_C_Files/louvain_multilevel
//  then run automatically from Python file
// note: debug variable in this code can be set to TRUE to get progress printouts