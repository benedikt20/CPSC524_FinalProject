// Barabási–Albert Graph Generator (like social networks)
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <unordered_set>
#include <algorithm>

using namespace std;

struct Edge {
    int u, v;
    Edge(int uu, int vv) : u(uu), v(vv) {}
};

class Graph {
public:
    int n;  // total number of nodes
    vector<vector<int>> adj;  // adjacency list
    vector<Edge> edges;       // list of edges

    Graph(int n) : n(n) {
        adj.resize(n);
    }
    
    // Add an undirected edge (avoiding self-loops)
    void add_edge(int u, int v) {
        if(u == v) return;
        // Optionally, you might check for duplicates.
        adj[u].push_back(v);
        adj[v].push_back(u);
        edges.push_back(Edge(u, v));
    }
    
    // Generate Barabási–Albert graph.
    // m0: number of nodes in the seed complete graph.
    // m: number of edges to attach from each new node.
    void generate_barabasi_albert(int m0, int m) {
        if(m0 < m) {
            cerr << "m0 must be at least m." << endl;
            return;
        }
        // Create a complete graph on m0 nodes.
        for (int i = 0; i < m0; i++) {
            for (int j = i + 1; j < m0; j++) {
                add_edge(i, j);
            }
        }
        
        // Prepare a list where each node appears as many times as its current degree.
        vector<int> nodeList;
        // For complete graph, degree for each node is (m0 - 1).
        for (int i = 0; i < m0; i++) {
            for (int j = 0; j < m0 - 1; j++) {
                nodeList.push_back(i);
            }
        }
        
        // Random number generator.
        mt19937 gen(chrono::system_clock::now().time_since_epoch().count());
        
        // For each new node from m0 to n-1, attach m edges.
        for (int new_node = m0; new_node < n; new_node++) {
            unordered_set<int> targets;
            while (targets.size() < static_cast<size_t>(m) && !nodeList.empty()) {
                int idx = gen() % nodeList.size();
                int candidate = nodeList[idx];
                targets.insert(candidate);
            }
            // Attach new_node to each selected target.
            for (int t : targets) {
                add_edge(new_node, t);
            }
            // Update nodeList: add new_node m times (its degree will be m).
            for (int i = 0; i < m; i++){
                nodeList.push_back(new_node);
            }
            // Also, for each target, add one extra occurrence to nodeList.
            for (int t : targets) {
                nodeList.push_back(t);
            }
        }
    }
    
    // Save the graph to file in edge-list format.
    // First line: <n> <num_edges>
    // Then, one edge per line: u v
    bool save_to_file(const string &filename) const {
        ofstream file(filename);
        if (!file.is_open()) {
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
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <num_nodes> <m0> <m> <output_file>" << endl;
        return 1;
    }
    
    int n = stoi(argv[1]);      // total number of nodes
    int m0 = stoi(argv[2]);     // initial complete graph nodes
    int m = stoi(argv[3]);      // edges per new node
    string output_file = argv[4];
    
    if(n <= m0) {
        cerr << "Error: Total nodes must be greater than m0." << endl;
        return 1;
    }
    if(m0 < m) {
        cerr << "Error: m0 must be at least m." << endl;
        return 1;
    }
    
    Graph g(n);
    g.generate_barabasi_albert(m0, m);
    cout << "Generated Barabási–Albert graph with " << n << " nodes." << endl;
    if (g.save_to_file(output_file)) {
        cout << "Graph saved to " << output_file << endl;
    } else {
        cerr << "Failed to save graph." << endl;
        return 1;
    }
    
    return 0;
}

// compile: g++ -std=c++17 -O2 GraphSim_barabasi.cpp -o barabasi_graph_generator
// then run:
//          creates a graph with 1000 nodes starting from a seed of 10 nodes (a complete graph) 
//          and attaches 3 edges from each new node
// ./barabasi_graph_generator 1000 10 3 barabasi_graph_1000_10_3.txt
