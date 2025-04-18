// File: graph_binary.cpp
// -- graph handling source
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void write_graph_to_binary(const string &filename) {
    // Define the graph parameters
    int nb_nodes = 4;
    vector<int> degrees = {2, 2, 2, 1}; // Cumulative degrees for nodes (for simplicity)
    vector<int> links = {1, 0, 2, 1, 3}; // The nodes each node is connected to
    vector<float> weights = {2.0f, 1.5f, 3.0f}; // Weights for the edges

    // Open a binary file to write
    ofstream outfile(filename, ios::binary);

    if (!outfile) {
        cerr << "Error opening file for writing." << endl;
        return;
    }

    // Write the number of nodes (4 bytes)
    outfile.write(reinterpret_cast<const char*>(&nb_nodes), sizeof(nb_nodes));

    // Write the cumulative degrees (8 * nb_nodes bytes)
    for (int i = 0; i < nb_nodes; i++) {
        outfile.write(reinterpret_cast<const char*>(&degrees[i]), sizeof(degrees[i]));
    }

    // Write the links (4 * sum_degrees bytes)
    for (int i = 0; i < links.size(); i++) {
        outfile.write(reinterpret_cast<const char*>(&links[i]), sizeof(links[i]));
    }

    // Write the weights (4 * sum_degrees bytes, for weighted graph)
    for (int i = 0; i < weights.size(); i++) {
        outfile.write(reinterpret_cast<const char*>(&weights[i]), sizeof(weights[i]));
    }

    // Close the file
    outfile.close();

    cout << "Graph written to binary file." << endl;
}

int main() {
    write_graph_to_binary("input_file.bin");
    return 0;
}
