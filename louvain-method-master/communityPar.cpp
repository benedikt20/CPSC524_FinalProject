// File: communityPar.cpp
// -- parallelisation of Lovains method
//-----------------------------------------------------------------------------
// Community detection
// Original code from: E. Lefebvre, adapted by J.-L. Guillaume (February 2008 version)
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Edits: Parallelized for multi-core processor machines
// Edited by: Benedikt Farag
// Email    : benedikt.farag@yale.edu
// Location : New Haven CT
// Time	    : April 2025
//-----------------------------------------------------------------------------
// To compile  : g++-12 -O2 -o community_par communityPar.cpp res/graph.cpp res/graph_binary.cpp res/main_community.cpp
// To execute  : ./community_par input_file.bin -w weights_file.bin -v -o output_communities.txt
//      - Only to execute a single file. 
//      - Use the pyton script run_louvain.py to execute for the SBM graphs.

#include "res/community.h"
#include <omp.h> // for parallism

using namespace std;

Community::Community(char * filename, char * filename_w, int type, int nbp, double minm) {
  g = Graph(filename, filename_w, type);
  size = g.nb_nodes;

  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;

  n2c.resize(size);
  in.resize(size);
  tot.resize(size);

  // Parallelization p4
  // -----------------------------------------------------
  // Parallelization of the data initialization
  #pragma omp parallel for
  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    tot[i] = g.weighted_degree(i);
    in[i]  = g.nb_selfloops(i);
  }
  // -----------------------------------------------------
  // End of arallelization p4

  nb_pass = nbp;
  min_modularity = minm;
}

Community::Community(Graph gc, int nbp, double minm) {
  g = gc;
  size = g.nb_nodes;

  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;

  n2c.resize(size);
  in.resize(size);
  tot.resize(size);

  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
    tot[i] = g.weighted_degree(i);
  }

  nb_pass = nbp;
  min_modularity = minm;
}

void
Community::init_partition(char * filename) {
  ifstream finput;
  finput.open(filename,fstream::in);

  // read partition
  while (!finput.eof()) {
    unsigned int node, comm;
    finput >> node >> comm;
    
    if (finput) {
      int old_comm = n2c[node];
      neigh_comm(node);

      remove(node, old_comm, neigh_weight[old_comm]);

      unsigned int i=0;
      for ( i=0 ; i<neigh_last ; i++) {
	unsigned int best_comm     = neigh_pos[i];
	float best_nblinks  = neigh_weight[neigh_pos[i]];
	if (best_comm==comm) {
	  insert(node, best_comm, best_nblinks);
	  break;
	}
      }
      if (i==neigh_last)
	insert(node, comm, 0);
    }
  }
  finput.close();
}


void
Community::display() {
  for (int i=0 ; i<size ; i++)
    cerr << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i] ;
  cerr << endl;
}

// Parallelization p2
// -----------------------------------------------------
// Parallelization of the modularity function
double Community::modularity() {
  double q = 0.;
  double m2 = (double)g.total_weight;
  
  #pragma omp parallel
  {
      double q_local = 0.0;
      
      #pragma omp for nowait
      for (int i=0; i<size; i++) {
          if (tot[i]>0)
              q_local += (double)in[i]/m2 - ((double)tot[i]/m2)*((double)tot[i]/m2);
      }
      
      #pragma omp atomic
      q += q_local;
  }
  return q;
}
// Original structure:
//
// double
// Community::modularity() {
//   double q  = 0.;
//   double m2 = (double)g.total_weight;

//   for (int i=0 ; i<size ; i++) {
//     if (tot[i]>0)
//       q += (double)in[i]/m2 - ((double)tot[i]/m2)*((double)tot[i]/m2);
//   }

//   return q;
// }
// -----------------------------------------------------
// End of parallelization 2

void
Community::neigh_comm(unsigned int node) {
  for (unsigned int i=0 ; i<neigh_last ; i++)
    neigh_weight[neigh_pos[i]]=-1;
  neigh_last=0;

  pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(node);

  unsigned int deg = g.nb_neighbors(node);

  neigh_pos[0]=n2c[node];
  neigh_weight[neigh_pos[0]]=0;
  neigh_last=1;

  for (unsigned int i=0 ; i<deg ; i++) {
    unsigned int neigh        = *(p.first+i);
    unsigned int neigh_comm   = n2c[neigh];
    double neigh_w = (g.weights.size()==0)?1.:*(p.second+i);
    
    if (neigh!=node) {
      if (neigh_weight[neigh_comm]==-1) {
	    neigh_weight[neigh_comm]=0.;
	    neigh_pos[neigh_last++]=neigh_comm;
      }
      neigh_weight[neigh_comm]+=neigh_w;
    }
  }
}

void
Community::partition2graph() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;


  for (int i=0 ; i<size ; i++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(i);

    int deg = g.nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
    }
  }
}

void
Community::display_partition() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  for (int i=0 ; i<size ; i++)
    cout << i << " " << renumber[n2c[i]] << endl;
}


Graph
Community::partition2graph_binary() {
  // Renumber communities
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  // Compute communities
  vector<vector<int> > comm_nodes(final);
  vector<vector<int> > communities(final);
  //printf("%s %d \n", __FILE__, __LINE__);
  for (int node=0 ; node<size ; node++) {
    //TODO add node handling
    vector<int>& comm = communities[renumber[n2c[node]]];  
    comm.insert(comm.end(), g.nodes[node].begin(), g.nodes[node].end());
    comm_nodes[renumber[n2c[node]]].push_back(node);
  }
  /*
  for (size_t i=0; i<g.nodes.size(); i++) {
    for (size_t j=0; j<g.nodes[i].size(); j++) {
      printf("%d n ", g.nodes[i][j]);
    }
    printf("n \n");
  }
  */
  /*
  for (size_t i=0; i<communities.size(); i++) {
    for (size_t j=0; j<communities[i].size(); j++) {
      printf("%d c ", communities[i][j]);
    }
    printf("i: %d \n", i);
  }
  printf("%s %d \n", __FILE__, __LINE__);
  */
  Graph g2(communities);  
  
  g2.nb_nodes = comm_nodes.size();
  g2.degrees.resize(comm_nodes.size());

  int comm_deg = comm_nodes.size();

  // Parallelization p3
  // -----------------------------------------------------
  // Parallelization of the new graph construction (with a necessary crititcal section)
  #pragma omp parallel for
  for (int comm=0 ; comm<comm_deg ; comm++) {
    map<int,float> m;
    map<int,float>::iterator it;

    int comm_size = comm_nodes[comm].size();
    for (int node=0 ; node<comm_size ; node++) {
      pair<vector<unsigned int>::iterator, vector<float>::iterator> p = g.neighbors(comm_nodes[comm][node]);
      int deg = g.nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
    int neigh        = *(p.first+i);
    int neigh_comm   = renumber[n2c[neigh]];
    double neigh_weight = (g.weights.size()==0)?1.:*(p.second+i);

    it = m.find(neigh_comm);
    if (it==m.end())
      m.insert(make_pair(neigh_comm, neigh_weight));
    else
      it->second+=neigh_weight;
      }
    }
    g2.degrees[comm]=(comm==0)?m.size():g2.degrees[comm-1]+m.size();

    // global variables need to be protected here
    #pragma omp critical
    {
      g2.nb_links+=m.size();

      for (it = m.begin() ; it!=m.end() ; it++) {
        g2.total_weight  += it->second;
        g2.links.push_back(it->first);
        g2.weights.push_back(it->second);
      }
    }
  }
  // -----------------------------------------------------
  // End of parallelization p3
  

  return g2;
}


bool
Community::one_level() {
  bool improvement=false ;
  int nb_moves;
  int nb_pass_done = 0;
  double new_mod   = modularity();
  double cur_mod   = new_mod;

  vector<int> random_order(size);
  for (int i=0 ; i<size ; i++)
    random_order[i]=i;
  for (int i=0 ; i<size-1 ; i++) {
    int rand_pos = rand()%(size-i)+i;
    int tmp      = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // repeat while 
  //   there is an improvement of modularity
  //   or there is an improvement of modularity greater than a given epsilon 
  //   or a predefined number of pass have been done
  do {
    cur_mod = new_mod;
    nb_moves = 0;
    nb_pass_done++;


    // for each node: remove the node from its community and insert it in the best community

    // Parallelization p0 and p1
    // -----------------------------------------------------
    // Parallelization of the node removal and insertion to antoher community
    // Parallized without (p0) and with dynamic scheduling (p1)
    #pragma omp parallel for schedule(dynamic)
    for (int node_tmp=0 ; node_tmp<size ; node_tmp++) {
//      int node = node_tmp;
      int node = random_order[node_tmp];
      int node_comm     = n2c[node];
      double w_degree = g.weighted_degree(node);

      // computation of all neighboring communities of current node
      neigh_comm(node);
      // remove node from its current community
      remove(node, node_comm, neigh_weight[node_comm]);

      // compute the nearest community for node
      // default choice for future insertion is the former community
      int best_comm        = node_comm;
      double best_nblinks  = 0.;
      double best_increase = 0.;
      for (unsigned int i=0 ; i<neigh_last ; i++) {
        double increase = modularity_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
        if (increase>best_increase) {
          best_comm     = neigh_pos[i];
          best_nblinks  = neigh_weight[neigh_pos[i]];
          best_increase = increase;
        }
      }

      // insert node in the nearest community
      insert(node, best_comm, best_nblinks);
     
      if (best_comm!=node_comm)
        nb_moves++;
    }

    double total_tot=0;
    double total_in=0;
    for (unsigned int i=0 ; i<tot.size() ;i++) {
      total_tot+=tot[i];
      total_in+=in[i];
    }

    new_mod = modularity();
    if (nb_moves>0)
      improvement=true;
    
  } while (nb_moves>0 && new_mod-cur_mod>min_modularity);

  return improvement;
}

