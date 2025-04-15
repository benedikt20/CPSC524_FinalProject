#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

// Include ParlayLib (adjust the path if needed)
#include "parlaylib/include/parlay/primitives.h"
#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/sequence.h"
#include "parlaylib/include/parlay/utilities.h"


// A simple 2D point structure
struct Point2D {
  double x, y;
  Point2D(double xx=0.0, double yy=0.0) : x(xx), y(yy) {}
};

// A helper to compute squared distance
inline double squared_distance(const Point2D& a, const Point2D& b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  return dx*dx + dy*dy;
}

// KD-Tree node
struct KDNode {
  int axis;          // 0 for x, 1 for y
  double splitValue; // coordinate pivot
  int pointIndex;    // index in the original array

  KDNode* left;
  KDNode* right;

  KDNode() : axis(0), splitValue(0.0), pointIndex(-1),
             left(nullptr), right(nullptr) {}
};

// DistIndex for storing (distance^2, index)
struct DistIndex {
  double dist;
  int index;
  DistIndex(double d=0, int i=0) : dist(d), index(i) {}
};

// For a max-heap, we want to put the largest distance on top
inline bool operator<(const DistIndex &a, const DistIndex &b) {
  return a.dist < b.dist;
}

// TODO: Implement a function to build the kd-tree
KDNode* build_kd_tree(
    parlay::slice<int*, int*> indices,
    const parlay::sequence<Point2D>& points,
    int depth = 0
) {
  // 1) Base cases: if 0 or 1 points, create a leaf node or return
  if (indices.size() == 0) return nullptr;
  if (indices.size() == 1) {
    KDNode* node = new KDNode();
    node->pointIndex = indices[0];
    return node;
  }

  // 2) Determine axis = (depth % 2)
  int axis = depth % 2;

  // 3) Sort indices by that axis (x or y)
  parlay::sort_inplace(indices, [&](int a, int b) {
    return axis==0 ? (points[a].x < points[b].x) : (points[a].y < points[b].y);
  });

  // 4) Find median index
  size_t mid_ind = indices.size() / 2;

  // 5) Create a node with that pivot
  KDNode* node = new KDNode();
    node->axis = axis;
    node->pointIndex = indices[mid_ind];
    node->splitValue = (axis == 0) ? 
        points[indices[mid_ind]].x : points[indices[mid_ind]].y;

  // 6) Recurse left and right in parallel
  auto left_ind = indices.cut(0, mid_ind);
  auto right_ind = indices.cut(mid_ind + 1, indices.size());

  parlay::par_do(
    [&]() {node->left = build_kd_tree(left_ind, points, depth + 1); },
    [&]() {node->right = build_kd_tree(right_ind, points, depth + 1); }
  );

  return node; 
  // return nullptr; // placeholder
}

// KNN Helper: holds a local max-heap of size k
class KNNHelper {
public:
  KNNHelper(const parlay::sequence<Point2D>& pts, int kk)
    : points(pts), k(kk) {
    best.reserve(k);
  }

  // Perform recursive search
  void search(const KDNode* node, const Point2D& q) {
    // TODO:
    if (node == nullptr) return; // base case 

    // 1) compute dist2 from q to node->pointIndex
    double dist2 = squared_distance(q, points[node->pointIndex]);

    // 2) update_best if needed
    update_best(dist2, node->pointIndex);

    // 3) compare q's coordinate to splitValue
    double delta = (node->axis == 0) ? (q.x - node->splitValue) : (q.y - node->splitValue);
    const KDNode* first = (delta < 0) ? node->left : node->right;
    const KDNode* second = (delta < 0) ? node->right : node->left;

    // 4) search near side, possibly search far side if needed
    search(first, q);

    if (best.size() < k || delta * delta < best.front().dist){
      search(second, q); // only search second if needed
    }
  }

  // Return final results sorted by ascending distance
  parlay::sequence<DistIndex> get_results() const {
    parlay::sequence<DistIndex> result(best.begin(), best.end());
    parlay::sort_inplace(result, [&](auto &a, auto &b){
      return a.dist < b.dist;
    });
    return result;
  }

private:
  const parlay::sequence<Point2D>& points;
  int k;
  std::vector<DistIndex> best; // will be a max-heap

  // If we have fewer than k, push. Otherwise compare with largest so far.
  void update_best(double dist2, int idx) {
    // TODO: use a max-heap for best and update best with dist2 and idx
    if (best.size() < k){
      best.push_back(DistIndex(dist2, idx)); //add elements if count<k
      if (best.size() == k){
        std::make_heap(best.begin(), best.end()); //once count=k: make a vector containing maxHeap
      }
    }
    else if (dist2 < best.front().dist){
      // if distance<max distance in maxHeap, replace it with this point
      std::pop_heap(best.begin(), best.end());
      best.back() = DistIndex(dist2, idx);
      std::push_heap(best.begin(), best.end());
    }
  }
};

// Parallel k-NN for all queries
parlay::sequence<parlay::sequence<DistIndex>>
knn_search_all(const KDNode* root,
               const parlay::sequence<Point2D>& data_points,
               const parlay::sequence<Point2D>& query_points,
               int k) {
  int Q = (int)query_points.size();
  parlay::sequence<parlay::sequence<DistIndex>> results(Q);

  parlay::parallel_for(0, Q, [&](int i){
    KNNHelper helper(data_points, k);
    helper.search(root, query_points[i]);
    results[i] = helper.get_results();
  });

  return results;
}

// A function to load points from a file
parlay::sequence<Point2D> load_points_from_file(const std::string &filename) {
  // TODO: open file, read N, read N lines of x y into a parlay::sequence
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return {};
  }
  
  int N;
  file >> N;
  
  parlay::sequence<Point2D> points(N);
  for (int i = 0; i < N; i++) {
    double x, y;
    file >> x >> y;
    points[i] = Point2D(x, y);
  }
  
  return points;
  //return {};
}

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <data_file> <query_file> <k>\n";
    return 1;
  }

  std::string data_file  = argv[1];
  std::string query_file = argv[2];
  int k = std::stoi(argv[3]);

  auto data_points = load_points_from_file(data_file);
  int N = (int)data_points.size();

  parlay::sequence<int> indices(N);
  parlay::parallel_for(0, N, [&](int i){ indices[i] = i; });
  KDNode* root = build_kd_tree(indices.cut(0, N), data_points, 0);

  auto query_points = load_points_from_file(query_file);
  int Q = (int)query_points.size();

  auto results = knn_search_all(root, data_points, query_points, k);

  // for (int q = 0; q < std::min(Q, 3); q++) { // initial setup
  for (int q = 0; q < Q; q++) {                 // loop through each query (to print all queries for autograder)
    std::cout << "Query " << q << " : ("
              << query_points[q].x << ", "
              << query_points[q].y << ")\n";
    std::cout << "  kNN: ";
    for (auto &di : results[q]) {
      std::cout << "(dist2=" << di.dist
                << ", idx=" << di.index << ") ";
    }
    std::cout << "\n";
  }

  return 0;
}

// parlaylib locations:
// /Users/BFF/parlaylib
// /Users/BFF/Documents/GitHub/CPSC424/homework5/parlaylib

// In current location:
// compile
// g++-12 -std=c++17 -O2 -fopenmp template.cpp -I /Users/BFF/parlaylib -o kdknn 

// test program with 3 queries
// ./kdknn data.txt queries.txt 3

