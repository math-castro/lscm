#pragma once

#include <set>
#include <utility>
#include <queue>

#include "halfedge/HalfedgeDS.hpp"

struct dfsResult {
  double sharpness;
  int second;
};

class Segmentation {
 public:

  Segmentation(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, igl::opengl::glfw::Viewer &viewer);
  void expandFeatureCurves();
  void colorInitialFeatures();
  void colorExpandedFeatures();

 private:

  void findInitialFeatures();
  void expandFeatureCurve(int start);
  double sharpness(int h);
  Eigen::Vector3d normal(int f);
  dfsResult dfs(int u, int depth, double sharp, int second, std::set<int> &detected_feature);
  std::vector<int> outNeighbors(int h);
  std::vector<int> inNeighbors(int h);
  void tagFeatures(std::vector<int> &features);
  void tagAsFeature(int h);
  void tagNeighborhoods(std::queue<std::pair<int,int>> &q, int depth);
  void tagAsNeighbor(int h);
  void colorEdge(int h);

  HalfedgeDS heds;
  const Eigen::MatrixXd &V;
  const Eigen::MatrixXi &F;
  igl::opengl::glfw::Viewer &viewer;
  std::vector<int> tag;
  const int max_string_length = 5;
  const int min_feature_length = 15;
  double threshold;
  std::set<std::pair<double,int>> top_features;
};