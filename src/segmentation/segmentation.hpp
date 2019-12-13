#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include <vector>
#include <utility>

#include "halfedge/HalfedgeDS.hpp"

class Segmentation {
 public:
  Segmentation(Eigen::MatrixXd &V_original, Eigen::MatrixXi &F_original, HalfedgeDS &mesh,
               igl::opengl::glfw::Viewer &viewer_);
  void colorSharpEdges();
  void colorFeatures();
  Eigen::Vector3d getNormal(int f);
  float getEdgeSharpness(int e);
  void getEdgeSharpnessMatrix();
  void setThreshold(float newThreshold);
  float getThreshold();
  std::vector<int> getInNeighbours(int v);
  std::vector<int> getOutNeighbours(int v);
  std::pair<float, std::vector<int>> DFS(int current, std::vector<int> &S,
                                         float sharpness, int length,
                                         int maxStringLength);
  void tagNeighbours(int v);
  void expandFeatureCurve(int startEdge);
  void print(std::vector<int> v);

 private:
  int debug;
  igl::opengl::glfw::Viewer *viewer;
  float threshold;
  int *EdgeMap;
  std::vector<int> sharpEdges;
  int *tag;
  /** Half-edge representation of the original input mesh */
  HalfedgeDS *he;
  Eigen::MatrixXd *V;  // vertex coordinates of the original input mesh
  Eigen::MatrixXi *F;  // REMARK: not needed if using the half-edge data structure
  Eigen::VectorXd *EdgeSharpness;

  int nEdges;
  // int nVertices, nFaces; // number of vertices, faces
};
