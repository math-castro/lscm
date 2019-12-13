#include "segmentation/segmentation.hpp"

#include <vector>

#include "halfedge/HalfedgeBuilder.hpp"

using namespace std;
using namespace Eigen;

Segmentation::Segmentation(const MatrixXd &V, const MatrixXi &F, igl::opengl::glfw::Viewer &viewer) : V(V), F(F), heds(0,0), viewer(viewer) {
  HalfedgeBuilder hb;
  heds = hb.createMesh(V.rows(), F);
  findInitialFeatures();
  threshold = top_features.begin()->first;
}

void Segmentation::findInitialFeatures() {
  cerr << "Finding Initial features..." << flush;
  top_features.clear();
  for(int i = 0; i < heds.sizeOfHalfedges(); i++) {
    top_features.insert({sharpness(i), i});
    if(top_features.size() > .05*heds.sizeOfHalfedges()) top_features.erase(top_features.begin());
  }
  cerr << "done" << endl;
}

void Segmentation::expandFeatureCurve(int start) {
  vector<int> detected_feature;

  // forwards
  int h = start;
  int hp = h;
  dfsResult best;
  do {
    best = dfs(hp, 1, 0, -1);
    hp = best.second;
    detected_feature.push_back(hp);
  } while (best.sharpness > max_string_length * threshold);

  // backwards
  h = heds.getOpposite(start);
  hp = h;
  do {
    best = dfs(hp, 1, 0, -1);
    hp = best.second;
    detected_feature.push_back(hp);
  } while (best.sharpness > max_string_length * threshold);

  if (detected_feature.size() > min_feature_length) {
    tagFeatures(detected_feature);
  }
}

void Segmentation::tagFeatures(vector<int> &features) {
  for(int f : features) {
    tagAsFeature(f);
    tagNeighborhood(f);
  }
}

void Segmentation::tagAsFeature(int h) {
  tag[h] = tag[heds.getOpposite(h)] = 1;
}

void Segmentation::tagNeighborhood(int h) {
  // tag neighboors as neighboors
  vector<int> neighbors = outNeighbors(h);
  for (int n : neighbors) {
    tagAsNeighbor(n);
  }
}

void Segmentation::tagAsNeighbor(int h) {
  if(!tag[h]) tag[h] = 2;
  h = heds.getOpposite(h);
  if(!tag[h]) tag[h] = 2;
}

dfsResult Segmentation::dfs(int h, int length, double sharp, int second) {
  if(length == 2) second = h;
  sharp += sharpness(h);

  if(length == max_string_length) {
    return {sharp, second};
  }

  vector<int> out = outNeighbors(h);
  dfsResult best{-1,-1};

  for(int nh : out) {
    if(nh == heds.getOpposite(h) || tag[nh])
      continue;
    dfsResult temp = dfs(nh, length+1, sharp, second);
    if(temp.sharpness > best.sharpness) {
      temp = best;
    }
  }

  return best;
}

double Segmentation::sharpness(int h) {
  int f1 = heds.getFace(h);
  int f2 = heds.getFace(heds.getOpposite(h));
  if(f1 == -1 or f2 == -1)
    return M_PI;
  Vector3d n1 = normal(f1);
  Vector3d n2 = normal(f2);
  return acos(n1.dot(n2));
}

Vector3d Segmentation::normal(int f) {
  Vector3d v1 = V.row(F(f, 1)) - V.row(F(f, 0));
  Vector3d v2 = V.row(F(f, 2)) - V.row(F(f, 0));
  return v1.cross(v2).normalized();
}

void Segmentation::colorInitialFeatures() {
  cerr << "Coloring Initial features..." << flush;
  for(auto &p : top_features)
    colorEdge(p.second);
  cerr << "done" << endl;
}

void Segmentation::colorEdge(int h) {
  RowVector3d v1 = V.row(heds.getTarget(h)); 
  RowVector3d v2 = V.row(heds.getTarget(heds.getOpposite(h)));
  viewer.data().add_edges(v1, v2, RowVector3d(1,0,0));
}

vector<int> Segmentation::outNeighbors(int h) {
  return vector<int>();
}

vector<int> Segmentation::inNeighbors(int h) {
  return vector<int>();
}