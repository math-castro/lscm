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
  tag.assign(heds.sizeOfHalfedges(), 0);
}

void Segmentation::expandFeatureCurves() {
  for(auto &p : top_features) {
    expandFeatureCurve(p.second);
  }
}

void Segmentation::colorExpandedFeatures() {
  cerr << "Coloring expanded features..." << flush;
  for(int i = 0; i < heds.sizeOfHalfedges(); i++)
    if(tag[i]==1)
      colorEdge(i);
  cerr << "done" << endl;
}

void Segmentation::colorInitialFeatures() {
  cerr << "Coloring Initial features..." << flush;
  for(auto &p : top_features)
    colorEdge(p.second);
  cerr << "done" << endl;
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
  if(tag[start]) return;

  vector<int> detected_feature;
  set<int> detected_vertices;
  detected_vertices.insert(heds.getTarget(start));
  detected_feature.push_back(start);

  // forwards
  int h = start;
  int hp = h;
  dfsResult best;
  do {
    best = dfs(hp, 1, 0, -1, detected_vertices);
    hp = best.second;
    if(hp != -1) {
      detected_vertices.insert(heds.getTarget(hp));
      detected_feature.push_back(hp);
    }
    else break;
  } while (best.sharpness > max_string_length * threshold);

  // backwards
  h = heds.getOpposite(start);
  hp = h;
  do {
    best = dfs(hp, 1, 0, -1, detected_vertices);
    hp = best.second;
    if(hp != -1) {
      detected_vertices.insert(heds.getTarget(hp));
      detected_feature.push_back(hp);
    }
    else break;
  } while (best.sharpness > max_string_length * threshold);


  if (detected_feature.size() > min_feature_length) {
    tagFeatures(detected_feature);
  }
}

dfsResult Segmentation::dfs(int h, int length, double sharp, int second, set<int> &detected_feature) {
  int prev_tag = tag[h];
  tag[h] = 3;

  if(length == 2) second = h;
  sharp += sharpness(h);

  if(length == max_string_length) {
    tag[h] = prev_tag;
    return {sharp, second};
  }

  vector<int> out = outNeighbors(h);

  dfsResult best{1,-1};

  for(int nh : out) {
    if(nh == heds.getOpposite(h) || tag[nh] || detected_feature.count(heds.getTarget(nh)))
      continue;
    dfsResult temp = dfs(nh, length+1, sharp, second, detected_feature);
    if(temp.sharpness > best.sharpness) {
      best = temp;
    }
  }

  tag[h] = prev_tag;
  return best;
}

void Segmentation::tagFeatures(vector<int> &features) {
  queue<pair<int,int>> q;
  for(int f : features) {
    tagAsFeature(f);
    q.push({f,0});
  }
  tagNeighborhoods(q,10);
}

void Segmentation::tagAsFeature(int h) {
  tag[h] = tag[heds.getOpposite(h)] = 1;
}

void Segmentation::tagNeighborhoods(queue<pair<int,int>> &q, int depth) {
  // tag neighboors as neighboors
  while(q.size()) {
    auto p = q.front(); q.pop();
    int h = p.first, d = p.second;

    if(d == depth) continue;

    vector<int> neighbors = outNeighbors(h);
    for (int n : neighbors) {
      if(!tag[n]) {
        tagAsNeighbor(n);
        q.push({n,d+1});
      }
    }
  }
}

void Segmentation::tagAsNeighbor(int h) {
  if(!tag[h]) tag[h] = 2;
  h = heds.getOpposite(h);
  if(!tag[h]) tag[h] = 2;
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

void Segmentation::colorEdge(int h) {
  RowVector3d v1 = V.row(heds.getTarget(h)); 
  RowVector3d v2 = V.row(heds.getTarget(heds.getOpposite(h)));
  viewer.data().add_edges(v1, v2, RowVector3d(1,0,0));
}

vector<int> Segmentation::outNeighbors(int h) {
  vector<int> neighbours;
  h = heds.getOpposite(h);
  int hn = h;
  do {
    neighbours.push_back(hn);
    hn = heds.getNext(heds.getOpposite(hn));
  } while(hn != h);
  return neighbours;
}

vector<int> Segmentation::inNeighbors(int h) {
  return vector<int>();
}