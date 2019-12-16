#include "segmentation/segmentation.hpp"

#include <vector>

#include "halfedge/HalfedgeBuilder.hpp"

using namespace std;
using namespace Eigen;

Segmentation::Segmentation(const MatrixXd &V, const MatrixXi &F,
                           igl::opengl::glfw::Viewer &viewer)
    : V(V), F(F), heds(0, 0), viewer(viewer) {
  HalfedgeBuilder hb;
  heds = hb.createMeshWithFaces(V.rows(), F);
}

map<int,vector<int>> Segmentation::run() {
  cerr << "Build incident egges" << endl;
  // buildIncidentEdges();
  cerr << "Find initial features" << endl;
  findInitialFeatures();
  threshold = top_features.begin()->first;
  // colorInitialFeatures();
  cerr << "Expand feature curves" << endl;
  expandFeatureCurves();
  // cerr << "Color expanded feature curves" << endl;
  colorExpandedFeatures();
  cerr << "Calculate distance to features" << endl;
  distanceToFeatures();
  cerr << "Expand charts" << endl;
  expandCharts();
  cerr << "Get charts" << endl;
  return getCharts();
  // return map<int,vector<int>>();
}

void Segmentation::expandCharts() {
  auto comp = [&] (int a, int b) {
    return distanceF[heds.getFace(a)] < distanceF[heds.getFace(b)];
  };
  priority_queue<int, vector<int>, decltype(comp)> pq(comp);

  set<int> chart_boundaries;
  for(int i = 0; i < heds.sizeOfHalfedges(); i++)
    chart_boundaries.insert(i);

  vector<int> max_facets = maximalFacets();
  init_unionfind();

  for (int f : max_facets) {
    par[f] = f;
    // int e = incidentEdge[f];
    int e = heds.getEdgeInFace(f);
    int ne = e;
    do {
      pq.push(ne);
      ne = heds.getNext(ne);
    } while (ne != e);
  }

  while (pq.size()) {
    int e = pq.top();
    pq.pop();
    int f = heds.getFace(e);
    int of = heds.getFace(heds.getOpposite(e));
    if(f == -1 or of == -1) continue;
    if (find(of) == -1) {
      join(f, of);
      chart_boundaries.erase(e);
      // TO-DO: remove non extremal edges (not necessary, will be removed in the course of the algorithm)
      // int oe = incidentEdge[of];
      int oe = heds.getEdgeInFace(of);
      int ne = oe;
      do {
        if (chart_boundaries.count(ne)) pq.push(ne);
        ne = heds.getNext(ne);
      } while (ne != oe);
    } else if (find(f) != find(of) and
               distanceC[find(f)] - distanceF[of] < eps and
               distanceC[find(of)] - distanceF[f] < eps) {
      join(f, of);
    }
  }
}
std::map<int, std::vector<int>> Segmentation::getCharts() {
  map<int,vector<int>> m;
  for(int f = 0; f < heds.sizeOfFaces(); f++) {
    m[find(f)].push_back(f);
  }
  return m;
}

int Segmentation::countCharts() {
  int r = 0;
  for(int f = 0; f < heds.sizeOfFaces(); f++)
    if(find(f)==f) r++;
  return r;
}

void Segmentation::buildIncidentEdges() {
  incidentEdge.assign(heds.sizeOfFaces(), 0);
  for(int i = 0; i < heds.sizeOfHalfedges(); i++)
    incidentEdge[heds.getFace(i)] = i;
}

void Segmentation::init_unionfind() {
  par.assign(heds.sizeOfFaces(), -1);
  size.assign(heds.sizeOfFaces(), 1);
}

int Segmentation::find(int u) {
  if (u == -1) return -1;
  return par[u] == u ? u : par[u] = find(par[u]);
}

void Segmentation::join(int u, int v) {
  u = find(u);
  if (find(v) == -1) {
    par[v] = u;
    size[u]++;
    distanceC[u] = max(distanceC[u], distanceC[v]);
  } else {
    v = find(v);
    if (size[u] < size[v]) swap(u, v);
    par[v] = u;
    size[u] += size[v];
    distanceC[u] = max(distanceC[u], distanceC[v]);
  }
}

vector<int> Segmentation::maximalFacets() {
  vector<int> r;
  for (int f = 0; f < heds.sizeOfFaces(); f++)
    if (isFacetLocalMaximum(f,5)) r.emplace_back(f);
  return r;
}

bool Segmentation::isFacetLocalMaximum(int f, int depth) {
  set<int> marc;
  queue<pair<int,int>> q;
  q.push({f,0});
  marc.insert(f);
  bool ok = true;
  while(q.size()) {
    auto p = q.front(); q.pop();
    int ff = p.first, d = p.second;
    ok &= (distanceF[f] >= distanceF[ff]);
    int e = heds.getEdgeInFace(ff);
    int ne = e;
    if(d==depth) continue;
    do {
      int of = heds.getFace(heds.getOpposite(ne));
      ne = heds.getNext(ne);
      if(!marc.count(of)) {
        q.push({of,d+1});
        marc.insert(of);
      }
    } while (ne != e);
  }
  return ok;
}

void Segmentation::distanceToFeatures() {
  distanceF.assign(heds.sizeOfFaces(), 0x3f3f3f3f);

  queue<pair<int, int>> q;
  vector<bool> marc(heds.sizeOfHalfedges());
  for (int i = 0; i < heds.sizeOfHalfedges(); i++)
    if (tag[i] == 1) {
      q.push({i, 0});
      marc[i] = true;
    }
  while (q.size()) {
    auto p = q.front();
    q.pop();
    int e = p.first, d = p.second;
    int f = heds.getFace(e);

    distanceF[f] = min(distanceF[f], d);

    vector<int> neighbors = outNeighbors(e);

    for (int ne : neighbors) {
      if (!marc[ne]) {
        q.push({ne, d + 1});
        marc[ne] = true;
        eps = max(eps, (d+1)/4);
      }
    }
  }

  distanceC = distanceF;
}

void Segmentation::expandFeatureCurves() {
  tag.assign(heds.sizeOfHalfedges(), 0);
  for (auto &p : top_features) {
    expandFeatureCurve(p.second);
  }
}

void Segmentation::colorExpandedFeatures() {
  for (int i = 0; i < heds.sizeOfHalfedges(); i++)
    if (tag[i] == 1) colorEdge(i);
}

void Segmentation::colorInitialFeatures() {
  for (auto &p : top_features) colorEdge(p.second);
}

void Segmentation::findInitialFeatures() {
  top_features.clear();
  for (int i = 0; i < heds.sizeOfHalfedges(); i++) {
    top_features.insert({sharpness(i), i});
    if (top_features.size() > .05 * heds.sizeOfHalfedges())
      top_features.erase(top_features.begin());
  }
}

void Segmentation::expandFeatureCurve(int start) {
  if (tag[start]) return;

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
    if (hp != -1) {
      detected_vertices.insert(heds.getTarget(hp));
      detected_feature.push_back(hp);
    } else
      break;
  } while (best.sharpness > max_string_length * threshold);

  // backwards
  h = heds.getOpposite(start);
  hp = h;
  do {
    best = dfs(hp, 1, 0, -1, detected_vertices);
    hp = best.second;
    if (hp != -1) {
      detected_vertices.insert(heds.getTarget(hp));
      detected_feature.push_back(hp);
    } else
      break;
  } while (best.sharpness > max_string_length * threshold);

  if (detected_feature.size() > min_feature_length) {
    tagFeatures(detected_feature);
  }
}

dfsResult Segmentation::dfs(int h, int length, double sharp, int second,
                            set<int> &detected_feature) {
  int prev_tag = tag[h];
  tag[h] = 3;

  if (length == 2) second = h;
  sharp += sharpness(h);

  if (length == max_string_length) {
    tag[h] = prev_tag;
    return {sharp, second};
  }

  vector<int> out = outNeighbors(h);

  dfsResult best{1, -1};

  for (int nh : out) {
    if (nh == heds.getOpposite(h) || tag[nh] ||
        detected_feature.count(heds.getTarget(nh)))
      continue;
    dfsResult temp = dfs(nh, length + 1, sharp, second, detected_feature);
    if (temp.sharpness > best.sharpness) {
      best = temp;
    }
  }

  tag[h] = prev_tag;
  return best;
}

void Segmentation::tagFeatures(vector<int> &features) {
  queue<pair<int, int>> q;
  for (int f : features) {
    tagAsFeature(f);
    q.push({f, 0});
  }
  tagNeighborhoods(q, 5);
}

void Segmentation::tagAsFeature(int h) {
  tag[h] = tag[heds.getOpposite(h)] = 1;
}

void Segmentation::tagNeighborhoods(queue<pair<int, int>> &q, int depth) {
  // tag neighboors as neighboors
  while (q.size()) {
    auto p = q.front();
    q.pop();
    int h = p.first, d = p.second;

    if (d == depth) continue;

    vector<int> neighbors = outNeighbors(h);
    for (int n : neighbors) {
      if (!tag[n]) {
        tagAsNeighbor(n);
        q.push({n, d + 1});
      }
    }
  }
}

void Segmentation::tagAsNeighbor(int h) {
  if (!tag[h]) tag[h] = 2;
  h = heds.getOpposite(h);
  if (!tag[h]) tag[h] = 2;
}

double Segmentation::sharpness(int h) {
  int f1 = heds.getFace(h);
  int f2 = heds.getFace(heds.getOpposite(h));
  if (f1 == -1 or f2 == -1) return M_PI;
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
  viewer.data().add_edges(v1, v2, RowVector3d(1, 0, 0));
}

vector<int> Segmentation::outNeighbors(int h) {
  vector<int> neighbours;
  h = heds.getOpposite(h);
  int hn = h;
  do {
    neighbours.push_back(hn);
    hn = heds.getNext(heds.getOpposite(hn));
  } while (hn != h);
  return neighbours;
}
