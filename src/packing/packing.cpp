#include "packing/packing.hpp"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <map>
#include <algorithm>
#include <chrono>

using namespace Eigen;
using namespace std;

void rescale(const MatrixXd &X, MatrixXd &U, const MatrixXi &T) {
  double s = sqrt(area(X,T)/area(U,T));
  U *= s;
}

double area(const MatrixXd &V, const MatrixXi &T) {
  double r = 0;
  for(int i = 0; i < T.rows(); i++) {
    double a = (V.row(T(i,1))-V.row(T(i,0))).norm();
    double b = (V.row(T(i,2))-V.row(T(i,1))).norm();
    double c = (V.row(T(i,0))-V.row(T(i,2))).norm();
    double p = (a+b+c)/2;
    r += sqrt(p*(p-a)*(p-b)*(p-c));
  }
  return r;
}

pair<int,int> findDiameter(const MatrixXd &V) {
  const int n = V.rows();
  double m = 0;
  pair<int,int> best{0,1};
  for(int i = 0; i < n; i++) {
    for(int j = i+1; j < n; j++) {
      double d = (V.row(i)-V.row(j)).norm();
      if(d>m) {
        m = d;
        best = pair<int,int>{i,j};
      }
    }
  }
  return best;
}

void alignVertical(Eigen::MatrixXd &U) {
  pair<int,int> diam = findDiameter(U);
  RowVectorXd u0 = U.row(diam.first);
  for(int i = 0; i < U.rows(); i++)
    U.row(i) -= u0;
  double a = acos(-1)/2 - atan2(U(diam.second,1), U(diam.second,0));
  MatrixXd R(2,2);
  R(0,0) = cos(a);
  R(0,1) = -sin(a);
  R(1,0) = sin(a);
  R(1,1) = cos(a);
  for(int i = 0; i < U.rows(); i++)
    U.row(i) = (R*(U.row(i).transpose())).transpose();
}

void alignBottomLeft(Eigen::MatrixXd &U) {
  // double min_u = DBL_MAX, min_v = DBL_MAX;
  double min_u = U.colwise().minCoeff()(0);
  double min_v = U.colwise().minCoeff()(1);
  for(int i = 0; i < U.rows(); i++) {
    U(i,0) -= min_u;
    U(i,1) -= min_v;
  }
}

Chart horizon(const MatrixXd &U, const MatrixXi &T, double resolution) {
  map<int, int> upper_horizon, lower_horizon;
  int max_y = 0;
  const int margin = 20;

  for(int i = 0; i < T.rows(); i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = j+1; k < 3; k++) {
        int x1 = (int)lround(U(T(i,j),0)/resolution);
        int y1 = (int)lround(U(T(i,j),1)/resolution);
        int x2 = (int)lround(U(T(i,k),0)/resolution);
        int y2 = (int)lround(U(T(i,k),1)/resolution);

        if(x1 > x2) {
          swap(x1,x2);
          swap(y1,y2);
        }

        for(int x = x1; x <= x2; x++) {
          int y = (x1==x2) ? y1 : y1 + (y2-y1)*(x-x1)/(x2-x1);
          if(lower_horizon.count(x)) {
            lower_horizon[x] = min(lower_horizon[x], y-margin);
            upper_horizon[x] = max(upper_horizon[x], y+margin);
          }
          else {
            lower_horizon[x] = y-margin;
            upper_horizon[x] = y+margin;
          }
        }
      }
    }
  }
  int max_x = lower_horizon.rbegin()->first;

  vector<int> lh(max_x+1+margin),uh(max_x+1+margin);

  // cout << max_x << endl;
  for(int i = 0; i <= max_x; i++) {
    if(lower_horizon.count(i)) {
      lh[i] = lower_horizon[i];
      uh[i] = upper_horizon[i];
    }
    else {
      auto next = lower_horizon.upper_bound(i);
      auto prev = next; --prev;
      lh[i] = prev->second + (next->second-prev->second)*(i-prev->first)/(next->first-prev->first);

      next = upper_horizon.upper_bound(i);
      prev = next; --prev;
      uh[i] = prev->second + (next->second-prev->second)*(i-prev->first)/(next->first-prev->first);
    }
    // if(max_x < 25)
    //   cout << i << " " << lh[i] << " " << uh[i] << endl;
  }
  for(int i = max_x+1; i < max_x+1+margin; i++) {
    lh[i] = lh[max_x];
    uh[i] = uh[max_x];
  }

  return {lh, uh, max_y, max_x};
}

vector<int> pack(vector<const MatrixXd*> Xs, vector<MatrixXd*> Us, vector<const MatrixXi*> Ts) {
  cout << "Started packing:..." << flush;
  auto start = chrono::high_resolution_clock::now();
  vector<double> d;
  vector<int> id;
  vector<Chart> charts;
  d.reserve(Xs.size());
  id.reserve(Xs.size());
  charts.reserve(Xs.size());

  for(int i = 0; i < Us.size(); i++) {
    const MatrixXd &X = *Xs[i];
    MatrixXd &U = *Us[i];
    const MatrixXi &T = *Ts[i];

    rescale(X,U,T);
    alignVertical(U);
    alignBottomLeft(U);
  }

  double resolution = DBL_MAX;

  for(auto pU : Us) {
    auto &U = *pU;
    resolution = min(resolution, U.colwise().maxCoeff()(0)/10);
  }

  for(int i = 0; i < Us.size(); i++) {
    MatrixXd &U = *Us[i];
    const MatrixXi &T = *Ts[i];

    charts.emplace_back(horizon(U, T, resolution));

    auto diam = findDiameter(U);
    d.emplace_back((U.row(diam.first) - U.row(diam.second)).norm());
    id.emplace_back(i);
  }

  auto comp = [&](int a, int b) {return d[a] > d[b];};
  sort(id.begin(), id.end(), comp);

  double size = sqrt(totalArea(Us, Ts));
  vector<int> hor;
  vector<int> dx, dy;
  FitResult fr;

  double l = sqrt(totalArea(Us, Ts));
  double r = 10*l;
  while(fabs(l-r)/l > 1e-2) {
    double m = (r+l)/2;
    // cout << l << " " << r << endl;
    fr = canFit(m, resolution, charts, id);
    if(fr.ok) r = m;
    else l = m;
  }
    fr = canFit(r, resolution, charts, id);

  for (int i = 0; i < id.size(); i++) {
    translateU(*Us[id[i]], fr.dx[i], fr.dy[i], resolution);
  }

  auto finish = chrono::high_resolution_clock::now();
  cout << "done: " << chrono::duration<double>(finish-start).count() << " s" << endl;
  return id;
}

FitResult canFit(double size, double resolution, vector<Chart> &charts, vector<int> &id) {
  bool ok = true;
  vector<int> hor(size / resolution);
  vector<int> dx, dy;
  // cout << int(size/resolution) << endl;

  for (int i : id) {
    pair<int, int> xy = calculateBestXY(hor, charts[i], size / resolution);
    int x = xy.first, y = xy.second;
    if (x == -1) {
      ok = false;
      break;
    }
    updateHorizon(hor, charts[i], x, y);
    dx.emplace_back(x);
    dy.emplace_back(y);
  }

  return {dx,dy,ok};
}

int calculateDY(vector<int> &hor, vector<int> &lh, int x) {
  int dy = 0;
  for(int i = 0; i < lh.size(); i++) {
    if(x+i < hor.size())
      dy = max(dy, hor[x+i]-lh[i]);
    else
      throw "Out of bounds";
  }
  return dy;
}

int calculateAreaUnder(vector<int> &hor, vector<int> &lh, int x, int y) {
  int a = 0;
  for(int i = 0; i < lh.size(); i++) {
    if(lh[i]+y-hor[x+i] < 0)
      throw "Fudeu";
    if(x+i < hor.size())
      a += lh[i]+y-hor[x+i];
    else
      throw "Out of bounds";
  }
  return a;
}

pair<int,int> calculateBestXY(vector<int> &hor, Chart &chart, int MAX_Y) {
  int min_area = 0x3f3f3f3f;
  int best_x = 0;
  int best_y = 0;
  vector<int> &lh = chart.lh;
  for(int i = 0; i + lh.size() <= hor.size(); i++) {
    int dy = calculateDY(hor, lh, i);
    if(dy+chart.max_y >= MAX_Y) continue;
    int a = calculateAreaUnder(hor, lh, i, dy);
    if(a < min_area) {
      min_area = a;
      best_x = i;
      best_y = dy;
    }
  }
  if(min_area == 0x3f3f3f3f) return {-1,-1};
  return {best_x, best_y};
}

void updateHorizon(vector<int> &hor, Chart &chart, int x, int y) {
  vector<int> &uh = chart.uh;
  for(int i = 0; i < uh.size(); i++) {
    if(x+i < hor.size())
      hor[x+i] = uh[i]+y;
    else
      throw "Out of bounds";
  }
}

void translateU(MatrixXd &U, int x, int y, double resolution) {
  double dx = resolution*x, dy = resolution*y;
  for(int i = 0; i < U.rows(); i++) {
    U(i,0) += dx;
    U(i,1) += dy;
  }
}

double totalArea(std::vector<Eigen::MatrixXd*> Us, std::vector<const Eigen::MatrixXi*> Ts) {
  double a = 0;
  for(int i = 0; i < Us.size(); i++)
    a += area(*Us[i], *Ts[i]);
  return a;
}
