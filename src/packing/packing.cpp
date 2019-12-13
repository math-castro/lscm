#include "packing/packing.hpp"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <map>

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

pair<vector<double>, vector<double>> horizon(const MatrixXd &U, double resolution) {
  double min_u = DBL_MAX, max_u = -DBL_MAX;
  map<int, double> upper_horizon, lower_horizon;
  for(int i = 0; i < U.rows(); i++) {
    min_u = min(min_u, U(i, 0));
  }
  for(int i = 0; i < U.rows(); i++) {
    int x = lround((U(i,0)-min_u)/resolution);
    if(lower_horizon.count(x)) {
      lower_horizon[x] = min(lower_horizon[x], U(i,1));
      upper_horizon[x] = max(upper_horizon[x], U(i,1));
    }
    else {
      lower_horizon[x] = U(i,1);
      upper_horizon[x] = U(i,1);
    }
  }
  int max_x = lower_horizon.rbegin()->first;

  vector<double> lh(max_x+1),uh(max_x+1);

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
  }

  return pair<vector<double>, vector<double>>(lh, uh);
}
