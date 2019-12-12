#include "parametrization/parametrization.hpp"

#include <vector>
#include <cmath>

using namespace Eigen;
using namespace std;

typedef Triplet<double> Td;

MatrixX2d parametrize(MatrixXd &V, MatrixXi &T) {
  const int np = T.rows(), n = V.rows();

  // List of triplets for sparse matrix creation
  vector<Td> ta, tb, tu;
  ta.reserve(np * 12);

  // For each triangle
  for (int i = 0; i < np; i++) {
    // Calculate x,y
    // First vertex as (0,0), second as (s,0), third as t(cos(a), sin(a))
    RowVectorXd p[3];
    double x[3], y[3];

    p[0] = V.row(T(i, 0));
    p[1] = V.row(T(i, 1));
    p[2] = V.row(T(i, 2));

    double a = angleBetweenSides(p[1] - p[0], p[2] - p[0]);
    double s = (p[1] - p[0]).norm();
    double t = (p[2] - p[0]).norm();

    x[0] = y[0] = 0;
    x[1] = s, y[1] = 0;
    x[2] = t * cos(a), y[2] = t * sin(a);

    double d = sqrt((x[0] * y[1] - y[0] * x[1]) + (x[1] * y[2] - y[1] * x[2]) + (x[2] * y[0] - y[2] * x[0]));
//    d=1;

    // Calculate W
    double Wr[3], Wi[3];
    Wr[0] = (x[2] - x[1]) / d, Wi[0] = (y[2] - y[1]) / d;
    Wr[1] = (x[0] - x[2]) / d, Wi[1] = (y[0] - y[2]) / d;
    Wr[2] = (x[1] - x[0]) / d, Wi[2] = (y[1] - y[0]) / d;


    // Push values to A and B
    if (T(i, 0) < n - 2) {
      int j = T(i, 0);
      ta.emplace_back(i, j, Wr[0]);
      ta.emplace_back(i + np, j + n - 2, Wr[0]);
      ta.emplace_back(i, j + n - 2, - Wi[0]);
      ta.emplace_back(i + np, j, Wi[0]);
    } else {
      int j = T(i, 0) - n + 2;
      tb.emplace_back(i, j, Wr[0]);
      tb.emplace_back(i + np, j + 2, Wr[0]);
      tb.emplace_back(i, j + 2, -Wi[0]);
      tb.emplace_back(i + np, j, Wi[0]);
    }
    if (T(i, 1) < n - 2) {
      int j = T(i, 1);
      ta.emplace_back(i, j, Wr[1]);
      ta.emplace_back(i + np, j + n - 2, Wr[1]);
      ta.emplace_back(i, j + n - 2, - Wi[1]);
      ta.emplace_back(i + np, j, Wi[1]);
    } else {
      int j = T(i, 1) - n + 2;
      tb.emplace_back(i, j, Wr[1]);
      tb.emplace_back(i + np, j + 2, Wr[1]);
      tb.emplace_back(i, j + 2, - Wi[1]);
      tb.emplace_back(i + np, j, Wi[1]);
    }
    if (T(i, 2) < n - 2) {
      int j = T(i, 2);
      ta.emplace_back(i, j, Wr[2]);
      ta.emplace_back(i + np, j + n - 2, Wr[2]);
      ta.emplace_back(i, j + n - 2, -Wi[2]);
      ta.emplace_back(i + np, j, Wi[2]);
    } else {
      int j = T(i, 2) - n + 2;
      tb.emplace_back(i, j, Wr[2]);
      tb.emplace_back(i + np, j + 2, Wr[2]);
      tb.emplace_back(i, j + 2, -Wi[2]);
      tb.emplace_back(i + np, j, Wi[2]);
    }
  }

  // Build A and B
  SparseMatrix<double> A(2 * np, 2 * (n - 2));
  SparseMatrix<double> B(2 * np, 4);
  A.setFromTriplets(ta.begin(), ta.end());
  ta.clear();
  B.setFromTriplets(tb.begin(), tb.end());
  tb.clear();

  // Build up 
  VectorXd up(4);
  up(0) = 0;
  up(1) = 1;
  up(2) = 0;
  up(3) = 1;

  VectorXd b = - B * up;

  // Solve least squares using QR decomposition
  A.makeCompressed();
  SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver(A);
  VectorXd uf = solver.solve(b);

  // Join uf and up
  MatrixX2d U(n, 2);
  for (int i = 0; i < n - 2; i++) {
    U(i, 0) = uf(i);
    U(i, 1) = uf(i + n - 2);
  }
  U(n - 2, 0) = U(n - 2, 1) = 0;
  U(n - 1, 0) = U(n - 1, 1) = 1;

  return U;
}

double angleBetweenSides(const RowVectorXd &a, const RowVectorXd &b) {
  return acos(a.dot(b) / a.norm() / b.norm());
}