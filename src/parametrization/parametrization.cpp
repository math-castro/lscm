#include "parametrization/parametrization.hpp"

#include <vector>
#include <cmath>

using namespace Eigen;
using namespace std;

typedef Triplet<double> Td;

void parametrize(MatrixXd &V, MatrixXi &T) {
  // List of triplets for sparse matrix creation
  vector<Td> ta, tb, tu;
  ta.reserve(T.rows()*6);
  tb.reserve(T.rows()*6);

  // For each triangle
  for(int i = 0; i < T.rows(); i++) {
    // Calculate x,y
    // First vertex as (0,0), second as (s,0), third as t(cos(a), sin(a))
    RowVectorXd p[3];
    double x[3], y[3];

    p[0] = V.row(T(i,0));
    p[1] = V.row(T(i,1));
    p[2] = V.row(T(i,2));
    
    x[0]=y[0]=0;
    x[1]=(p[1]-p[0]).norm(), y[1] = 0;
    double a = angleBetweenSides(p[1]-p[0], p[2]-p[0]);
    double t = (p[2]-p[0]).norm();
    x[2] = t*cos(a), y[2] = t*sin(a);
    double d = (x[0]*y[1]-y[0]*x[1]) + (x[1]*y[2]-y[1]*x[2]) + (x[2]*y[0]-y[2]*x[0]);

    // Calculate W
    double Wr[3], Wi[3];
    Wr[0] = x[2]-x[1]; 
    Wr[1] = x[0]-x[2]; 
    Wr[2] = x[1]-x[0];
    Wi[0] = y[2]-y[1];
    Wi[1] = y[0]-y[2];
    Wi[2] = y[1]-y[0];

    const int np = T.rows(), n = V.rows();

    // Push values to A and B
    if(T(i,0) < V.rows()-2) {
      int j = T(i,0);
      ta.push_back(Td(i, j, Wr[0]/d));
      ta.push_back(Td(i+np, j+n-2, Wr[0]/d));
      ta.push_back(Td(i, j+n-2, -Wi[0]/d));
      ta.push_back(Td(i+np, j, Wi[0]/d));
    } 
    else {
      int j = T(i,0)-n+2;
      tb.push_back(Td(i, j, Wr[0]/d));
      tb.push_back(Td(i+np, j+2, Wr[0]/d));
      tb.push_back(Td(i, j+2, -Wi[0]/d));
      tb.push_back(Td(i+np, j, Wi[0]/d));
    }
    if(T(i,1) < V.rows()-2) {
      int j = T(i,1);
      ta.push_back(Td(i, j, Wr[1]/d));
      ta.push_back(Td(i+np, j+n-2, Wr[1]/d));
      ta.push_back(Td(i, j+n-2, -Wi[1]/d));
      ta.push_back(Td(i+np, j, Wi[1]/d));
    }
    else {
      int j = T(i,1)-n+2;
      tb.push_back(Td(i, j, Wr[1]/d));
      tb.push_back(Td(i+np, j+2, Wr[1]/d));
      tb.push_back(Td(i, j+2, -Wi[1]/d));
      tb.push_back(Td(i+np, j, Wi[1]/d));
    }
    if(T(i,2) < V.rows()-2) {
      int j = T(i,2);
      ta.push_back(Td(i, j, Wr[2]/d));
      ta.push_back(Td(i+np, j+n-2, Wr[2]/d));
      ta.push_back(Td(i, j+n-2, -Wi[2]/d));
      ta.push_back(Td(i+np, j, Wi[2]/d));
    }
    else {
      int j = T(i,2)-n+2;
      tb.push_back(Td(i, j, Wr[2]/d));
      tb.push_back(Td(i+np, j+2, Wr[2]/d));
      tb.push_back(Td(i, j+2, -Wi[2]/d));
      tb.push_back(Td(i+np, j, Wi[2]/d));
    }
  }
  
  // Build A and B
  SparseMatrix<double> A(2*T.rows(), V.rows()-2);
  SparseMatrix<double> B(2*T.rows(), 4);
  A.setFromTriplets(ta.begin(), ta.end());
  ta.clear();
  B.setFromTriplets(tb.begin(), tb.end());
  tb.clear();

  // Build Up 
  VectorXd Up(4);
  Up(0) = 0;
  Up(1) = 1;
  Up(2) = 0;
  Up(3) = 1;

  VectorXd b = B*Up;

  // Solve least squares using QR decomposition
  A.makeCompressed();
  SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver(A);
  VectorXd u = solver.solve(b);
  
}

double angleBetweenSides(const RowVectorXd& a, const RowVectorXd& b) {
  return acos(a.dot(b)/a.norm()/b.norm());
}