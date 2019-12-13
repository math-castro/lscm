#pragma once

#include <igl/opengl/glfw/Viewer.h>

std::vector<Eigen::MatrixXd> parametrize(std::vector<const Eigen::MatrixXd*> &Vs, std::vector<const Eigen::MatrixXi*> &Ts);
Eigen::MatrixXd parametrize(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
double angleBetweenSides(const Eigen::RowVectorXd &a, const Eigen::RowVectorXd &b);
std::pair<int, int> approximateDiameter(const Eigen::MatrixXd &V);
void emplaceAB(std::vector<Eigen::Triplet<double>> &ta,
               std::vector<Eigen::Triplet<double>> &tb,
               int i,
               int j,
               int n,
               int np,
               double wr,
               double wi,
               std::pair<int, int> &diam);
