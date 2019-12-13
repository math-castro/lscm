#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include <vector>
#include <utility>

void rescale(const Eigen::MatrixXd &X, Eigen::MatrixXd &U, const Eigen::MatrixXi &T);
double area(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
std::pair<int, int> findDiameter(const Eigen::MatrixXd &V);
void alignVertical(Eigen::MatrixXd &U);
void alignBottomLeft(Eigen::MatrixXd &U);
std::pair<std::vector<int>, std::vector<int>> horizon(const Eigen::MatrixXd &U, double resolution);
