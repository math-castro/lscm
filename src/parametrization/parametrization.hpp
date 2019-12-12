#pragma once

#include <igl/opengl/glfw/Viewer.h>

void parametrize(Eigen::MatrixXd &V, Eigen::MatrixXi &T);
double angleBetweenSides(const Eigen::RowVectorXd& a, const Eigen::RowVectorXd& b);