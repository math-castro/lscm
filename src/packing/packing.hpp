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
void pack(std::vector<const Eigen::MatrixXd*>Xs, std::vector<Eigen::MatrixXd*> Us, std::vector<const Eigen::MatrixXi*> Ts);
int calculateDY(std::vector<int> &hor, std::vector<int> &lh, int x);
int calculateAreaUnder(std::vector<int> &hor, std::vector<int> &lh, int x, int y);
std::pair<int,int> calculateBestXY(std::vector<int> &hor, std::vector<int> &lh);
void updateHorizon(std::vector<int> &hor, std::vector<int> &uh, int x, int y);
void translateU(Eigen::MatrixXd &U, int x, int y, double resolution);
