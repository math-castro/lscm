#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include <vector>
#include <utility>

struct Chart {
  std::vector<int> lh;
  std::vector<int> uh;
  int max_y;
  int max_x;
};

struct FitResult {
  std::vector<int> dx;
  std::vector<int> dy;
  bool ok;
};

void rescale(const Eigen::MatrixXd &X, Eigen::MatrixXd &U, const Eigen::MatrixXi &T);
double area(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);
std::pair<int, int> findDiameter(const Eigen::MatrixXd &V);
void alignVertical(Eigen::MatrixXd &U);
void alignBottomLeft(Eigen::MatrixXd &U);
Chart horizon(const Eigen::MatrixXd &U, double resolution);
std::vector<int> pack(std::vector<const Eigen::MatrixXd*>Xs, std::vector<Eigen::MatrixXd*> Us, std::vector<const Eigen::MatrixXi*> Ts);
FitResult canFit(double size, double resolution, std::vector<Chart> &charts, std::vector<int> &id);
int calculateDY(std::vector<int> &hor, std::vector<int> &lh, int x);
int calculateAreaUnder(std::vector<int> &hor, std::vector<int> &lh, int x, int y);
std::pair<int,int> calculateBestXY(std::vector<int> &hor, Chart &chart, int MAX_Y);
void updateHorizon(std::vector<int> &hor, Chart &chart, int x, int y);
void translateU(Eigen::MatrixXd &U, int x, int y, double resolution);
double totalArea(std::vector<Eigen::MatrixXd*> Us, std::vector<const Eigen::MatrixXi*> Ts);
