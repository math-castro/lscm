#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <iostream>

#include "halfedge/HalfedgeBuilder.hpp"
#include "packing/packing.hpp"
#include "parametrization/parametrization.hpp"
#include "segmentation/segmentation.hpp"

using namespace Eigen;  // to use the classes provided by Eigen library
using namespace std;

// ------------ main program ----------------
int main(int argc, char *argv[]) {
  int mode = atoi(argv[argc-1]);

  MatrixXd V;
  MatrixXi F;

  if (argc < 2) {
    // igl::readOFF("../src/data/quad.off", V, F);
    // igl::readPLY("../src/data/bunny.ply", V, F);
    igl::readOBJ("../src/data/LSCM_bunny.obj", V, F);
  } else if (argc == 2) {
	std::cout << "reading input file: " << argv[1] << std::endl;
    igl::readOFF(argv[1], V, F);
  } else {
    std::cout << "reading input file: " << argv[1] << std::endl;
    igl::readOBJ(argv[1], V, F);
  }

  igl::opengl::glfw::Viewer viewer;  // create the 3d viewer

  // SEGMENTATION 

  Segmentation segmentation(V, F, viewer);
  auto m = segmentation.run(mode);
  cout << "Number of charts: " << m.size() << endl;

  // BUILD CHARTS

  vector<MatrixXd> Vs;
  vector<MatrixXi> Fs;
  Vs.reserve(m.size());
  Fs.reserve(m.size());

  for(auto &p : m) {
    auto &v = p.second;
    map<int,int> id;
    for(int f : v) {
      for(int i = 0; i < 3; i++) {
        int x = F(f,i);
        if(id.count(x) == 0)
          id[x] = id.size();
      }
    }
    MatrixXd Vp(id.size(), 3);
    MatrixXi Fp(v.size(), 3);
    for(int j = 0; j < v.size(); j++) {
      int f = v[j];
      for(int i = 0; i < 3; i++)
        Fp(j,i) = id[F(f,i)];
    }
    for(auto &pp : id) {
      int x = pp.first, xp = pp.second;
      Vp.row(xp) = V.row(x);
    }
    Vs.push_back(Vp);
    Fs.push_back(Fp);
    v.clear();
  }

  // SHOW CHARTS

  if(mode == 2) {
    for(int i = 0; i < Vs.size(); i++) {
      viewer.append_mesh(true);
      viewer.data().set_mesh(Vs[i], Fs[i]);
      auto color = RowVector3d::Random();
      viewer.data().set_colors(color);
      viewer.data().show_lines = false;
    }
    viewer.launch();
    exit(EXIT_SUCCESS);
  }

  // PARAMETRIZATION

  vector<const MatrixXd *> pVs;
  vector<const MatrixXi *> pFs;
  pVs.reserve(Vs.size());
  pFs.reserve(Fs.size());
  for (auto &V : Vs) pVs.emplace_back(&V);
  for (auto &F : Fs) pFs.emplace_back(&F);

  vector<MatrixXd> Us = parametrize(pVs, pFs);

  vector<MatrixXd *> pUs;
  pUs.reserve(Us.size());
  for (auto &U : Us) pUs.emplace_back(&U);

  // BUILD ISO U,V LINES 

  typedef Matrix<unsigned char, -1, -1> MatrixXc;
  MatrixXc R, G, B;
  const int n = 100;
  R = 255 * MatrixXc::Ones(n, n);
  G = 255 * MatrixXc::Ones(n, n);
  B = 255 * MatrixXc::Ones(n, n);

  for (int i = 0; i < n; i += 10)
    for (int j = 0; j < n; j++) G(i, j) = B(i, j) = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j += 10) G(i, j) = R(i, j) = 0;

  
  // SHOW ISO CURVES
  if(mode == 3) {
    for(int i = 0; i < Vs.size(); i++) {
      viewer.append_mesh(true);
      viewer.data().set_mesh(Vs[i], Fs[i]);
      viewer.data().set_colors(RowVector3d(1,1,1));
      viewer.data().set_texture(R,G,B);
      viewer.data().set_uv(Us[i]);
      viewer.data().show_texture = true;
      viewer.data().show_lines = false;
    }
    viewer.launch();
    exit(EXIT_SUCCESS);
  }

  // PACKING

  vector<int> id = pack(pVs, pUs, pFs);
  reverse(id.begin(), id.end());

  for (int j = 0; j < Us.size(); j++) {
    int i = id[j];
    auto &U = Us[i];
    auto &T = Fs[i];

    viewer.append_mesh();

    MatrixXd new_U = 0*MatrixXd::Ones(U.rows(), 3);
    new_U.col(0) = U.col(0);
    new_U.col(1) = U.col(1);

    viewer.data().set_mesh(new_U,T);
    viewer.data().set_colors(RowVector3d(1,1,1));
  }

  // viewer.core().align_camera_center(V, F);
  viewer.launch();
}
