#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>

#include "halfedge/HalfedgeBuilder.hpp"
#include "segmentation/segmentation.hpp"
#include "parametrization/parametrization.hpp"
#include "packing/packing.hpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

// ------------ main program ----------------
int main(int argc, char *argv[]) {
  MatrixXd V;
  MatrixXi F;

  if (argc < 2) {
    // igl::readOFF("../src/data/quad.off", V, F);
    //igl::readPLY("../src/data/bunny.ply", V, F);
		igl::readOBJ("../src/data/LSCM_bunny.obj", V, F);
  } else {
    std::cout << "reading input file: " << argv[1] << std::endl;
    igl::readOFF(argv[1], V, F);
  }

  igl::opengl::glfw::Viewer viewer; // create the 3d viewer

//  viewer.data().set_mesh(V, F);
//  viewer.append_mesh(true);

//	HalfedgeBuilder *builder = new HalfedgeBuilder();
//	HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
//	Segmentation *segmentation = new Segmentation(V, F, he, viewer);
  vector<const MatrixXd*> Vs = {&V};
  vector<const MatrixXi*> Fs = {&F};
  vector<MatrixXd> Us = parametrize(Vs, Fs);

 viewer.data().set_mesh(*Vs[0], *Fs[0]);
 viewer.data().set_uv(Us[0]);
 viewer.data().show_texture = true;

  typedef Matrix<unsigned char, -1, -1> MatrixXc;
  MatrixXc R, G, B;
  const int n = 10;
  R = 255*MatrixXc::Ones(n,n);
  G = 255*MatrixXc::Ones(n,n);
  B = 255*MatrixXc::Ones(n,n);

  for(int i = 0; i < n; i+=10)
    for(int j = 0; j < n; j++)
      G(i,j) = B(i, j) = 0;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j+=10)
      G(i,j) = R(i, j) = 0;


 viewer.data().set_texture(R,G,B);
 viewer.data().set_colors(RowVector3d(1,1,1));
 viewer.core().lighting_factor = 0;


  vector<MatrixXd*> pUs;
  pUs.reserve(Us.size());
  for(auto & U : Us) pUs.emplace_back(&U);

  pack(Vs, pUs, Fs);

  for(int i = 0; i < Us.size(); i++) {
    auto & U = Us[i];
    auto & T = *Fs[i];

    viewer.append_mesh();

    MatrixXd new_U = MatrixXd::Zero(U.rows(), 3);
    new_U.col(0) = U.col(0);
    new_U.col(1) = U.col(1);

    viewer.data().set_mesh(new_U,T);
//    viewer.data().set_texture(R,G,B);
//    viewer.data().show_texture = true;
    viewer.data().set_colors(RowVector3d(1,1,1));
  }

  viewer.core(0).align_camera_center(V, F);
  viewer.launch(); // run the viewer
}
