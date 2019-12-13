#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>

#include "halfedge/HalfedgeBuilder.cpp"
#include "Segmentation.cpp"
#include "parametrization/parametrization.hpp"
#include "packing/packing.hpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

// ------------ main program ----------------
int main(int argc, char *argv[]) {
  MatrixXd V;
  MatrixXi F;

  if (argc < 2) {
    igl::readOFF("../src/data/quad.off", V, F);
    //igl::readPLY("../src/data/bunny.ply", V, F);
//		igl::readOBJ("../src/data/LSCM_bunny.obj", V, F);
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
  cout << Us[0] << endl;

  viewer.data().set_mesh(*Vs[0], *Fs[0]);
  viewer.data().set_uv(Us[0]);
  viewer.data().show_texture = true;

//  MatrixX3d new_V = MatrixXd::Ones(V.rows(), 3);
//  rescale(V, U, F);
//  alignVertical(U);
//  alignBottomLeft(U);
//  new_V.col(0) = U.col(0);
//  new_V.col(1) = U.col(1);
//  viewer.data().set_mesh(new_V, F);
//  cout << U << endl;

  viewer.core(0).align_camera_center(V, F);
  viewer.launch(); // run the viewer
}
