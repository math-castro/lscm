#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/writeOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include "halfedge/HalfedgeBuilder.cpp"
#include "Segmentation.cpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
	if (key == '1') {
		HalfedgeBuilder *builder = new HalfedgeBuilder();
		HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
		//std::cout << "V(" << V.rows() << ", " << V.cols() << "): " << std::endl;
		//std::cout << V << std::endl;
		//std::cout << "F(" << F.rows() << ", " << F.cols() << "): " << std::endl;
		//std::cout << F << std::endl;
		return true;
	}
	return false;
}


// ------------ main program ----------------
int main(int argc, char *argv[]) {

	if (argc < 2) {
		//igl::readOFF("../src/data/bunny.off", V, F);
		//igl::readPLY("../src/data/bunny.ply", V, F);
		igl::readOBJ("../src/data/LSCM_bunny.obj", V, F);
	} else {
		std::cout << "reading input file: " << argv[1] << std::endl;
		igl::readOFF(argv[1], V, F);
	}

	igl::opengl::glfw::Viewer viewer; // create the 3d viewer

	viewer.callback_key_down = &key_down;
	viewer.data().set_mesh(V, F);

	HalfedgeBuilder *builder = new HalfedgeBuilder();
	HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
	Segmentation *segmentation = new Segmentation(V, F, he, viewer);

	viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer

}
