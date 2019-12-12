#include <igl/opengl/glfw/Viewer.h>
#include <math.h>

#ifndef HALFEDGE_DS_HEADER
	#define HALFEDGE_DS_HEADER
	#include "halfedge/HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;


class Segmentation {

public:
	/**
	 * Initialize the data structures
	 **/
	Segmentation(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh) {
		threshold = 0.05;
		he = &mesh;
		V = &V_original;
		F = &F_original;
		nEdges = he->sizeOfHalfedges() / 2; // number of edges
		//int n = V_original.rows();		   // number of vertices
		EdgeSharpness = new VectorXd(nEdges);
		getEdgeSharpnessMatrix();
	}

	Vector3d getNormal(int f) {
		// returns normal vector to face f
		Vector3d v1 = V->row((*F)(f, 1)) - V->row((*F)(f, 0));
		Vector3d v2 = V->row((*F)(f, 2)) - V->row((*F)(f, 0));
		return v1.cross(v2).normalized();
	}

	float getEdgeSharpness(int e) {
		// returns angle between normals to the faces adjacent to edge e
		int f1 = he->getFace(e);
		int f2 = he->getFace(he->getOpposite(e));
		if (f1 == -1 || f2 == -1) {
			return M_PI;
		}
		Vector3d n1 = getNormal(f1);
		Vector3d n2 = getNormal(f2);
		return acos(n1.dot(n2));
	}

	void getEdgeSharpnessMatrix() {
		//Assign a number, between 1..nEdges, to all halfedges
		int* EdgeMap = new int[2*nEdges]();
		for (int e=0; e<2*nEdges; e++) {
			EdgeMap[e] = -1;
		}
		int i = 0;
		for (int e=0; e<2*nEdges; e++) {
			if (EdgeMap[e] != -1) {
				continue;
			}
			EdgeMap[e] = i;
			EdgeMap[he->getOpposite(e)] = i;
			i++;
		}
		//Fill the EdgeSharpness Matrix with the sharpness criterion (angle between normals) for each edge
		for (int e=0; e<2*nEdges; e++) {
			int i = EdgeMap[e];
			EdgeSharpness->row(i) << getEdgeSharpness(e);
		}
		//cout << *EdgeSharpness << endl;
	}

	void setThreshold(float newThreshold) { threshold = newThreshold; }
	float getThreshold() { return threshold; }

	vector<int> getNeighbours(int v) {
		return vector<int>();
	}

	void expandFeatureCurve(int startEdge) {
		int* tag = new int[2*nEdges](); // tag[i] = 0 (nothing), 1 (feature), 2 (feature neighbor)
		vector<int> detectedFeature;
		int maxStringLength = 5;
		int minFeatureLength = 15;
		int edges[] = {startEdge, he->getOpposite(startEdge)};
		for (const int &edge : edges) {
    		int current = edge;
			float sharpness = 0;
			vector<int> S;
			do {
				S.push_back(current);
				sharpness += getEdgeSharpness(current);
				//DFS(current, S, sharpness);
				current = S[1];
				detectedFeature.push_back(current);

			} while (sharpness > maxStringLength * threshold);
		}
		if (detectedFeature.size() > minFeatureLength) {
			// tag elements of detectedFeature as feature
			for (int i=0; i<detectedFeature.size(); i++) {
				tag[detectedFeature[i]] = 1;
			}
			// tag neighbours of detectedFeature as feature neighbor
			int start = he->getTarget(he->getOpposite(detectedFeature[0]));
			vector<int> neighbours = getNeighbours(start);
			for (int j=0; j<neighbours.size(); j++) {
				if (tag[neighbours[j]] == 0) {
					tag[neighbours[j]] = 2;
				}
			}
			for (int i=0; i<detectedFeature.size(); i++) {
				int current = he->getTarget(detectedFeature[i]);
				vector<int> neighbours = getNeighbours(start);
				for (int j=0; j<neighbours.size(); j++) {
					if (tag[neighbours[j]] == 0) {
						tag[neighbours[j]] = 2;
					}
				}
			}
		}
	}

private:
	float threshold;
	/** Half-edge representation of the original input mesh */
	HalfedgeDS *he;
	MatrixXd *V; // vertex coordinates of the original input mesh
	MatrixXi *F; // REMARK: not needed if using the half-edge data structure
	VectorXd *EdgeSharpness;

	int nEdges;
	//int nVertices, nFaces; // number of vertices, faces
};
