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
	Segmentation(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh, igl::opengl::glfw::Viewer &viewer_) {
		viewer = &viewer_;
		//threshold = 0.7; // low poly bunny;
		threshold = 0.838; // for the high resolution bunny, allows 5% of the edges.
		he = &mesh;
		V = &V_original;
		F = &F_original;
		nEdges = he->sizeOfHalfedges() / 2; // number of edges
		//int n = V_original.rows();		   // number of vertices
		EdgeSharpness = new VectorXd(nEdges);
		getEdgeSharpnessMatrix();
		//cout << *EdgeSharpness << endl;
		//cout << nEdges << endl;
		colorSharpEdges();
	}

	void colorSharpEdges() {
		//viewer.append_mesh();
		cout << "show_lines " << viewer->data().show_lines << endl;
		viewer->data().show_lines = false;
		cout << "show_lines " << viewer->data().show_lines << endl;

		/*cout << "point_size " << viewer->data().point_size << endl;
		viewer->data().point_size = 60.0f;
		cout << "point_size " << viewer->data().point_size << endl;

		cout << "linewidth " << viewer->data().line_width << endl;
		viewer->data().line_width = 500.0f;
		cout << "linewidth " << viewer->data().line_width << endl;*/
		MatrixXd e1(2*nEdges, 3);
		MatrixXd e2(2*nEdges, 3);
		int i = 0;
		for (int e=0; e<2*nEdges; e++) {
			if ((*EdgeSharpness)(EdgeMap[e]) < threshold) {
				continue;
			}

			e1.row(i) << V->row(he->getTarget(e));
			e2.row(i) << V->row(he->getTarget(he->getOpposite(e)));
			i++;
		}
		//cout << "linewidth " << viewer->data().line_width << endl;
		viewer->data().add_edges(
			e1,
			e2,
			Eigen::RowVector3d(1, 0, 0));

		/*MatrixX3d Vc;
		Vc << -1, -1, -1,
			  -1, -1,  1,
			  -1,  1, -1,
	  		  -1,  1,  1,
			   1, -1, -1,
	  		   1, -1,  1,
	  		   1,  1, -1,
	  	  	   1,  1,  1;
		MatrixXi Fc;
		Fc << 0, 1, 2,
			  1, 3, 2,
			  4, 0, 6,
			  0, 2, 6,
			  5, 4, 7,
			  4, 6, 7,
			  1, 5, 3,
			  5, 7, 3,
			  2, 3, 6,
			  3, 7, 6,
			  4, 5, 0,
			  5, 1, 0;
		viewer->data().*/

		/*cout << viewer->data(0).line_width << endl;
		viewer->data(1).line_width = 1000.0f;
		cout << viewer->data(1).line_width << endl;
		cout << viewer->data(10).line_width << endl;
		cout << viewer->data(10000000).line_width << endl;
		*/
		//cout << "linewidth " << viewer->data().line_width << endl;
		//viewer->data().line_width = 300.0f;
		//cout << "linewidth " << viewer->data().line_width << endl;
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
		// Assign a number, between 1..nEdges, to all halfedges
		EdgeMap = new int[2*nEdges]();
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
		// Fill the EdgeSharpness Matrix with the sharpness criterion (angle between normals) for each edge
		for (int e=0; e<2*nEdges; e++) {
			int i = EdgeMap[e];
			EdgeSharpness->row(i) << getEdgeSharpness(e);
		}
		//cout << *EdgeSharpness << endl;
	}

	void setThreshold(float newThreshold) { threshold = newThreshold; }
	float getThreshold() { return threshold; }

	vector<int> getInNeighbours(int v) {// in-going halfedges incident to vertex v
		vector<int> neighbours = vector<int>();
		int startEdge = he->getEdge(v);
		neighbours.push_back(startEdge);
		int edge = he->getOpposite(he->getNext(startEdge));
		while (edge != startEdge) {
			edge = he->getOpposite(he->getNext(edge));
			neighbours.push_back(edge);
		}
		return neighbours;
	}

	vector<int> getOutNeighbours(int v) {// out-going halfedges incident from vertex v
		vector<int> neighbours = vector<int>();
		int startEdge = he->getOpposite(he->getEdge(v));
		neighbours.push_back(startEdge);
		int edge = he->getNext(he->getOpposite(startEdge));
		while (edge != startEdge) {
			edge = he->getNext(he->getOpposite(edge));
			neighbours.push_back(edge);
		}
		return neighbours;
	}

	pair<float, vector<int>> DFS(int current, vector<int> &S, float sharpness, int* tag, int length, int maxStringLength) {
		if (length == maxStringLength) {
			return pair<float, vector<int>>(sharpness, S);
		}
		vector<int> bestS;
		float bestSharpness = sharpness;
		vector<int> outNeighbours = getOutNeighbours(he->getTarget(current));
		for (int i=0; i<outNeighbours.size(); i++) {
			int next = outNeighbours[i];
			if (next == he->getOpposite(current) || tag[next] == 2) {
				continue;
			}
			S.push_back(next);
			pair<float, vector<int>> result = DFS(next, S, sharpness + (*EdgeSharpness)(EdgeMap[next]), tag, length+1, maxStringLength);
			S.pop_back();
			if (result.first > bestSharpness) {
				bestSharpness = result.first;
				bestS = result.second;
			}
		}
		return pair<float, vector<int>>(bestSharpness, bestS);
	}

	void tagNeighbours(int v, int* tag) {
		vector<int> neighbours = getInNeighbours(v);
		for (int j=0; j<neighbours.size(); j++) {
			if (tag[neighbours[j]] == 0) {
				tag[neighbours[j]] = 2;
			}
		}
		neighbours = getOutNeighbours(v);
		for (int j=0; j<neighbours.size(); j++) {
			if (tag[neighbours[j]] == 0) {
				tag[neighbours[j]] = 2;
			}
		}
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
				sharpness += getEdgeSharpness(EdgeMap[current]);
				pair<float, vector<int>> result = DFS(current, S, sharpness, tag, 1, maxStringLength);
				S = result.second;
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
			tagNeighbours(start, tag);
			for (int i=0; i<detectedFeature.size(); i++) {
				int current = he->getTarget(detectedFeature[i]);
				tagNeighbours(current, tag);
			}
		}
	}

private:
	igl::opengl::glfw::Viewer *viewer;
	float threshold;
	int* EdgeMap;
	/** Half-edge representation of the original input mesh */
	HalfedgeDS *he;
	MatrixXd *V; // vertex coordinates of the original input mesh
	MatrixXi *F; // REMARK: not needed if using the half-edge data structure
	VectorXd *EdgeSharpness;

	int nEdges;
	//int nVertices, nFaces; // number of vertices, faces
};
