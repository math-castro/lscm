/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

/**
 * @author Luca Castelli Aleardi (2019)
 * Minimal array-based implementation of a the Half-edge data structure for representing polygonal meshes<br>
 * Features:<br>
 * -) References are encoded with integer values (32 bits), stores in an array of int <br>
 * -) vertex coordinates are not stored <br>
 * -) no implementation of dynamic updates <br>
 */
class HalfedgeDS
{

public:
	/** 
	 * Allocate the memory space for storing the data structure.
	 * This version of the constructor does NOT allocate face indices (for reducing memory requirements)
	 **/
	HalfedgeDS(int n, int h)
	{
		nVertices = n;
		nHalfedges = h;
		T = new int[nHalfedges * sizeT];
		for (int i = 0; i < nHalfedges * sizeT; i++)
			T[i] = -1; // "-1" means that reference is NOT DEFINED (null)

		incidentEdge = new int[nVertices];
	}

	/** 
	 * Allocate the memory space for storing the data structure.
	 * This version of the constructor stores face/halfedge incidence relations
	 **/
	HalfedgeDS(int n, int h, int f)
	{
		nVertices = n;
		nHalfedges = h;
		nFaces = f;
		T = new int[nHalfedges * sizeT];
		for (int i = 0; i < nHalfedges * sizeT; i++)
			T[i] = -1; // "-1" means that reference is NOT DEFINED (null)

		incidentEdge = new int[nVertices];
		faces=new int[nFaces];
	}

	/** 
	 * Set the opposite half-edge 
	 **/
	void setOpposite(int e, int eOpposite)
	{
		T[e * sizeT] = eOpposite;
	}

	/** 
	 * Set the next half-edge: the half-edge 'eNext' following 'e' in the same face 
	 **/
	void setNext(int e, int eNext)
	{
		T[sizeT * e + 1] = eNext;
	}

	/** 
	 * Set the (target) incident vertex of 'e' 
	 **/
	void setVertex(int e, int v)
	{
		T[sizeT * e + 2] = v;
	}

	/** 
	 * Set the face containing the given half-edge
	 **/
	void setFace(int e, int f)
	{
		T[sizeT * e + 3] = f;
	}

	/** 
	 * Set the face containing the given half-edge
	 **/
	void setPrev(int e, int ePrev)
	{
		T[sizeT * e + 4] = ePrev;
	}

	void setEdge(int v, int e)
	{
		incidentEdge[v] = e;
	}

	/** 
	 * Set the half-edge 'e' incident to the given face 'f'
	 **/
	int setEdgeInFace(int f, int e)
	{
		faces[f]=e;
	}

	//--- methods for accessing data and navigating between mesh elements ---

	/** 
	 * Return the following edge in ccw orientation around the same face
	 **/
	int getNext(int e)
	{
		return T[sizeT * e + 1];
	}

	/** 
	 * Return the opposite edge in the neighboring face (having opposite direction)
	 **/
	int getOpposite(int e)
	{
		return T[e * sizeT];
	}

	/** 
	 * Return the previous edge in ccw orientation around the same face
	 **/
	int getPrev(int e)
	{
		return T[sizeT * e + 4];
	}

	/** 
	 * Return the target vertex incident to the half-edge (its target vertex)
	 **/
	int getTarget(int e)
	{
		return T[sizeT * e + 2];
	}

	/** 
	 * Return the face containing the half-edge
	 **/
	int getFace(int e)
	{
		return T[sizeT * e + 3];
	}

	/** 
	 * Return a half edge incident to the vertex
	 **/
	int getEdge(int v)
	{
		return incidentEdge[v];
	}

	/** 
	 * Return a half edge incident to the face
	 **/
	int getEdgeInFace(int f)
	{
		return faces[f];
	}

	/** 
	 * Return the number of vertices
	 **/
	int sizeOfVertices()
	{
		return nVertices;
	}

	/** 
	 * Return the number of half-edges
	 **/
	int sizeOfHalfedges()
	{
		return nHalfedges;
	}

	/** 
	 * Return the number of faces
	 **/
	int sizeOfFaces()
	{
		return nFaces;
	}

	/** 
	 * Print the array T[] storing all references
	 **/
	void print()
	{
		for (int i = 0; i < nHalfedges; i++)
		{
			cout << "he" << i << ": \t" << T[sizeT * i] << "\t" << T[sizeT * i + 1] << "\t" << T[sizeT * i + 2] << "\t" << T[sizeT * i + 3] << endl;
		}

		cout << "face list: " << nFaces << endl;
		if(faces!=NULL) {
		for (int i = 0; i < nFaces; i++)
		{
			cout << "f" << i << ": \t";
			int e1=getEdgeInFace(i);
			cout << "incident edge: e" << e1;
			int e2=getNext(e1);
			int e3=getNext(e2);
			cout << "\t v" << getTarget(e3) << ", v" << getTarget(e1) << ", v" << getTarget(e2);
			cout << endl;
		}

		}
	}

private:
	int nVertices, nHalfedges, nFaces; // number of vertices, halfedges and faces in the mesh

	int *T;			   // a table for storing references between half-edges: each halfedge is represented with three integer references
	int *incidentEdge; // for each vertex we store an incident (ingoing) halfedge

	int *faces; // for each face we store an incident halfedge

	const int sizeT = 5;
};
