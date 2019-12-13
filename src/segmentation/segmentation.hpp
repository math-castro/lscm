#pragma once

#include <igl/opengl/glfw/Viewer.h>

void colorSharpEdges();
void colorFeatures();
Vector3d getNormal(int f);
float getEdgeSharpness(int e);
void getEdgeSharpnessMatrix();
void setThreshold(float newThreshold);
float getThreshold();
vector<int> getInNeighbours(int v);
vector<int> getOutNeighbours(int v);
pair<float, vector<int>> DFS(int current, vector<int> &S, float sharpness, int length, int maxStringLength);
void tagNeighbours(int v);
void expandFeatureCurve(int startEdge)
void print(vector<int> v);
