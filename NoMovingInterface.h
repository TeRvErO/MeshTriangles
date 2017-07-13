//#include "Mesh triangles.h"

// Integrating fuction
double functionForIntegral(Point A);
// Equation of interface
double EqInterface(Point A);
// Simple check of intersection interface with edge
bool InterSimple(Edge E);
// Simple check of intersection for vector of edges
bool CheckIntersection(vector<Edge> AllEdges);
bool TriangIntoSurface(Triangle T);
// Creation of vector of little edges from initial edge and increment "grid_step"
vector<Edge> DividingEdge(Edge E, double grid_step);
// Search of intersections between mesh triangles and interface
void SearchIntersect(vector<Triangle> TriangleMesh, double grid_step);

void DrawEdgesMatlab(vector<Edge> EdgesMesh);