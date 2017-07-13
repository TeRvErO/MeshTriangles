#pragma once

// Equation of moving interface (let "time" = 0 for no move)
double EqInterfaceAnim(Point A, double time);
// Explicite equation of x for zero level of moving interface (let "time" = 0 for no move)
double EqInterfaceAnimExpliciteY(double x, double time);
// Explicite equation of y for zero level of moving interface (let "time" = 0 for no move)
vector<double> EqInterfaceAnimExpliciteX(double y, double time);

// Simple check of intersection edge with moving interface (let "time" = 0 for no move)
bool InterSimpleAnim(Edge E, double time);
// Simple check of intersection vector of edges with moving interface (let "time" = 0 for no move)
bool CheckIntersectionAnim(vector<Edge> AllEdges, double time);

// Discrete distance of Hausdorf between two curves with given points
double DiscrHausdDist(vector<Point> firstcurve, vector<Point> secondcurve);
Point getPointInThisEdge(Edge E, Point A, Point B);
vector<Point> PointsOnInterfaceInRectang(Rectangular R, int Num_points, double time);

vector<Point> PointsOnCurveIntoTriang(Triangle T, int N, double time);
vector<Point> PointsOnCurve(Point Start, Point Finish, int N, double time);
vector<Point> PointsOnEdge(Edge E, int N);
void SearchIntersectEdgesAnim(vector<Edge> &EdgesMesh, double grid_step, double time);
void SearchIntersectTriangAnim(vector<Triangle> &TrianglesMesh, double grid_step, double time);
// Add edges from current vector of triangles T which are not yet in mesh of edges
void AddEdgesFromTriangInMeshEdges(vector<Triangle> TrianglesMesh, vector<Edge> &EdgesMesh);
// Write script code for MATLAB with mesh of triangles and moving interface
void DrawTrianglesMatlabAnimation(vector<Triangle> TrianglesMesh, double grid_step);
// Write script code for MATLAB with mesh of edges and moving interface
void DrawEdgesMatlabAnimation(vector<Edge> EdgesMesh, double grid_step);