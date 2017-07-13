#pragma once

// Quadrature points in current triangle T (only for (k_q, l_q) = (1,1) and (2,3))
vector<Point> QuadraturePointsInTriang(Triangle T, int k_q, int l_q);
// Value of integral on the triangle T with using quadrature rules 
double QuadrMethod(Triangle T, int k_q, int l_q);
// Change of order of vertices in triangle when it intersect by interface (P_1 is the vertice when the value of level set function has the different sign from the other two vertices)
Triangle RightTriang(Triangle Triang, double time);
// Points of intersection triangle with an interface (good case). Return pair (left_point, right_point) of ordered points.   
vector<Point> IntersectInterfaceTriangle(Triangle Triang, double time);
// Pair of vector triangles (DownTriangles, UpTriangles) obtained by method of connection points on interface with a vertice.
pair<vector<Triangle>, vector<Triangle>> DivisionTriangleWithCurveOneVertice(Triangle Triang, int NumPointsOnCurve, double time);
// Pair of vector triangles (DownTriangles, UpTriangles) obtained by method of projection points on interface on an edge. 
pair<vector<Triangle>, vector<Triangle>> DivisionTriangleWithCurveProjection(Triangle Triang, int NumPointsOnCurve, double time);
// Pair of vector triangles (DownTriangles, UpTriangles) obtained by method of connection points on interface with barycenters.
pair<vector<Triangle>, vector<Triangle>> DivisionTriangleWithCurveBarycenter(Triangle Triang, int NumPointsOnCurve, double time);
// Pair of values of integrals for down and up parts of triangle T divided by moving interface (cases = 1 for method of one vertice, 2 - for projection, 3 - for barycenter)
pair<double, double> QuadrMethodWithInterface(Triangle T, int k_q, int l_q, int N, int cases, double time);
// Write script for MATLAB with one triangle which intersect an interface 
// and the corresponding down and up triangles using different methods of those creation (cases = 1 for method of one vertice, 2 - for projection, 3 - for barycenter)
void DrawTriangle(Triangle T, int NumPointsOnCurve, int cases, double time);

