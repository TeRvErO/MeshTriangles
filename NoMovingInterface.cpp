#include "Mesh triangles.h"
#include "NoMovingInterface.h"
double functionForIntegral(Point A) {
	//return pow(A.gety() - 0.1, 2) + pow(A.getx() - 0.1, 2) - 0.05; //circle
	//return A.gety() + A.getx() - 0.5; // straight
	return A.gety()*A.getx(); // parabola
							  //return A.gety() - 4*A.getx() + 4.; // straight
}
// Equation of interface
double EqInterface(Point A) {
	//return pow(A.gety() - 0.1, 2) + pow(A.getx() - 0.1, 2) - 0.05; //circle
	//return A.gety() + A.getx() - 0.5; // straight
	return A.gety() - A.getx()*(1 - A.getx()); // parabola
											   //return A.gety() - 4*A.getx() + 4.; // straight
}
// Simple check of intersection interface with edge
bool InterSimple(Edge E) {
	return (EqInterface(E.getA()) >= 0 && EqInterface(E.getB()) <= 0) || (EqInterface(E.getA()) <= 0 && EqInterface(E.getB()) >= 0);
}
bool CheckIntersection(vector<Edge> AllEdges) {
	for (int i = 0; i < AllEdges.size(); i++) {
		bool Intersect = InterSimple(AllEdges.at(i));
		if (Intersect == 1)
			return 1;
	}
	return 0;
}
bool TriangIntoSurface(Triangle T) {
	return EqInterface(T.getD1().getA()) <= 0 && EqInterface(T.getD1().getB()) <= 0 && EqInterface(T.getD3().getA()) <= 0;
}
// Creation of vector of little edges from initial edge and increment "grid_step"
vector<Edge> DividingEdge(Edge E, double grid_step) {
	vector<Edge> AllEdges;
	vector<Point> PointsOnEdge;
	PointsOnEdge.push_back(E.getA());
	int i = 1;
	do {
		Point current(E.getA().getx() + i*grid_step*(E.getB().getx() - E.getA().getx()) / E.length(), E.getA().gety() + i*grid_step*(E.getB().gety() - E.getA().gety()) / E.length());
		PointsOnEdge.push_back(current);
		AllEdges.push_back(Edge(PointsOnEdge.at(i - 1), PointsOnEdge.at(i)));
		i++;
	} while (Edge(PointsOnEdge.at(i - 1), E.getB()).length() > grid_step);
	AllEdges.push_back(Edge(PointsOnEdge.at(i - 1), E.getB()));
	return AllEdges;
}
// Search of intersections between mesh triangles and interface
void SearchIntersect(vector<Triangle> TriangleMesh, double grid_step) {
	for (int i = 0; i < TriangleMesh.size(); i++) {
		bool intertrian = InterSimple(TriangleMesh.at(i).getD1()) || InterSimple(TriangleMesh.at(i).getD2()) || InterSimple(TriangleMesh.at(i).getD3());
		if (intertrian == 1) {
			cout << i + 1 << ": intersect" << endl;
			TriangleMesh.at(i).ChangeColorTriang('r');
		}
		else {
			bool interwithdiv = CheckIntersection(DividingEdge(TriangleMesh.at(i).getD1(), grid_step)) || CheckIntersection(DividingEdge(TriangleMesh.at(i).getD2(), grid_step)) ||
				CheckIntersection(DividingEdge(TriangleMesh.at(i).getD3(), grid_step));
			bool into = TriangIntoSurface(TriangleMesh.at(i));
			if (interwithdiv == 1) {
				cout << i + 1 << ": found intersection with dividing algorithm" << endl;
				TriangleMesh.at(i).ChangeColorTriang('r');
			}
			else
				if (into == 1)
					cout << i + 1 << ": into" << endl;
				else cout << i + 1 << ": extra" << endl;
		}
	}
}
void DrawEdgesMatlab(vector<Edge> EdgesMesh) {
	ofstream fout;
	fout.open("DrawingEdges.m");
	for (int i = 0; i < EdgesMesh.size(); i++) {
		fout << "drawLine([" << EdgesMesh.at(i).getA().getx() << "," << EdgesMesh.at(i).getA().gety() << "], " << "[" << EdgesMesh.at(i).getB().getx() << "," << EdgesMesh.at(i).getB().gety() << "], '" << EdgesMesh.at(i).getcolour() << "');" << endl;
		fout << "hold on" << endl;
	}
	fout << "x = 0:0.01:1;" << endl;
	fout << "y = x.*(1-x);" << endl;
	fout << "plot(x, y, 'r')" << endl;
	fout.close();
}