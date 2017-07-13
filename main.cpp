#include "Mesh triangles.h"
#include "Efficient integration.h"
#include "NoMovingInterface.h"
#include "MovingInterface.h"
#include "IntegrationQuadratures.h"

int main() {
	//real mesh data
	const int Size = 40000;
	float data[Size];
	ifstream input("tri01.mesh2d");

	for (int i = 0; i < Size; i++)
		input >> data[i];
	int NumPoints = data[0];
	vector<Point> PointsMesh(NumPoints);
	for (int i = 0; i < PointsMesh.size(); i++) {
		PointsMesh.at(i) = Point(data[2 * i + 1], data[2 * i + 2]);
		//cout << i+1 << ": "<< PointsMesh.at(i) << endl;
	}
	int NumTriang = data[2 * NumPoints + 1];
	vector<Triangle> TrianglesMesh(NumTriang);
	vector<Edge> AllEdgesMesh(3 * NumTriang);
	//cout<<NumTriang;
	for (int i = 0; i < TrianglesMesh.size(); i++) {
		TrianglesMesh.at(i) = Triangle(PointsMesh.at(data[2 * NumPoints + 4 * i + 3] - 1), PointsMesh.at(data[2 * NumPoints + 4 * i + 4] - 1), PointsMesh.at(data[2 * NumPoints + 4 * i + 5] - 1));
		AllEdgesMesh.at(3 * i) = TrianglesMesh.at(i).getD1();
		AllEdgesMesh.at(3 * i + 1) = TrianglesMesh.at(i).getD2();
		AllEdgesMesh.at(3 * i + 2) = TrianglesMesh.at(i).getD3();
		//cout << i+1 << ": "<< TrianglesMesh.at(i) << endl;
	}
	int NumEdgeBorder = data[2 * NumPoints + 4 * NumTriang + 2];
	//cout<<NumEdgeBorder;
	vector<Edge> EdgeBorderMesh(NumEdgeBorder);
	for (int i = 0; i < EdgeBorderMesh.size(); i++) {
		EdgeBorderMesh.at(i) = Edge(PointsMesh.at(data[2 * NumPoints + 4 * NumTriang + 3 * i + 4] - 1), PointsMesh.at(data[2 * NumPoints + 4 * NumTriang + 3 * i + 5] - 1));
		//cout << i+1 << ": "<< EdgeBorderMesh.at(i) << endl;
	}

	double h = 0.05;

	//SearchIntersect(TrianglesMesh, h);
	//SearchIntersectEdgesAnim(AllEdgesMesh, h);
	//DrawEdgesMatlab(AllEdgesMesh);

	int N = 10;
	//cout << DiscrHausdDist(PointsOnCurve(edge.getA(), edge.getB(), N), PointsOnEdge(edge, N));
	//cout << DiscrHausdDist(PointsOnEdge(edge1, N), PointsOnEdge(edge2, N));
	//DrawEdgesMatlabAnimation(AllEdgesMesh, h);

	Triangle Tik(Point(0.1, 0.5), Point(0.9, 0.5), Point(0.5, 0.9));
	Triangle Tak(Point(0.1, 0.5), Point(0.1, 0.9), Point(0.5, 0.9));
	Triangle Tuk(Point(0.9, 0.9), Point(0.9, 0.5), Point(0.5, 0.9));
	vector<Triangle> Ttr;
	Ttr.push_back(Tik);
	Ttr.push_back(Tak);
	Ttr.push_back(Tuk);
	//DrawTrianglesMatlabAnimation(Ttr, h);

	//DrawTrianglesMatlabAnimation(TrianglesMesh, h);
	Triangle T(Point(0, 0.), Point(0.5, 1.), Point(1., 0));
	//cout << T.barycenter() << " "<< T.square() << endl;
	//cout << QuadrMethod(T, 2,3) <<endl;
	double down = 0.0371813;
	cout << "real values of down integral = " << down << " and up integral = " << QuadrMethod(T, 2, 3) - down << "; sum of it = " << QuadrMethod(T, 2, 3) << endl;
	int NumberPoints = 20;
	int cases;
	double time = 0.24;
	pair<double, double> DownAndUpIntegrals;
	DownAndUpIntegrals = QuadrMethodWithInterface(T, 2, 3, NumberPoints, cases = 1, time);
	cout << "First case: 'Number' of points = " << NumberPoints << "; value  of down integral = " << DownAndUpIntegrals.first << "; up integral = " << DownAndUpIntegrals.second << "; sum of it = " << DownAndUpIntegrals.first + DownAndUpIntegrals.second << endl;
	DownAndUpIntegrals = QuadrMethodWithInterface(T, 2, 3, NumberPoints, cases = 2, time);
	cout << "Second case (with projection): 'Number' of points = " << NumberPoints << "; value  of down integral = " << DownAndUpIntegrals.first << "; up integral = " << DownAndUpIntegrals.second << "; sum of it = " << DownAndUpIntegrals.first + DownAndUpIntegrals.second << endl;
	DownAndUpIntegrals = QuadrMethodWithInterface(T, 2, 3, NumberPoints, cases = 3, time);
	cout << "Third case (with barycenters): 'Number' of points = " << NumberPoints << "; value  of down integral = " << DownAndUpIntegrals.first << "; up integral = " << DownAndUpIntegrals.second << "; sum of it = " << DownAndUpIntegrals.first + DownAndUpIntegrals.second << endl;

	//Draw different cases for division on triangles with points on interface
	NumberPoints = 10;
	// Draw triangles with the first method naive
	DrawTriangle(T, NumberPoints, cases = 1, time);

	// Draw triangles with the second method with projections of points 
	DrawTriangle(T, NumberPoints, cases = 2, time);

	// Draw triangles with the third method with two barycenters
	DrawTriangle(T, NumberPoints, cases = 3, time);

	//efficient integration
	Triangle Tr(Point(0, 2), Point(2, 0), Point(-2, 0));
	PrintSituation(Tr);
	cout << "Value of integral = " << EffIntegration(Tr) << endl;
	getch();
}