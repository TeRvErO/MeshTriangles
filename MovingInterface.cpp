#include "Mesh triangles.h"
#include "NoMovingInterface.h"
#include "MovingInterface.h"

//equation of moving interface
double EqInterfaceAnim(Point A, double time) {
	//return pow(A.gety() - 0.1 - time, 2) + pow(A.getx() - 0.1 - time, 2) - 0.05; //circle
	//return A.gety() + A.getx() - 0.5 - time; // straight
	return A.gety() - (A.getx() + 0.5)*(0.5 - A.getx()) - time; // parabola
}
double EqInterfaceAnimExpliciteY(double x, double time) {
	//return -x + 0.5; // straight
	return (x + 0.5)*(0.5 - x) + time; // parabola
}
vector<double> EqInterfaceAnimExpliciteX(double y, double time) {
	vector<double> x;
	//x.push_back(-y + 0.5 - time) // straight 
	x.push_back(1. + sqrt(1 - 4 * (y - time)) / 2.); x.push_back(1. - sqrt(1 - 4 * (y - time)) / 2.); // parabola																				
	return x; // parabola

}

bool InterSimpleAnim(Edge E, double time) {
	return (EqInterfaceAnim(E.getA(), time) >= 0 && EqInterfaceAnim(E.getB(), time) <= 0) || (EqInterfaceAnim(E.getA(), time) <= 0 && EqInterfaceAnim(E.getB(), time) >= 0);
}
bool CheckIntersectionAnim(vector<Edge> AllEdges, double time) {
	for (int i = 0; i < AllEdges.size(); i++) {
		bool Intersect = InterSimpleAnim(AllEdges.at(i), time);
		if (Intersect == 1)
			return 1;
	}
	return 0;
}

double DiscrHausdDist(vector<Point> firstcurve, vector<Point> secondcurve) {
	double haus_dist = 0;
	for (int i = 0; i < firstcurve.size(); i++) {
		double shortest = std::numeric_limits<double>::infinity();
		for (int j = 0; j < secondcurve.size(); j++) {
			double dij = firstcurve.at(i).dist(secondcurve.at(j));
			if (dij < shortest)
				shortest = dij;
		}
		if (shortest > haus_dist)
			haus_dist = shortest;
	}
	for (int j = 0; j < secondcurve.size(); j++) {
		double shortest = std::numeric_limits<double>::infinity();
		for (int i = 0; i < firstcurve.size(); i++) {
			double dij = secondcurve.at(j).dist(firstcurve.at(i));
			if (dij < shortest)
				shortest = dij;
		}
		if (shortest > haus_dist)
			haus_dist = shortest;
	}
	return haus_dist;
}
Point getPointInThisEdge(Edge E, Point A, Point B) {
	if (E.In(A) == 1)
		return A;
	else if (E.In(B) == 1)
		return B;
	else {
		Point P(-1000., -1000.); //fictive
		return P;
	}
}
vector<Point> PointsOnInterfaceInRectang(Rectangular R, int Num_points, double time) {
	Point Startx = R.getPoints()[0];
	Point Finishx = R.getPoints()[3];
	Point Starty = R.getPoints()[0];
	Point Finishy = R.getPoints()[1];

	vector<Point> Points;
	for (int i = 0; i < Num_points + 1; i++) {
		double next_x = Startx.getx() + i*(Finishx.getx() - Startx.getx()) / Num_points;
		double next_y = EqInterfaceAnimExpliciteY(next_x, time);
		if (next_y >= Starty.gety() && next_y <= Finishy.gety())
			Points.push_back(Point(next_x, next_y));
	}
	/*for (int i = Points.size()/3; i < Points.size(); ++(++i)){
	double next_y = Starty.gety() + i*h*(Finishy.gety() - Starty.gety())/Starty.dist(Finishy);
	Points.at(i) = Point(EqInterfaceAnimExpliciteX(next_y, time).at(0), next_y);
	Points.at(i + 1) = Point(EqInterfaceAnimExpliciteX(next_y, time).at(1), next_y);
	}*/
	return Points;
}

vector<Point> PointsOnCurveIntoTriang(Triangle T, int N, double time) {
	Rectangular R(T);
	vector<Point> PointsInRect = PointsOnInterfaceInRectang(R, N, time);
	vector<Point> PointsInTriang;
	for (int i = 0; i < PointsInRect.size(); i++)
		if (T.pointInTriangle(PointsInRect.at(i)))
			PointsInTriang.push_back(PointsInRect.at(i));
	return PointsInTriang;
}
vector<Point> PointsOnCurve(Point Start, Point Finish, int N, double time) {
	double h = Start.dist(Finish) / N;
	vector<Point> Points(3 * (N + 1));
	for (int i = 0; i < Points.size() / 3; i++) {
		double next_x = Start.getx() + i*h*(Finish.getx() - Start.getx()) / Start.dist(Finish);
		Points.at(i) = Point(next_x, EqInterfaceAnimExpliciteY(next_x, time));
	}
	for (int i = Points.size() / 3; i < Points.size(); ++(++i)) {
		double next_y = Start.gety() + i*h*(Finish.gety() - Start.gety()) / Start.dist(Finish);
		Points.at(i) = Point(EqInterfaceAnimExpliciteX(next_y, time).at(0), next_y);
		Points.at(i + 1) = Point(EqInterfaceAnimExpliciteX(next_y, time).at(1), next_y);
	}
	return Points;
}
vector<Point> PointsOnEdge(Edge E, int N) {
	vector<Point> Points(N + 1);
	for (int i = 0; i < Points.size(); i++)
		Points.at(i) = Point(E.getA().getx() + i*(E.getB().getx() - E.getA().getx()) / N, E.getA().gety() + i*(E.getB().gety() - E.getA().gety()) / N);
	return Points;
}
void SearchIntersectEdgesAnim(vector<Edge> &EdgesMesh, double grid_step, double time) {
	for (int i = 0; i < EdgesMesh.size(); i++) {
		bool intertrian = InterSimpleAnim(EdgesMesh.at(i), time);
		if (intertrian == 1) {
			int N = 10;
			double h = DiscrHausdDist(PointsOnCurve(EdgesMesh.at(i).getA(), EdgesMesh.at(i).getB(), N, time), PointsOnEdge(EdgesMesh.at(i), N));
			if (h < 0.05) {
				EdgesMesh.at(i).ChangeColorEdge('r');
			}
			else {
				//cout << i+1 << ": intersect" << endl;
				EdgesMesh.at(i).ChangeColorEdge('y');
			}
		}
		else {
			bool interwithdiv = CheckIntersectionAnim(DividingEdge(EdgesMesh.at(i), grid_step), time);
			if (interwithdiv == 1) {
				int N = 100;
				double h = DiscrHausdDist(PointsOnCurve(EdgesMesh.at(i).getA(), EdgesMesh.at(i).getB(), N, time), PointsOnEdge(EdgesMesh.at(i), N));
				if (h < 0.05)
					EdgesMesh.at(i).ChangeColorEdge('r');
				else {
					//cout << i+1 << ": found intersection with dividing algorithm" << endl;
					EdgesMesh.at(i).ChangeColorEdge('y');
				}
			}
		}
	}
}
void SearchIntersectTriangAnim(vector<Triangle> &TrianglesMesh, double grid_step, double time) {
	for (int i = 0; i < TrianglesMesh.size(); i++) {
		bool intertrian = InterSimpleAnim(TrianglesMesh.at(i).getD1(), time) || InterSimpleAnim(TrianglesMesh.at(i).getD2(), time) || InterSimpleAnim(TrianglesMesh.at(i).getD3(), time);
		if (intertrian == 1) {
			int N = 50;
			vector<Point> PointsOnCurve = PointsOnCurveIntoTriang(TrianglesMesh.at(i), N, time);
			double h1 = DiscrHausdDist(PointsOnCurve, PointsOnEdge(TrianglesMesh.at(i).getD1(), N));
			double h2 = DiscrHausdDist(PointsOnCurve, PointsOnEdge(TrianglesMesh.at(i).getD2(), N));
			double h3 = DiscrHausdDist(PointsOnCurve, PointsOnEdge(TrianglesMesh.at(i).getD3(), N));
			if (min(h1, min(h2, h3)) < 0.05) {
				//cout << i+1 << ": bad intersect " << endl;
				TrianglesMesh.at(i).ChangeColorTriang('r');
			}
			else {
				//cout << i+1 << ": good intersect" << endl;
				TrianglesMesh.at(i).ChangeColorTriang('y');
			}
		}
		else {
			bool interwithdiv = CheckIntersectionAnim(DividingEdge(TrianglesMesh.at(i).getD1(), grid_step), time) || CheckIntersectionAnim(DividingEdge(TrianglesMesh.at(i).getD2(), grid_step), time) ||
				CheckIntersectionAnim(DividingEdge(TrianglesMesh.at(i).getD3(), grid_step), time);
			if (interwithdiv == 1) {
				//cout << i+1 << ": bad intersect " << endl;
				TrianglesMesh.at(i).ChangeColorTriang('r');
			}
		}
	}
}
void AddEdgesFromTriangInMeshEdges(vector<Triangle> TrianglesMesh, vector<Edge> &EdgesMesh) {
	for (int i = 0; i < TrianglesMesh.size(); i++) {
		bool AlreadyInEdgesD1 = false;
		bool AlreadyInEdgesD2 = false;
		bool AlreadyInEdgesD3 = false;
		for (int j = 0; j < EdgesMesh.size(); j++) {
			if (TrianglesMesh.at(i).getD1() == EdgesMesh.at(j))
				AlreadyInEdgesD1 = true;
			if (TrianglesMesh.at(i).getD2() == EdgesMesh.at(j))
				AlreadyInEdgesD2 = true;
			if (TrianglesMesh.at(i).getD3() == EdgesMesh.at(j))
				AlreadyInEdgesD3 = true;
		}
		if (AlreadyInEdgesD1 == false)
			EdgesMesh.push_back(TrianglesMesh.at(i).getD1());
		if (AlreadyInEdgesD2 == false)
			EdgesMesh.push_back(TrianglesMesh.at(i).getD2());
		if (AlreadyInEdgesD3 == false)
			EdgesMesh.push_back(TrianglesMesh.at(i).getD3());
	}
}
void DrawTrianglesMatlabAnimation(vector<Triangle> TrianglesMesh, double grid_step) {
	const int N = 100;
	double t[N];
	vector<Edge> AllEdges(1);
	vector<Triangle> RedTriang(1);
	vector<Triangle> YellowTriang(1);
	vector<Triangle> BlueTriang(1);
	ofstream fout;
	fout.open("DrawingEvolutionTriangles.m");
	for (int i = 0; i < N; i++) {
		vector<Triangle> CurrentMesh(TrianglesMesh);
		t[i] = i*1. / N;
		SearchIntersectTriangAnim(CurrentMesh, grid_step, t[i]);

		/*for(int iter = 0; iter < CurrentMesh.size(); iter++){
		AllEdges.push_back(CurrentMesh.at(iter).getD1());
		AllEdges.push_back(CurrentMesh.at(iter).getD2());
		AllEdges.push_back(CurrentMesh.at(iter).getD3());
		}*/

		for (int iter = 0; iter < CurrentMesh.size(); iter++) {
			if (CurrentMesh.at(iter).getcolour() == 'r')
				RedTriang.push_back(CurrentMesh.at(iter));
			else
				if (CurrentMesh.at(iter).getcolour() == 'y')
					YellowTriang.push_back(CurrentMesh.at(iter));
				else BlueTriang.push_back(CurrentMesh.at(iter));
		}
		AddEdgesFromTriangInMeshEdges(BlueTriang, AllEdges);

		AddEdgesFromTriangInMeshEdges(YellowTriang, AllEdges);

		AddEdgesFromTriangInMeshEdges(RedTriang, AllEdges);

		fout << "clf()" << endl;
		for (int j = 0; j < AllEdges.size(); j++) {
			fout << "drawLine([" << AllEdges.at(j).getA().getx() << "," << AllEdges.at(j).getA().gety() << "], " << "[" << AllEdges.at(j).getB().getx() << ","
				<< AllEdges.at(j).getB().gety() << "], '" << AllEdges.at(j).getcolour() << "');" << endl;
			fout << "hold on" << endl;
		}
		fout << "x = 0:0.01:1;" << endl;
		fout << "y = x.*(1-x) + " << t[i] << ";" << endl;
		fout << "plot(x, y, 'r');" << endl;
		fout << "axis([0 1 0 1.2]);" << endl;
		fout << "pause(0.5)" << endl;
		RedTriang.clear();
		YellowTriang.clear();
		BlueTriang.clear();
		AllEdges.clear();
		CurrentMesh = TrianglesMesh;
	}
	fout.close();
}
void DrawEdgesMatlabAnimation(vector<Edge> EdgesMesh, double grid_step) {
	const int N = 100;
	double t[N];
	vector<Edge> CurrentMesh(EdgesMesh);
	ofstream fout;
	fout.open("DrawingEvolutionEdges.m");
	for (int j = 0; j < N; j++) {
		t[j] = j*1. / N;
		SearchIntersectEdgesAnim(CurrentMesh, grid_step, t[j]);
		fout << "clf()" << endl;
		for (int i = 0; i < CurrentMesh.size(); i++) {
			fout << "drawLine([" << CurrentMesh.at(i).getA().getx() << "," << CurrentMesh.at(i).getA().gety() << "], " << "[" << CurrentMesh.at(i).getB().getx() << ","
				<< CurrentMesh.at(i).getB().gety() << "], '" << CurrentMesh.at(i).getcolour() << "');" << endl;
			fout << "hold on" << endl;
		}
		fout << "x = 0:0.01:1;" << endl;
		fout << "y = x.*(1-x) + " << t[j] << ";" << endl;
		fout << "plot(x, y, 'r');" << endl;
		fout << "axis([0 1 0 1.2]);" << endl;
		fout << "pause(0.5)" << endl;
		CurrentMesh.clear();
		CurrentMesh = EdgesMesh;
	}
	fout.close();
}