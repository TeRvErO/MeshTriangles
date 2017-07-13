#include "Mesh triangles.h"
#include "NoMovingInterface.h"
#include "MovingInterface.h"
#include "IntegrationQuadratures.h"

vector<Point> QuadraturePointsInTriang(Triangle T, int k_q, int l_q) {
	vector<Point> PointsForQuadratures;
	if (k_q == 1 && l_q == 1)
		PointsForQuadratures.push_back(T.barycenter());
	if (k_q == 2 && l_q == 3) {
		PointsForQuadratures.push_back(Point(1. / 2.*T.getP1().getx() + 1. / 2.*T.getP2().getx() + 0 * T.getP3().getx(), 1. / 2.*T.getP1().gety() + 1. / 2.*T.getP2().gety() + 0 * T.getP3().gety()));
		PointsForQuadratures.push_back(Point(1. / 2.*T.getP1().getx() + 0 * T.getP2().getx() + 1. / 2.*T.getP3().getx(), 1. / 2.*T.getP1().gety() + 0 * T.getP2().gety() + 1. / 2.*T.getP3().gety()));
		PointsForQuadratures.push_back(Point(0 * T.getP1().getx() + 1. / 2.*T.getP2().getx() + 1. / 2.*T.getP3().getx(), 0 * T.getP1().gety() + 1. / 2.*T.getP2().gety() + 1. / 2.*T.getP3().gety()));
	}
	return PointsForQuadratures;
}
double QuadrMethod(Triangle T, int k_q, int l_q) {
	vector<Point> PointsForQuadratures = QuadraturePointsInTriang(T, k_q, l_q);
	if (k_q == 1 && l_q == 1)
		return T.square()*functionForIntegral(PointsForQuadratures.at(0));
	if (k_q == 2 && l_q == 3) {
		double result = 0;
		for (int i = 0; i < PointsForQuadratures.size(); i++)
			result += T.square()*functionForIntegral(PointsForQuadratures.at(i)) / 3.;
		return result;
	}
}
Triangle RightTriang(Triangle Triang, double time) {
	vector<Point> Positive;
	vector<Point> Negative;
	Triangle T;
	if (EqInterfaceAnim(Triang.getP1(), time) <= 0)
		Negative.push_back(Triang.getP1());
	else Positive.push_back(Triang.getP1());
	if (EqInterfaceAnim(Triang.getP2(), time) <= 0)
		Negative.push_back(Triang.getP2());
	else Positive.push_back(Triang.getP2());
	if (EqInterfaceAnim(Triang.getP3(), time) <= 0)
		Negative.push_back(Triang.getP3());
	else Positive.push_back(Triang.getP3());

	if (Negative.size() == 1)
		T = Triangle(Negative.at(0), Positive.at(0), Positive.at(1));
	else T = Triangle(Positive.at(0), Negative.at(0), Negative.at(1));
	return T;
}
vector<Point> IntersectInterfaceTriangle(Triangle Triang, double time) {
	vector<Point> Points;
	int Number = 10;
	Triangle T = RightTriang(Triang, time);
	vector<Point> PointsOn12;
	vector<Point> PointsOn13;
	for (int i = 0; i < Number + 1; i++) {
		PointsOn12.push_back(Point(T.getP1().getx() + i*(T.getP2().getx() - T.getP1().getx()) / Number, T.getP1().gety() + i*(T.getP2().gety() - T.getP1().gety()) / Number));
		PointsOn13.push_back(Point(T.getP1().getx() + i*(T.getP3().getx() - T.getP1().getx()) / Number, T.getP1().gety() + i*(T.getP3().gety() - T.getP1().gety()) / Number));
	}
	pair<Point, Point> ChangeSignOn12;
	pair<Point, Point> ChangeSignOn13;
	for (int i = 0; i < Number; i++)
		if ((EqInterfaceAnim(PointsOn12.at(i), time) >= 0 && EqInterfaceAnim(PointsOn12.at(i + 1), time) <= 0) || (EqInterfaceAnim(PointsOn12.at(i), time) <= 0 && EqInterfaceAnim(PointsOn12.at(i + 1), time) >= 0)) {
			ChangeSignOn12 = make_pair(PointsOn12.at(i), PointsOn12.at(i + 1));
			break;
		}
	for (int i = 0; i < Number; i++)
		if ((EqInterfaceAnim(PointsOn13.at(i), time) >= 0 && EqInterfaceAnim(PointsOn13.at(i + 1), time) <= 0) || (EqInterfaceAnim(PointsOn13.at(i), time) <= 0 && EqInterfaceAnim(PointsOn13.at(i + 1), time) >= 0)) {
			ChangeSignOn13 = make_pair(PointsOn13.at(i), PointsOn13.at(i + 1));
			break;
		}
	double x_first = ChangeSignOn12.first.getx() - EqInterfaceAnim(ChangeSignOn12.first, time)*(ChangeSignOn12.second.getx() - ChangeSignOn12.first.getx()) / (EqInterfaceAnim(ChangeSignOn12.second, time) - EqInterfaceAnim(ChangeSignOn12.first, time));
	double y_first = ChangeSignOn12.first.gety() - EqInterfaceAnim(ChangeSignOn12.first, time)*(ChangeSignOn12.second.gety() - ChangeSignOn12.first.gety()) / (EqInterfaceAnim(ChangeSignOn12.second, time) - EqInterfaceAnim(ChangeSignOn12.first, time));
	double x_last = ChangeSignOn13.first.getx() - EqInterfaceAnim(ChangeSignOn13.first, time)*(ChangeSignOn13.second.getx() - ChangeSignOn13.first.getx()) / (EqInterfaceAnim(ChangeSignOn13.second, time) - EqInterfaceAnim(ChangeSignOn13.first, time));
	double y_last = ChangeSignOn13.first.gety() - EqInterfaceAnim(ChangeSignOn13.first, time)*(ChangeSignOn13.second.gety() - ChangeSignOn13.first.gety()) / (EqInterfaceAnim(ChangeSignOn13.second, time) - EqInterfaceAnim(ChangeSignOn13.first, time));
	Point First(x_first, y_first);
	Point Last(x_last, y_last);
	if (First < Last) {
		Points.push_back(First);
		Points.push_back(Last);
	}
	else {
		Points.push_back(Last);
		Points.push_back(First);
	}
	return Points;
}
pair<vector<Triangle>, vector<Triangle>> DivisionTriangleWithCurveProjection(Triangle Triang, int NumPointsOnCurve, double time) {
	vector<Triangle> DownTriangles;
	vector<Triangle> UpTriangles;
	vector<Point> Points;
	vector<Point> Points_Add;
	vector<Point> Points_Proj;
	Triangle T = RightTriang(Triang, time);

	Points_Add = PointsOnCurveIntoTriang(T, NumPointsOnCurve, time);

	vector<Point> IntersectInterface = IntersectInterfaceTriangle(T, time);
	Points.push_back(IntersectInterface.at(0));
	for (int i = 0; i < Points_Add.size(); i++)
		Points.push_back(Points_Add.at(i));
	Points.push_back(IntersectInterface.at(1));

	for (int i = 0; i < Points.size(); i++)
		Points_Proj.push_back(T.getD2().Projection(Points.at(i)));

	Point PBetween12 = getPointInThisEdge(T.getD1(), Points.at(0), Points.at(Points.size() - 1));
	Point PBetween13 = getPointInThisEdge(T.getD3(), Points.at(0), Points.at(Points.size() - 1));

	if (T.getP2().dist(Points_Proj.at(0)) < T.getP2().dist(Points_Proj.at(Points_Proj.size() - 1)))
		DownTriangles.push_back(Triangle(T.getP2(), Points_Proj.at(0), PBetween12));
	else DownTriangles.push_back(Triangle(T.getP2(), Points_Proj.at(Points_Proj.size() - 1), PBetween12));
	for (int i = 0; i < Points.size() - 1; i++) {
		DownTriangles.push_back(Triangle(Points.at(i), Points.at(i + 1), Points_Proj.at(i)));
		DownTriangles.push_back(Triangle(Points_Proj.at(i), Points_Proj.at(i + 1), Points.at(i + 1)));
		UpTriangles.push_back(Triangle(Points.at(i), Points.at(i + 1), T.getP1()));
	}
	if (T.getP3().dist(Points_Proj.at(0)) < T.getP3().dist(Points_Proj.at(Points_Proj.size() - 1)))
		DownTriangles.push_back(Triangle(Points_Proj.at(0), T.getP3(), PBetween13));
	else DownTriangles.push_back(Triangle(Points_Proj.at(Points_Proj.size() - 1), T.getP3(), PBetween13));
	return make_pair(DownTriangles, UpTriangles);
}
pair<vector<Triangle>, vector<Triangle>> DivisionTriangleWithCurveOneVertice(Triangle Triang, int NumPointsOnCurve, double time) {
	vector<Triangle> UpTriangles;
	vector<Triangle> DownTriangles;
	vector<Point> Points;
	vector<Point> Points_Add;
	Triangle T = RightTriang(Triang, time);
	Points_Add = PointsOnCurveIntoTriang(T, NumPointsOnCurve, time);

	vector<Point> IntersectInterface = IntersectInterfaceTriangle(T, time);

	Points.push_back(IntersectInterface.at(0));
	for (int i = 0; i < Points_Add.size(); i++)
		Points.push_back(Points_Add.at(i));
	Points.push_back(IntersectInterface.at(1));

	Point PBetween13 = getPointInThisEdge(T.getD3(), Points.at(0), Points.at(Points.size() - 1));
	for (int i = 0; i < Points.size() - 1; i++) {
		DownTriangles.push_back(Triangle(Points.at(i), Points.at(i + 1), T.getP2()));
		UpTriangles.push_back(Triangle(Points.at(i), Points.at(i + 1), T.getP1()));
	}
	DownTriangles.push_back(Triangle(T.getP3(), T.getP2(), PBetween13));
	return make_pair(DownTriangles, UpTriangles);
}
pair<vector<Triangle>, vector<Triangle>> DivisionTriangleWithCurveBarycenter(Triangle Triang, int NumPointsOnCurve, double time) {
	vector<Triangle> UpTriangles;
	vector<Triangle> DownTriangles;
	vector<Point> Points;
	vector<Point> Points_Add;
	Triangle T = RightTriang(Triang, time);
	Points_Add = PointsOnCurveIntoTriang(T, NumPointsOnCurve, time);

	vector<Point> IntersectInterface = IntersectInterfaceTriangle(T, time);
	Points.push_back(IntersectInterface.at(0));
	for (int i = 0; i < Points_Add.size(); i++)
		Points.push_back(Points_Add.at(i));
	Points.push_back(IntersectInterface.at(1));

	double x_db = 0;
	double y_db = 0;
	double x_ub = 0;
	double y_ub = 0;
	double x = 0;
	double y = 0;
	for (int i = 0; i < Points.size(); i++) {
		x += Points.at(i).getx();
		y += Points.at(i).gety();
	}
	x /= Points.size();
	y /= Points.size();
	if (T.isVertice(Points.at(0)))
		if (T.isVertice(Points.at(Points.size() - 1))) {
			x_db = x;
			y_db = y;
		}
		else {
			x_db = (x + T.getP2().getx()) / 2;
			y_db = (y + T.getP2().gety()) / 2;
		}
	else {
		x_db = (x + T.getP3().getx() + T.getP2().getx()) / 3;
		y_db = (y + T.getP3().gety() + T.getP2().gety()) / 3;
	}
	Point DownBarycenter(x_db, y_db);
	Point UpBarycenter((x + T.getP1().getx()) / 2, (y + T.getP1().gety()) / 2);

	Point PBetween12 = getPointInThisEdge(T.getD1(), Points.at(0), Points.at(Points.size() - 1));
	Point PBetween13 = getPointInThisEdge(T.getD3(), Points.at(0), Points.at(Points.size() - 1));

	for (int i = 0; i < Points.size() - 1; i++) {
		DownTriangles.push_back(Triangle(Points.at(i), Points.at(i + 1), DownBarycenter));
		UpTriangles.push_back(Triangle(Points.at(i), Points.at(i + 1), UpBarycenter));
	}
	DownTriangles.push_back(Triangle(T.getP2(), DownBarycenter, PBetween12));
	DownTriangles.push_back(Triangle(T.getP3(), DownBarycenter, PBetween13));
	DownTriangles.push_back(Triangle(T.getP2(), DownBarycenter, T.getP3()));
	UpTriangles.push_back(Triangle(T.getP1(), UpBarycenter, Points.at(0)));
	UpTriangles.push_back(Triangle(T.getP1(), UpBarycenter, Points.at(Points.size() - 1)));
	return make_pair(DownTriangles, UpTriangles);
}
pair<double, double> QuadrMethodWithInterface(Triangle T, int k_q, int l_q, int NumPointsOnCurve, int cases, double time) {
	vector<Triangle> UpTriangles;
	vector<Triangle> DownTriangles;
	pair<vector<Triangle>, vector<Triangle>> DownAndUpTriangles;
	if (cases == 1)
		DownAndUpTriangles = DivisionTriangleWithCurveOneVertice(T, NumPointsOnCurve, time);
	else if (cases == 2)
		DownAndUpTriangles = DivisionTriangleWithCurveProjection(T, NumPointsOnCurve, time);
	else if (cases == 3)
		DownAndUpTriangles = DivisionTriangleWithCurveBarycenter(T, NumPointsOnCurve, time);
	else {
		cout << "wrong number of cases";
		return make_pair(0, 0);
	}

	DownTriangles = DownAndUpTriangles.first;
	UpTriangles = DownAndUpTriangles.second;
	double integ_down = 0;
	for (int i = 0; i < DownTriangles.size(); i++)
		integ_down += QuadrMethod(DownTriangles.at(i), k_q, l_q);
	double integ_up = 0;
	for (int i = 0; i < UpTriangles.size(); i++)
		integ_up += QuadrMethod(UpTriangles.at(i), k_q, l_q);
	return make_pair(integ_down, integ_up);
}
void DrawTriangle(Triangle T, int NumPointsOnCurve, int cases, double time) {
	vector<Edge> AllEdges;
	vector<Triangle> AllTriangles;
	pair<vector<Triangle>, vector<Triangle>> DownAndUpTriangles;

	ofstream fout;
	if (cases == 1) {
		DownAndUpTriangles = DivisionTriangleWithCurveOneVertice(T, NumPointsOnCurve, time);
		fout.open("DrawTriangleOneVertice.m");
	}
	else if (cases == 2) {
		DownAndUpTriangles = DivisionTriangleWithCurveProjection(T, NumPointsOnCurve, time);
		fout.open("DrawTriangleProjection.m");
	}
	else if (cases == 3) {
		DownAndUpTriangles = DivisionTriangleWithCurveBarycenter(T, NumPointsOnCurve, time);
		fout.open("DrawTriangleBarycenter.m");
	}
	else {
		cout << "wrong number of cases" << endl;
	}

	for (int i = 0; i < DownAndUpTriangles.first.size(); i++) {
		DownAndUpTriangles.first.at(i).ChangeColorTriang('g');
		AllTriangles.push_back(DownAndUpTriangles.first.at(i));
	}
	for (int i = 0; i < DownAndUpTriangles.second.size(); i++)
		AllTriangles.push_back(DownAndUpTriangles.second.at(i));

	for (int iter = 0; iter < AllTriangles.size(); iter++) {
		AllEdges.push_back(AllTriangles.at(iter).getD1());
		AllEdges.push_back(AllTriangles.at(iter).getD2());
		AllEdges.push_back(AllTriangles.at(iter).getD3());
	}
	fout << "clf()" << endl;
	for (int j = 0; j < AllEdges.size(); j++) {
		fout << "drawLine([" << AllEdges.at(j).getA().getx() << "," << AllEdges.at(j).getA().gety() << "], " << "[" << AllEdges.at(j).getB().getx() << ","
			<< AllEdges.at(j).getB().gety() << "], '" << AllEdges.at(j).getcolour() << "');" << endl;
		fout << "hold on" << endl;
	}
	fout << "x = 0:0.01:1;" << endl;
	fout << "y = (x+0.5).*(0.5-x) + " << time << ";" << endl;
	fout << "plot(x, y, 'r');" << endl;
	fout << "axis([0 1 0 1.2]);" << endl;
	fout.close();
}

