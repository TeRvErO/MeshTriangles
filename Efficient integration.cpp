#include "Mesh triangles.h"
#include "Efficient integration.h"
// Examples of interfaces
double Circle(Point A) {
	return pow(A.gety(), 2) + pow(A.getx(), 2) - 2;
}
double CircleYofXpositive(Point A) {
	return sqrt(2 - pow(A.getx(), 2));
}
double Parabola(Point A) {
	return A.gety() - pow(A.getx(), 2);
}
double ParabolaYofX(double x) {
	return pow(x, 2);
}

// The equation of edge y = a*x + b
double EquationOfEdge(double x, double coeff[2]) {
	return x*coeff[0] + coeff[1];
}
// Integrating function for curvilinear integral (minus primitive functon by y of initial one) 
double IntegFunc(double x, double y) {
	return -y;
}
// Corresponding point to z for standart edge [-1,1]
double CorrespPointToStEdge(Edge E, double z) {
	return z*(E.getB().getx() - E.getA().getx()) / 2. + (E.getB().getx() + E.getA().getx()) / 2.;
}
// Gaussian quadrature with 5 points for calculating integral on edge
double GaussianQuadrForEdge(Edge E, double * PointsOnStEdge, double * weights, int NumPoints) {
	double coeff[2];
	coeff[0] = (E.getB().gety() - E.getA().gety()) / (E.getB().getx() - E.getA().getx());
	coeff[1] = E.getB().gety() - (E.getB().getx()*(E.getB().gety() - E.getA().gety()) / ((E.getB().getx() - E.getA().getx())));
	
	double a = 0;
	for (int i = 0; i < NumPoints; i++) {
		a += weights[i] * IntegFunc(CorrespPointToStEdge(E, PointsOnStEdge[i]), EquationOfEdge(CorrespPointToStEdge(E, PointsOnStEdge[i]), coeff));
	}
	//double a = (5. / 9.*IntegFunc(CorrespPointToStEdge(E, -sqrt(0.6)), EquationOfEdge(CorrespPointToStEdge(E, -sqrt(0.6)), coeff)) + 8. / 9.*IntegFunc(CorrespPointToStEdge(E, 0), EquationOfEdge(CorrespPointToStEdge(E, 0), coeff)) + 5. / 9.*IntegFunc(CorrespPointToStEdge(E, sqrt(0.6)), EquationOfEdge(CorrespPointToStEdge(E, sqrt(0.6)), coeff)));
	return (E.getB().getx() - E.getA().getx()) / 2.* a;
}
// Gaussian quadrature  with 5 points for calculating integral on part of interface which intersect triangle in points A and B
double GaussianQuadrForCurve(Point A, Point B, double * PointsOnStEdge, double * weights, int NumPoints) {
	Edge E(A, B);
	double a = 0;
	for (int i = 0; i < NumPoints; i++) {
		a += weights[i] * IntegFunc(CorrespPointToStEdge(E, PointsOnStEdge[i]), ParabolaYofX(CorrespPointToStEdge(E, PointsOnStEdge[i])));
	}
	//double a = (5. / 9.*IntegFunc(CorrespPointToStEdge(E, -sqrt(0.6)), ParabolaYofX(CorrespPointToStEdge(E, -sqrt(0.6)))) + 8. / 9.*IntegFunc(CorrespPointToStEdge(E, 0), ParabolaYofX(CorrespPointToStEdge(E, 0))) + 5. / 9.*IntegFunc(CorrespPointToStEdge(E, sqrt(0.6)), ParabolaYofX(CorrespPointToStEdge(E, sqrt(0.6)))));
	return (B.getx() - A.getx()) / 2.*a;
}
// Value of integral for curve triangle (need of improvement - for the arbitrary triangle)
double EffIntegration(Triangle T) {
	//intersection points
	Point one(-1, 1);
	Point two(1, 1);
	Edge E1(T.getP1(), one);
	Edge E2(two, T.getP1());
	double I1, I2, I3;
	const int NumPoints = 5;
	double PointsOnStEdge[NumPoints];
	PointsOnStEdge[0] = 0;
	PointsOnStEdge[1] = sqrt(5. - 2.*sqrt(10. / 7.)) / 3.;
	PointsOnStEdge[2] = -sqrt(5. - 2.*sqrt(10. / 7.)) / 3.;
	PointsOnStEdge[3] = sqrt(5. + 2.*sqrt(10. / 7.)) / 3.;
	PointsOnStEdge[4] = -sqrt(5. + 2.*sqrt(10. / 7.)) / 3.;
	double weights[NumPoints];
	weights[0] = 128. / 225.;
	weights[1] = (322. + 13.*sqrt(70)) / 900.;
	weights[2] = (322. + 13.*sqrt(70)) / 900.;
	weights[3] = (322. - 13.*sqrt(70)) / 900.;
	weights[4] = (322. - 13.*sqrt(70)) / 900.;
	I1 = GaussianQuadrForEdge(E1, PointsOnStEdge, weights, NumPoints);
	I2 = GaussianQuadrForCurve(E1.getB(), E2.getA(), PointsOnStEdge, weights, NumPoints);
	I3 = GaussianQuadrForEdge(E2, PointsOnStEdge, weights, NumPoints);
	return I1 + I2 + I3;
}
// Draw triangle and interface
void PrintSituation(Triangle T) {
	vector<Edge> AllEdges;
	AllEdges.push_back(T.getD1());
	AllEdges.push_back(T.getD2());
	AllEdges.push_back(T.getD3());
	ofstream fout;
	fout.open("DrawTriangleAndCircle.m");
	fout << "clf()" << endl;
	for (int j = 0; j < AllEdges.size(); j++) {
		fout << "drawLine([" << AllEdges.at(j).getA().getx() << "," << AllEdges.at(j).getA().gety() << "], " << "[" << AllEdges.at(j).getB().getx() << ","
			<< AllEdges.at(j).getB().gety() << "], '" << AllEdges.at(j).getcolour() << "');" << endl;
		fout << "hold on" << endl;
	}
	fout << "x = -1.5:0.01:1.5;" << endl;
	fout << "y = x.*x;" << endl;
	fout << "plot(x, y, 'r');" << endl;
	fout << "axis([-3 3 0 3]);" << endl;
	fout.close();
}