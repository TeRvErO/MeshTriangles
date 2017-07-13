#include "Mesh triangles.h"

Point::Point (): x(0), y(0) {}
Point::Point (double x1, double y1): x(x1), y(y1) {}
void Point::operator= (const Point P){
	x = P.x;
	y = P.y;
}
Point Point::operator- (Point &p){
	return Point (x - p.x, y - p.y);
}
Point operator* (double s, Point &p){
	return Point (s*p.x, s*p.y);
}
bool Point::operator== (Point &p){
	return (x == p.x) && (y == p.y);
}
bool Point::operator!= (Point &p){
	return !(*this == p);
}
bool Point::operator< (Point &p){
	return ((x < p.x) || ((x == p.x) && (y<p.y)));
}
bool Point::operator> (Point &p){
	return ((x>p.x) || ((x == p.x) && (y > p.y)));
}
double Point::length() {
	return sqrt(pow(x, 2) + pow(y, 2));
}
int Point::classify(Point &p0, Point &p1){
	Point p2 = *this;
	Point a = p1 - p0;
	Point b = p2 - p0;
	double sa = a. x * b.y - b.x * a.y;
	if (sa > 0.0)
		return LEFT;
	if (sa < 0.0)
		return RIGHT;
	if ((a.x * b.x < 0.0) || (a.y * b.y < 0.0))
		return BEHIND;
	if (a.length() < b.length())
		return BEYOND;
	if (p0 == p2)
		return ORIGIN;
	if (p1 == p2)
		return DESTINATION;
	return BETWEEN;
}
double Point::dist(const Point P){
	return sqrt(pow(getx() - P.getx(),2) + pow(gety() - P.gety(),2));
}

Edge::Edge (): colour('b') {}
Edge::Edge(Point P, Point Q): A(P), B(Q) {}
void Edge::operator= (const Edge example){
	A = example.A;
	B = example.B;
}
bool Edge::operator== (Edge &E){
	return (this->getA() == E.getA()) && (this->getB() == E.getB());
}
bool Edge::In(Point P) {
	double value = (B.gety() - A.gety()) / (B.getx() - A.getx())*(P.getx() - this->A.getx()) + this->A.gety();
	return (abs(value - P.gety()) < 0.0000001);
}
double Edge::length(){
	return sqrt(pow(A.getx() - B.getx(),2) + pow(A.gety() - B.gety(),2)); 
}
void Edge::ChangeColorEdge(char newcolour){
	colour = newcolour;
}
Point Edge::Projection(Point P){
	double a, b, c, inter_x, inter_y;
	Point intersection;
	if(A.getx() == B.getx())
		intersection = Point(A.getx(), P.gety());
	else if (A.gety() == B.gety())
		intersection = Point(P.getx(), A.gety());
	else {
		a = -1./(B.getx() - A.getx());
		b = 1./(B.gety() - A.gety());
		c = A.getx()/(B.getx() - A.getx()) - A.gety()/(B.gety() - A.gety());
		inter_x = 1/(2*a)*(a*P.getx() - b*P.gety() - c);
		inter_y = 1/(2*b)*(b*P.gety() - a*P.getx() - c);
		intersection = Point(inter_x, inter_y);
	}
	return intersection;
}

Triangle::Triangle (): color('b') {}
Triangle::Triangle (Edge one, Edge two, Edge three){
	Point X = one.getA();
	Point Y = one.getB();
	Point Z;
	if (two.getA() == X || two.getA() == Y)
		Z = two.getB();
	else Z = two.getA();
	P1 = X;
	Point XY = Y - X;
	Point XZ = Z - X;
	if (XY.getx() * XZ.gety() - XY.gety() * XZ.getx() > 0){
		P2 = Y;
		P3 = Z;
	}
	else {
		P2 = Z;
		P3 = Y;
	}
	D1 = Edge(P1, P2);
	D2 = Edge(P2, P3);
	D3 = Edge(P3, P1);
	color = 'b';
}
Triangle::Triangle(Point one, Point two, Point three) {
	*this = Triangle(Edge(one, two), Edge(two, three), Edge(three, one));
}
void Triangle::operator= (const Triangle T){
	this->D1 = T.D1;
	this->D2 = T.D2;
	this->D3 = T.D3;
	this->P1 = T.P1;
	this->P2 = T.P2;
	this->P3 = T.P3;
	this->color = T.color;
}
bool Triangle::pointInTriangle(Point p) {
	return ((p.classify(D1.getA(), D1.getB()) != RIGHT) && (p.classify(D2.getA(), D2.getB()) != RIGHT) && (p.classify(D3.getA(), D3.getB()) != RIGHT));
}
bool Triangle::isVertice(Point a) {
	return (a == P1) || (a == P2) || (a == P3);
}


void Triangle::ChangeColorTriang(char color) {
	D1.ChangeColorEdge(color);
	D2.ChangeColorEdge(color);
	D3.ChangeColorEdge(color);
	this->color = color;
}

Point Triangle::barycenter(){
	return Point(1./3*(this->P1.getx() + this->P2.getx() + this->P3.getx()),1./3* (this->P1.gety() + this->P2.gety() + this->P3.gety()));
}
double Triangle::square(){
	return 1./2.*((P1.getx() - P3.getx())*(P2.gety() - P3.gety()) - (P2.getx() - P3.getx())*(P1.gety() - P3.gety()));
}


Rectangular::Rectangular () {}
Rectangular::Rectangular (Point * P){
	for (int i = 0; i < 4; i++)
		vertices[i] = P[i];
}
Rectangular::Rectangular (Triangle T){
	vertices[0] = Point(min(min(T.getD1().getA().getx(), T.getD1().getB().getx()), T.getD2().getB().getx()), min(min(T.getD1().getA().gety(), T.getD1().getB().gety()), T.getD2().getB().gety()));
	vertices[1] = Point(min(min(T.getD1().getA().getx(), T.getD1().getB().getx()), T.getD2().getB().getx()), max(max(T.getD1().getA().gety(), T.getD1().getB().gety()), T.getD2().getB().gety()));
	vertices[2] = Point(max(max(T.getD1().getA().getx(), T.getD1().getB().getx()), T.getD2().getB().getx()), max(max(T.getD1().getA().gety(), T.getD1().getB().gety()), T.getD2().getB().gety()));
	vertices[3] = Point(max(max(T.getD1().getA().getx(), T.getD1().getB().getx()), T.getD2().getB().getx()), min(min(T.getD1().getA().gety(), T.getD1().getB().gety()), T.getD2().getB().gety()));
}	


ostream& operator<<(std::ostream& out, const Point& P){
	return out <<"("<<P.getx()<<","<<P.gety()<<")";
}
ostream& operator<<(std::ostream& out, Edge& E){
	return out <<E.getA()<<"--"<<E.getB();
}
ostream& operator<<(std::ostream& out, Triangle& T){
	return out <<T.getD1()<<" + "<<T.getD2()<<" + "<<T.getD3();
}
ostream& operator<<(std::ostream& out, Rectangular& R){
	Point* vertices = R.getPoints();
	return out <<vertices[0]<<" -- "<<vertices[1]<<" -- "<<vertices[2]<<" -- "<<vertices[3];
}

