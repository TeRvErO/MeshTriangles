#include <iostream>
#include <conio.h>
#include <string>
#include<math.h>
#include <algorithm>
#include <fstream>
#include<vector>
#include<cmath>
#include<limits>
using namespace std;
enum {LEFT,  RIGHT,  BEYOND,  BEHIND, BETWEEN, ORIGIN, DESTINATION};
class Point{
	double x;
	double y;
public:
	Point ();
	Point (double x1, double y1);
	void operator= (const Point P);
	Point operator- (Point &p);
	friend Point operator* (double s, Point &p);
	bool operator== (Point &p);
	bool operator!= (Point &p);
	bool operator< (Point &p);
	bool operator> (Point &p);
	double length();
	int classify(Point &p0, Point &p1);
	double getx() const { return x; }
	double gety() const { return y; }
	double dist(const Point P);
};

class Edge{
	Point A;
	Point B;
	char colour;
public:
	Edge ();
	Edge(Point P, Point Q);
	void operator= (const Edge example);
	bool operator== (Edge &E);
	bool In(Point P);
	Point getA(){ return A;}
	Point getB(){return B;}
	char getcolour() {return colour;}
	double length();
	void ChangeColorEdge(char newcolour);
	Point Projection(Point P);
};

class Triangle{
	Edge D1;
	Edge D2;
	Edge D3;
	Point P1;
	Point P2;
	Point P3;
	char color;
public:
	Triangle ();
	Triangle (Edge one, Edge two, Edge three);
	Triangle(Point one, Point two, Point three);
	void operator= (const Triangle T);
	bool pointInTriangle(Point p);
	bool isVertice(Point a);
	Edge getD1(){ return D1;}
	Edge getD2(){ return D2;}
	Edge getD3(){ return D3;}
	Point getP1(){ return P1;}
	Point getP2(){ return P2;}
	Point getP3(){ return P3;}
	char getcolour() { return color; }
	void ChangeColorTriang(char color);
	
	Point barycenter();
	double square();
	
};

class Rectangular{
	Point vertices[4];
public:
	Rectangular ();
	Rectangular (Point * P);
	Rectangular (Triangle T);	
	Point* getPoints(){ return vertices;}
};

ostream& operator<<(std::ostream& out, const Point& P);
ostream& operator<<(std::ostream& out, Edge& E);
ostream& operator<<(std::ostream& out, Triangle& T);
ostream& operator<<(std::ostream& out, Rectangular& R);


/*class Eleve{
	string nom;
	int note;
public:
	Eleve(string eleve, int note_eleve){
		nom = eleve;
		note = note_eleve;
	}
	string operator()(int i) const{
		if (i==0) return nom;
	}
	int operator[](int i) const{
		if (i==1) return note;
	}
	string getnom(){
		return nom;
	}
	int getnote(){
		return note;
	}
};

class CompareNom{
	int* cmp;
public:
	CompareNom( int* i){
		cmp = i;
	}
	bool operator ()(Eleve first, Eleve second) const{
		++*cmp;
		return first.getnom() < second.getnom();
	}
	int getiter(){
		return *cmp;
	}
};
class CompareNote{
	int* compt;
public:
	CompareNote( int* i){
		compt = i;
	}
	bool operator ()(Eleve first, Eleve second) const{
		++*compt;
		return first.getnote() > second.getnote();
	}
	int getiter(){
		return *compt;
	}
};
ostream& operator<<(ostream& str, const Eleve& A){
		str<< A(0) << "---" << A[1]<<endl;
		return str;
	}

template <class T>
void print(const vector<T>& V){
	vector<T>::const_iterator it = V.begin();
	for(; it != V.end(); ++it)
		cout<< *it;
	cout<<endl;
}*/