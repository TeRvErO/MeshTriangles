// Examples of interfaces
double Circle(Point A);
double CircleYofXpositive(Point A);
double Parabola(Point A);
double ParabolaYofX(double x);

// The equation of edge y = a*x + b
double EquationOfEdge(double x, double coeff[2]);
// Integrating function for curvilinear integral (minus primitive functon by y of initial one) 
double IntegFunc(double x, double y);
// Corresponding point to z for standart edge [-1,1] 
double CorrespPointToStEdge(Edge E, double z);
// Gaussian quadrature with 5 points for calculating integral on edge
double GaussianQuadrForEdge(Edge E, double * PointsOnStEdge, double * weights, int NumPoints);
// Gaussian quadrature  with 5 points for calculating integral on part of interface which intersect triangle in points A and B
double GaussianQuadrForCurve(Point A, Point B, double * PointsOnStEdge, double * weights, int NumPoints);
// Value of integral for curve triangle (need of improvement - for the arbitrary triangle)
double EffIntegration(Triangle T);
// Draw triangle and interface
void PrintSituation(Triangle T);