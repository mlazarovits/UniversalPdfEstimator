#ifndef BezierCurve_HH
#define BezierCurve_HH

#include <cmath>

class BezierCurve{
	public:
		BezierCurve(double *cp, int n_cps, double* r, int nSamples);
		virtual ~BezierCurve();
		double BernsteinPolynomial(double t, int n, int k);
		double FindExtremum(double *x, int n, bool max);
		void CalculateCurve(double* x);
		void GetInputArray(double *r);

		//control points
		double* _cp;
		//number of control points/highest order of Bernstein polynomial
		int _n; //see if i can make this const so nothing in this class can change it once it's assigned
		int _nSamples; //samples along input axis
		double* _r; //input axis
	private:
		double BinomialCoeff(int n, int k);









};


#endif
