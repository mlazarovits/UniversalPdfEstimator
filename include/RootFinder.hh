#ifndef RootFinder_HH
#define RootFinder_HH

#include <cmath>
#include "NumDeriv.hh"


class RootFinder{
	public:
		RootFinder(double* x, double* y, int n);
		virtual ~RootFinder();
		double NewtonsMethod(int nIt, double x0);
		void SetTolerance(double tol);

		
	private:
		double* _x;
		double* _y;
		int _n;
		double _tol = pow(10,-10);
		NumDeriv _deriv;


};







#endif
