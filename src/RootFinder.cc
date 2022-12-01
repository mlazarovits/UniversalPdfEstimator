#include "RootFinder.hh"
#include <iostream>

using std::cout;
using std::endl;

//TODO: needs to be debugged - make sure root finder works

RootFinder::RootFinder(double* x, double* y, int n){
	_x = x;
	_y = y;
	_n = n;
	_deriv = NumDeriv(_x, _y, _n);
}

RootFinder::~RootFinder(){ }


void RootFinder::SetTolerance(double tol){
	_tol = tol;
}


double RootFinder::NewtonsMethod(int nIt, double x0){
	double xNew;
	double xOld = x0;
	double arr[nIt];

	//numerically calculate derivative
	double fPrime[_n];
       	_deriv.FiniteDiff(fPrime);
	arr[0] = x0;

	for(int i = 1; i < nIt; i++){
		xNew = xOld - ( _y[i-1] / fPrime[i-1] );
		arr[i] = xNew;
		if( fabs(xNew - xOld) < _tol){
			cout << "Newton's Method has converged with tolerance: " << _tol << endl;
			return xNew;
		}
		xOld = xNew;
	}	
	cout << "Newton's method failed after " << nIt << " iterations." << endl;
	return 0;
}
