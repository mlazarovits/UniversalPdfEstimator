#include "BezierCurve.hh"

#include <iostream>
using std::cout;
using std::endl;


BezierCurve::BezierCurve(double *cp, int n_cps){
	_cp = cp;
	_n = n_cps;
//	for(int i = 0; i < n_cps; i++) cout << _cp[i] << endl;
}


BezierCurve::~BezierCurve(){ }


double BezierCurve::BinomialCoeff(int n,int k){
	if(k > _n){
		cout << "k must be less than n." << endl;
		return -999;
	}
	if (k == 0 || k == _n)
		return 1;
	return BinomialCoeff(_n - 1, k - 1) + BinomialCoeff(_n - 1, k);
}


double BezierCurve::FindExtremum(double *x, int n, bool max){
	double ret;
	if(max) ret = -999;
	else ret = 999;
	for(int i = 0; i < n; i++){
		if(max){ if(x[i] > ret) ret = x[i];}
		else{ if(x[i] < ret) ret = x[i];}
	}
	return ret;
}


double BezierCurve::BernsteinPolynomial(double* x){
	//map x to t from [0,1]
	int max = FindExtremum(x,_n,true);

	//k runs from 0 to _n
	int k = 1;
	double bc = BinomialCoeff(_n,k);

	return bc;//*pow(1-x,n-k)*pow(x,k);
}

