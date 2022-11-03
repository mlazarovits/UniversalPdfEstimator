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
	if(k > n){
		cout << "k must be less than n." << endl;
		return -999;
	}
	if (k == 0 || k == n)
		return 1;
	return BinomialCoeff(n - 1, k - 1) + BinomialCoeff(n - 1, k);
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


double BezierCurve::BernsteinPolynomial(double t, int n, int k){
	//k should run from 0 to n
	if(k > n || k < 0){ cout << "Error: k must be 0 < k < n. k = " << k << endl; return 0; }
	
	//check t is between [0,1]
	if(t < 0 || t > 1){ cout << "Error: t must be 0 < t < 1. t = " << t << endl; return 0.; }
	
	double coeff = BinomialCoeff(n,k);

	return coeff*pow((1 - t),n-k)*pow(t,k);
}


void BezierCurve::CalculateCurve(double* x, int nSamples, double* r){
	//map x to t from [0,1]
	double xmax = FindExtremum(x,nSamples,true);
	double xmin = FindExtremum(x,nSamples,false);
	double xint = (xmax - xmin); 
	double t[nSamples];
	for(int i = 0; i < nSamples; i++){
		t[i] = (x[i] - xmin)/xint;
		r[i] = 0.;
		//sum over n control points
		for(int k = 0; k < _n+1; k++){
			r[i] += _cp[k]*BernsteinPolynomial(t[i],_n,k);
		}
	}
//	return r;
}

