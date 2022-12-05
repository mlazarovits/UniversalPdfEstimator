#include "BezierCurve.hh"

#include <iostream>
using std::cout;
using std::endl;


BezierCurve::BezierCurve(double *cp, int n_cps, double* t, int nSamples){
	_cp = cp;
	//order of Bezier Curve: n = n_cps - 1
	_n = n_cps - 1;
	_nSamples = nSamples;
	for(int i = 0; i < _nSamples+1; i++){
		t[i] = double(i)/double(_nSamples);
		if(t[i] > 1 || t[i] < 0){
			cout << "Invalid t input: " << t[i] << endl;
			break;
		}
	}
	_t = t;
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


void BezierCurve::CalculateCurve(double* x){
	//sample r from 0 to 1 nSamples times
	for(int i = 0; i < _nSamples+1; i++){
		x[i] = 0.;
		//sum over n control points
		for(int k = 0; k < _n+1; k++){
			x[i] += _cp[k]*BernsteinPolynomial(_t[i],_n,k);
			//if(i == _nSamples || i == _nSamples - 1) cout << "	n:" << _n << " k:" << k << " x:" <<  x[i] << " cp:" << _cp[k] << " poly:" << BernsteinPolynomial(_t[i],_n,k) << endl;
		}
		//cout << "i: " << i << " r: " << _t[i] << " x: " << x[i] << endl;
	}

}

void BezierCurve::GetInputArray(double *t){
	t = _t;
}

void BezierCurve::CalculateWeightedCurve(double* x, double* w){
	//sample r from 0 to 1 nSamples times
//	cout << "order: " << _n << endl;
	double norm;
	for(int i = 0; i < _nSamples+1; i++){
		x[i] = 0.;
		norm = 0;
		//sum over n control points
		for(int k = 0; k < _n+1; k++){
			x[i] += _cp[k]*w[k]*BernsteinPolynomial(_t[i],_n,k);
			norm += w[k]*BernsteinPolynomial(_t[i],_n,k);
			//if(i == _nSamples || i == _nSamples - 1) cout << "	n:" << _n << " k:" << k << " x:" <<  x[i] << " cp:" << _cp[k] << " poly:" << BernsteinPolynomial(_t[i],_n,k) << endl;
		}
		x[i] /= norm;
		//cout << "i: " << i << " r: " << _t[i] << " x: " << x[i] << endl;
	}
}




void BezierCurve::CalculateWeightedCurve_MultiOrder(double* x, double* w, double o){
	//sample r from 0 to 1 nSamples times
	//cout << "order: " << _n*o << " with " << _n+1 << " terms" <<  endl;
	double norm;
	for(int i = 0; i < _nSamples+1; i++){
		x[i] = 0.;
		norm = 0;
		//sum over n control points
		for(int k = 0; k < _n+1; k++){
			x[i] += _cp[k]*w[k]*BernsteinPolynomial(_t[i],_n*o,k*o);
			norm += w[k]*BernsteinPolynomial(_t[i],_n*o,k*o);
			//if(i == _nSamples || i == _nSamples - 1) cout << "	n:" << _n << " k:" << k << " x:" <<  x[i] << " cp:" << _cp[k] << " poly:" << BernsteinPolynomial(_t[i],_n,k) << endl;
		}
		x[i] /= norm;
		//cout << "i: " << i << " r: " << _t[i] << " x: " << x[i] << endl;
	}
}






