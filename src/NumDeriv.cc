#include "NumDeriv.hh"
#include <iostream>

using std::cout;
using std::endl;

NumDeriv::NumDeriv(double* x, double* y, int n){
	_x = x;
	_y = y;
	_n = n;
//	for(int i = 0; i < _n; i++) cout << i << " x: " << _x[i] << " y: " << _y[i] << " x+1: " << _x[i+1] << " y+1: " << _y[i+1] << endl;
}

NumDeriv::~NumDeriv(){ }



void NumDeriv::FiniteDiff(double* deriv){
	for(int i = 0; i < _n; i++){
		deriv[i] = (_y[i+1]-_y[i])/(_x[i+1]-_x[i]);
	}
}

void NumDeriv::MapToInterval(double new_min, double new_max){
	double xmax = -999;
	double xmin = 999;
	for(int i = 0; i < _n; i++){
		if(_x[i] > xmax) xmax = _x[i];
		if(_x[i] < xmin) xmin = _x[i];
	}
//check that given input is within desired interval
	if(xmax <= new_max && xmin >= new_min) return;
	
	for(int i = 0; i < _n; i++){
	       	_x[i] = new_min + ((new_max - new_min)/(xmax - xmin))*(_x[i] - xmin);
	}
}
