#include "NumDeriv.hh"
#include <iostream>

using std::cout;
using std::endl;

NumDeriv::NumDeriv(double* x, double* y, int n){
	_x = x;
	_y = y;
	_n = n;
}

NumDeriv::~NumDeriv(){ }



void NumDeriv::FiniteDiff(double* deriv){
	cout << "FiniteDiff" << endl;
	for(int i = 0; i < _n-1; i++){
		deriv[i] = (_x[i+1]-_x[i])/(_y[i+1]-_y[i]);
		cout << "deriv: " << deriv[i] << " x: " << _x[i] << " y: " << _y[i] << endl;
	}
}

void NumDeriv::MapToInterval(double new_min, double new_max){
	double xmax = -999;
	double xmin = 999;
	for(int i = 0; i < _n; i++){
		if(_x[i] > xmax) xmax = _x[i];
		if(_x[i] < xmin) xmin = _x[i];
	}
	for(int i = 0; i < _n; i++) _x[i] = new_min + ((new_max - new_min)/(xmax - xmin))*(_x[i] - xmin);

}
