#include <iostream>
#include "RandomSample.hh"
#include "BezierCurve.hh"
#include "NumDeriv.hh"

#include <TH1D.h>
#include <TFile.h>

using std::cout;
using std::endl;
using std::string;
int main(int argc, char *argv[]){
	string fname;
	bool viz = false;
	if(argc > 1){
		viz = true;
		fname = string(argv[1]);
	}
	
	cout << "Hi welcome to Chili's!" << endl;
	
	unsigned long long seed = 1234;
	int nSample = 3;
	double x[nSample];
	double mu = 0.;
	double sigma = 1.;
	double xmax = 5;
	double xmin = -xmax;
	
cout << "A: " << nSample << endl;
	RandomSample rs(seed);
	rs.SetRange(xmin,xmax);
	rs.SampleGaussian(mu, sigma, nSample, x);
	
cout << "B: " << nSample << endl;
	int nCPs = 3;
	double cps[3] = {-1.5, 0., 1.0};
	double r[nSample];
	BezierCurve bc(cps, nCPs);
	bc.CalculateCurve(x, nSample, r);
cout << "C: " << nSample << endl;
	double deriv[nSample-1];
	NumDeriv derivate(x,r,nSample);
	derivate.MapToInterval(0.,1.);
	derivate.FiniteDiff(deriv);
	
	for(int i = 0; i < nSample; i++) cout << "point #: " << i << " x: " << x[i] << " r: " << r[i] << " derivative: " << deriv[i] << endl;
	if(viz){
		TH1D* hist = new TH1D("hist","hist",20,xmin,xmax);	
		TFile* f = TFile::Open(("test/"+fname).c_str(),"RECREATE");
		for(auto i : x) hist->Fill(i);
		f->cd();
		hist->Write();
		f->Close();
	}







}
