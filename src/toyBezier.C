#include <iostream>
#include "RandomSample.hh"
#include "BezierCurve.hh"

#include <TH1D.h>
#include <TFile.h>

using std::cout;
using std::endl;
using std::string;
int main(int argc, char *argv[]){
	
	cout << "Hi welcome to Chili's!" << endl;
	
	unsigned long long seed = 1234;
	int nSample = 3;
	double x[nSample];
	double mu = 0.;
	double sigma = 1.;
	double xmax = 5;
	double xmin = -xmax;
	bool viz = false;
	string fname;
	if(argc > 1){
		viz = true;
		fname = string(argv[1]);
	}
	
	RandomSample rs(seed);
	double cps[] = {-1.5, 0., 1.0};
	BezierCurve bc(cps, 3);

	rs.SetRange(xmin,xmax);

	rs.SampleGaussian(mu, sigma, nSample, x);
	
	if(viz){
		TH1D* hist = new TH1D("hist","hist",20,xmin,xmax);	
		TFile* f = TFile::Open(("test/"+fname).c_str(),"RECREATE");
		for(auto i : x) hist->Fill(i);
		f->cd();
		hist->Write();
		f->Close();
	}
	else
		for(auto i : x) cout << i << endl;







}
