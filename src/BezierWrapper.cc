#include "BezierWrapper.hh"
#include "RandomSample.hh"
#include "BezierCurve.hh"
#include "NumDeriv.hh"
#include <iostream>


#include <TLatex.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TFile.h>

using std::cout;
using std::endl;

BezierWrapper::BezierWrapper(const int nCP){ 
	//set intermediate (x,r) points, but (0,0) and (1,1) are set
	_nDataPts = nCP;
	_nDataPts += 2;
	
	Sample();

}


BezierWrapper::~BezierWrapper(){ }


void BezierWrapper::Sample(unsigned long long seed, double xmin, double xmax, double mu, double sigma){
	//initialize member array pointers
	double tmp_global[_nSamples];
	_x_global = tmp_global;
	_r_global = tmp_global;
	double tmp[_nDataPts];
	_x = tmp;
	_r = tmp;
	//simulate data for control points	
	RandomSample rs(seed);
	rs.SetRange(xmin,xmax);
	rs.SampleGaussian(mu, sigma, _nSamples, _x_global);
	rs.SampleGaussian(mu,sigma, _nSamples, _r_global);
	
	//select data points
	rs.SelectPoints(_x_global, _nSamples, _x, _nDataPts);
	rs.SelectPoints(_r_global, _nSamples, _r, _nDataPts);
	_x[0] = 0.;
	_r[0] = 0.;
	_x[_nDataPts-1] = 1.;
	_r[_nDataPts-1] = 1.;
	for(int i = 0; i < _nDataPts; i++) cout << i << " x: " << _x[i] << " r: " << _r[i] << endl;
}



void BezierWrapper::CalculateCurve(){
	//x and r are functions of t
	double tmp[_nSamples+1];
	_x_approx = tmp; //filled by reference in CalculateCurve
	_r_approx = tmp; //filled by reference in CalculateCurve
	_t = tmp; //filled by reference in CalculateCurve

	BezierCurve x_bc(_x, _nDataPts, _t, _nSamples);
	BezierCurve r_bc(_r, _nDataPts, _t, _nSamples);
	x_bc.CalculateCurve(_x_approx);
	r_bc.CalculateCurve(_r_approx);
	cout << "calculated Bezier curves for x and r - order " << _nDataPts-2 << endl;


}


void BezierWrapper::Transform(int n){
	_nTransforms = n;
	//x_old = x
	
	for(int i = 0; i < _nTransforms; i++){
	       	//x_new = transform(x_old);
		//x_old = x_new
	}
}

void BezierWrapper::SetNTransforms(int n){
	_nTransforms = n;
}

void BezierWrapper::_transform(){
	double x_transf[_nDataPts];
	x_transf[0] = 0.;
	x_transf[_nDataPts-1] = 1.;
	for(int c = 1; c < _nDataPts-1; c++){
		double diff = 999.;
		int p = 0; //idx of t we want
		for(int i = 0; i < _nSamples; i++){
			if( fabs(_r_approx[i] - _r[c]) < diff){
			       	diff = fabs(_r_approx[i] - _r[c]);
				p = i;
			}	
		}
		x_transf[c] = _x_approx[p];
	}


}

void BezierWrapper::CalculateJacobian(){
	//GET JACOBIAN
	double jac[_nSamples];
	double dr_dt[_nSamples];
	double dt_dx[_nSamples];
	NumDeriv drdt(_r_approx, _t, _nSamples);
	NumDeriv dtdx(_t, _x_approx, _nSamples);
	
	//make sure domain is [0,1]
	drdt.MapToInterval(0,1);
	dtdx.MapToInterval(0,1);

	drdt.FiniteDiff(dr_dt);
	dtdx.FiniteDiff(dt_dx);

	for(int i = 0; i < _nSamples; i++) jac[i] = dr_dt[i]*dt_dx[i];


}


void BezierWrapper::Visualize(){
	cout << "# data pts: " << _nDataPts-2 << endl;
	TGraph* gr = new TGraph(_nSamples+1, _r_approx, _x_approx);
	
	gr->SetLineWidth(1);
	gr->SetLineColor(kRed-7);
	gr->SetMarkerStyle(15);

	cout << "Plots stored in: " << "test/" << _fname << endl;
	TFile* f = TFile::Open(("test/"+_fname).c_str(),"RECREATE");
	TGraph* gr_data = new TGraph(_nDataPts, _r, _x);
	TLine* line = new TLine(0,0,1,1);
       	


	gr_data->SetMarkerSize(0.95);
       	gr_data->SetLineWidth(0);
	gr_data->SetMarkerStyle(20);
       	gr_data->SetMarkerColor(kViolet+7);
	
	
	_mg->Add(gr_data);
	_mg->Add(gr);
	
	_leg->AddEntry(gr_data, "Original data/control points");
	_leg->AddEntry(gr,"3rd order Bezier approximation");
	
	
	
	
	_cv->cd();
	_cv->SetRightMargin(0.09);
	_cv->SetLeftMargin(0.15);
	_cv->SetBottomMargin(0.15);
	

	_mg->SetTitle("");
	_mg->Draw("ALP");
	_mg->GetXaxis()->SetTitle("r - rank space");
	_mg->GetYaxis()->SetTitle("x - coord. space");
	_leg->Draw("same");

	line->SetLineStyle(7);
	line->Draw("ALsame");

	TLatex l;
	l.SetNDC();
	l.SetTextAlign(11);
	l.SetTextSize(0.03);
	l.DrawLatex(0.18,0.92,"Bezier Curves to Independently Approximate x and r");
	
	f->cd();
	_cv->Write();
	f->Close();
	

}


void BezierWrapper::VisualizeTransform(double* x_transf){
	TGraph* gr_transf = new TGraph(_nDataPts, _r, x_transf);
	gr_transf->SetMarkerSize(0.95);
       	gr_transf->SetLineWidth(0);
	gr_transf->SetMarkerStyle(20);
       	gr_transf->SetMarkerColor(kGreen+7);
	_mg->Add(gr_transf);
	_leg->AddEntry(gr_transf, "Transformed data/control points");


}
