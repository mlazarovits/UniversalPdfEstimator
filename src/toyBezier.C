#include <iostream>
#include "RandomSample.hh"
#include "BezierCurve.hh"
#include "NumDeriv.hh"

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
using std::string;
int main(int argc, char *argv[]){
	string fname = "test.root";
	bool viz = false;
	unsigned long long seed = 123;
	int nDataPts = 1;
	double mu = 0.5;
	double sigma = 0.25;
	double xmax = 1;
	double xmin = 0;
	bool hprint = false;
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i],"--help", 6) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"-h", 2) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"--output", 8) == 0){
     			i++;
    	 		fname = string(argv[i]);
   		}
		if(strncmp(argv[i],"-o", 2) == 0){
     			i++;
    	 		fname = string(argv[i]);
   		}
		if(strncmp(argv[i],"--viz",5) == 0){
			viz = true;
		}
		if(strncmp(argv[i],"-n", 2) == 0){
     			i++;
    	 		nDataPts = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--nDataPts", 10) == 0){
     			i++;
    	 		nDataPts = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--mu", 4) == 0){
     			i++;
    	 		mu = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--sigma", 7) == 0){
     			i++;
    	 		sigma = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--xmax", 6) == 0){
     			i++;
    	 		xmax = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--xmin", 6) == 0){
     			i++;
    	 		xmin = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--seed", 6) == 0){
     			i++;
    	 		seed = std::stoi(argv[i]);
   		}
	
	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)          print options" << endl;
   		cout << "   --ouput(-o) [file]  output root file (in test/)" << endl;
   		cout << "   --viz       	turns on visualization" << endl;
   		cout << "   --nDataPts(-n) [n]  sets number of data points/control points (order of Bezier curve)" << endl;
   		cout << "   --mu                sets mean of Gaussian to sample data from" << endl;
   		cout << "   --sigma             sets sigma of Gaussian to sample data from" << endl;
   		cout << "   --xmax              sets x-maximum of Gaussian to sample data from (xmax > xmin)" << endl;
   		cout << "   --xmin              sets x-minimum of Gaussian to sample data from (xmax > xmin)" << endl;
   		cout << "Example: ./toyBezier.x -n 10 --mu 0 --sigma 1 --xmax 3 --xmin -3 --viz -o testViz.root" << endl;

   		return 0;
  }
	
	cout << "Hi welcome to Chili's!" << endl;
	
//set intermediate (x,r) points, but (0,0) and (1,1) are set
	nDataPts += 2;
	int nSamples = 1000;//number of samples for t NOT including endpoints and x_global and r_global
	double x[nDataPts];
	double r[nDataPts];
	double x_global[nSamples];
	double r_global[nSamples];
	TGraph* gr;

//simulate data for control points	
	RandomSample rs(seed);
	rs.SetRange(xmin,xmax);
	rs.SampleGaussian(mu, sigma, nSamples, x_global);
	rs.SampleGaussian(mu,sigma, nSamples, r_global);
	
	
	//test SelectPoints
	rs.SelectPoints(x_global, nSamples, x, nDataPts);
	rs.SelectPoints(r_global, nSamples, r, nDataPts);
	x[0] = 0.;
	r[0] = 0.;
	x[nDataPts-1] = 1.;
	r[nDataPts-1] = 1.;
	for(int i = 0; i < nDataPts; i++) cout << i << " x: " << x[i] << " r: " << r[i] << endl;
	//sort sim x for plotting
	//std::sort(x,x+nDataPts); - don't need to since sampling from larger population

	


	//x and r are functions of t
	double x_approx[nSamples+1]; //filled by reference in CalculateCurve
	double r_approx[nSamples+1]; //filled by reference in CalculateCurve
	double t[nSamples+1];
	BezierCurve x_bc(x, nDataPts, t, nSamples);
	BezierCurve r_bc(r, nDataPts, t, nSamples);
	x_bc.CalculateCurve(x_approx);
	r_bc.CalculateCurve(r_approx);
	cout << "calculated Bezier curves for x and r - order " << nDataPts-2 << endl;


	if(viz){	
		gr = new TGraph(nSamples+1, r_approx, x_approx);
		
		gr->SetLineWidth(1);
		gr->SetLineColor(kRed-7);
		gr->SetMarkerStyle(15);

		cout << "Plots stored in: " << "test/" << fname << endl;
		TFile* f = TFile::Open(("test/"+fname).c_str(),"RECREATE");
		TGraph* gr_data = new TGraph(nDataPts, r, x);
		TLine* line = new TLine(0,0,1,1);
		TMultiGraph* mg = new TMultiGraph();		
		TLegend* leg = new TLegend(0.54,0.23,0.89,0.35);
		TCanvas* cv = new TCanvas("cv");
       	


		gr_data->SetMarkerSize(0.95);
       		gr_data->SetLineWidth(0);
		gr_data->SetMarkerStyle(20);
       		gr_data->SetMarkerColor(kViolet+7);
		
		mg->Add(gr_data);
		mg->Add(gr);
		
		leg->AddEntry(gr_data, "Original data/control points");
		leg->AddEntry(gr,"3rd order Bezier approximation");
		
		cv->cd();
		cv->SetRightMargin(0.09);
		cv->SetLeftMargin(0.15);
		cv->SetBottomMargin(0.15);
		

		mg->SetTitle("");
		mg->Draw("ALP");
		mg->GetXaxis()->SetTitle("r - rank space");
		mg->GetYaxis()->SetTitle("x - coord. space");
		leg->Draw("same");

		line->SetLineStyle(7);
		line->Draw("ALsame");

		TLatex l;
		l.SetNDC();
		l.SetTextAlign(11);
		l.SetTextSize(0.03);
		l.DrawLatex(0.18,0.92,"Bezier Curves to Independently Approximate x and r");
		
		f->cd();
		cv->Write();
		f->Close();
	}






}
