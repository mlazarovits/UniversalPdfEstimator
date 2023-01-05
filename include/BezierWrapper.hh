#ifndef BezierWrapper_HH
#define BezierWrapper_HH

#include <iostream>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>

using std::string;

class BezierWrapper{
	public:
		BezierWrapper(){ };
		BezierWrapper(int nCP);
		virtual ~BezierWrapper();

		void Sample(unsigned long long seed = 123, double xmin = 1, double xmax = 0, double mu = 0.5, double sigma = 0.25);
		void CalculateCurve();
		void SetNTransforms(int n);
		void Transform(int n = 0); //external transformation interface
		void CalculateJacobian();
		void Visualize();
		void VisualizeTransform(double* x_transf);
	private:
		string _fname = "test.root";
		int _nDataPts; //number of control points/data points
		int _nSamples = 1000;//number of samples for t NOT including endpoints and x_global and r_global
		int _nTransforms = 0; //number of times to transform data
		void _transform(); //internal transformation function
		
		//arrays
		double* _x;
		double* _r;
		double* _x_global;
		double* _r_global;
		double* _x_approx; //filled by reference in CalculateCurve
		double* _r_approx; //filled by reference in CalculateCurve
		double* _t;

		//plotting
		TMultiGraph* _mg = new TMultiGraph();		
		TLegend* _leg = new TLegend(0.54,0.23,0.89,0.35);
		TCanvas* _cv = new TCanvas("cv");
};







#endif
