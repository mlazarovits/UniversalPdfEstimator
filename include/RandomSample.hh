#ifndef RandomSample_HH
#define RandomSample_HH

#include <cmath>

class RandomSample{
	public:
		RandomSample(unsigned long long seed);
		virtual ~RandomSample();
		unsigned long long int64();
		unsigned int int32();
		double rand();
		double FlatGausScaled();
		double Gaussian(double x, double mu, double sigma);
		double SampleFlat();
		void SetRange(double xmin, double xmax);
		double SampleGaussian(double mean, double sigma, int Nsample, double* samples);
		
		double _xmax = 5;
		double _xmin = -5;

	private:
		unsigned long long _a;
		unsigned long long _b;
		unsigned long long _c;





};







#endif
