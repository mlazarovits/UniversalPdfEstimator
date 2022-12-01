#include "RandomSample.hh"
#include <iostream>

using std::cout;
using std::endl;

RandomSample::RandomSample(unsigned long long seed){
	_b = 4101842887655102017LL;
	_c = 1;
	_a = seed ^ _b;

	int64();
	_b = _a;
	int64();
	_c = _b;
	int64();



}


RandomSample::~RandomSample() { }

//returns pseudo-random 64bit integer
unsigned long long RandomSample::int64(){
	_a = _a * 2862933555777941757LL + 7046029254386353087LL;
	_b ^= _b >> 17;
	_b ^= _b << 31;
	_b ^= _b >> 8;
	_c = 4294957665U*(_c & 0xffffffff) + (_c >> 32);

	unsigned long long x = _a ^ (_a << 21);
	x ^= x >> 35;
	x ^= x << 4;
	return (x + _b) ^ _c;
}


//return random 32bit integer
unsigned int RandomSample::int32(){
	return (unsigned int)int64();
}

//return random double between (0,1) (uniform distribution)
double RandomSample::rand(){
	return 5.42101086242752217E-20 * int64();
}

//flat distribution scaled to Gaussian max
double RandomSample::FlatGausScaled(){
	return 1./sqrt(2.*acos(-1));
}

//normal distribution with specified mean and sigma
double RandomSample::Gaussian(double x, double mu, double sigma){
	double x_scaled = x - mu;
	return 1./(sigma*sqrt(2.*acos(-1)))*exp(-(x_scaled*x_scaled)/(2.*sigma*sigma));
}

//get random x value according to flat distribuion
double RandomSample::SampleFlat(){
	return _xmin + (_xmax - _xmin)*rand();
}

void RandomSample::SetRange(double xmin, double xmax){
	if(xmin >= xmax) return;
	_xmin = xmin;
	_xmax = xmax;
}

//return random double (xmin, xmax) according to Gaussian distribuion
double RandomSample::SampleGaussian(double mean, double sigma, int Nsample, double* samples){
	if(sigma < 0){
		cout << "Please input valid sigma > 0" << endl;
		return 0.;
	}

	double X, ran, R;
	//double samples[Nsample];
	int Ntrial = 0;
	for(int i = 0; i < Nsample; i++){
		Ntrial += 1;
		X = SampleFlat();
		R = Gaussian(X,mean,sigma)/FlatGausScaled();
		ran = rand();
		if(ran > R){ // reject if "die thrown" (random number) is above distribution
			i -= 1; // decrease i and try again
			continue;
		} else // accept if "die" is below (within) distribution
			samples[i] = X;
	}
	return X;
}




