#ifndef NumDeriv_HH
#define NumDeriv_HH


class NumDeriv{
	public:
		NumDeriv(){ };
		NumDeriv(double* x, double* y, int n);
		virtual ~NumDeriv();
		void FiniteDiff(double* deriv);
		void MapToInterval(double new_min, double new_max);

	private:
		double* _x;
		double* _y;
		int _n;

};




#endif
