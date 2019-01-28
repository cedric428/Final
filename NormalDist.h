// Normal distribution lib
// dhan, 20170917

#ifndef NORMALDIST
#define NORMALDIST

class NormalDist{
public:
	double CDF(double x); // normal cdf
	double InverseCDF(double x); // inverse normal cdf
	double PDF(double x); // normal pdf
};


#endif