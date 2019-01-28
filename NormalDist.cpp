// NormalDist.cpp
// dhan, 20170917
#include "NormalDist.h"
#include <cmath>
using namespace std;


double NormalDist::PDF(double x){
	// normal pdf
	double pi = atan(1) * 4;
    return (1/sqrt(2*pi)) * exp(-0.5 * x * x);
}

double NormalDist::CDF(double x){
	// normal cdf
	// Hastings approximation, maximum abs error less than 7.5*10e-8
	const double b1 = .31938153;
	const double b2 = -.356563782;
	const double b3 = 1.781477937;
	const double b4 = -1.821255978;
	const double b5 = 1.330274429;
	const double p = .2316419;
	const double c = 0.918938533204672;
	double a = fabs(x);
	double t = 1.0/(1.0+a*p);
	double s = (((((b5)*t+b4)*t+b3)*t+b2)*t+b1)*t;
	double y = s*exp(-0.5*x*x-c);
	y = (x>0)? 1-y:y;
	return y;

}

double NormalDist::InverseCDF(double u){
	// inverse cdf of normal
	// Beasley-Springer-Moro algo for inverse normal
	const double a0 = 2.50662823884;
    const double a1 = -18.61500062529;
    const double a2 = 41.39119773534;
    const double a3 = -25.44106049637;

    const double b0 = -8.47351093090;
    const double b1 = 23.08336743743;
    const double b2 = -21.06224101826;
    const double b3 = 3.13082909833;

    const double c0 = 0.3374754822726147;
    const double c1 = 0.9761690190917186;
    const double c2 = 0.1607979714918209;
    const double c3 = 0.0276438810333863;
    const double c4 = 0.0038405729373609;
    const double c5 = 0.0003951896511919;
    const double c6 = 0.0000321767881768;
    const double c7 = 0.0000002888167364;
    const double c8 = 0.0000003960315187;

    double y = u-.5;
    double r,x;
    if (fabs(y)<.42){
    	r = y*y;
    	x = y*(((a3*r+a2)*r+a1)*r+a0)/(((((b3*r+b2)*r+b1)*r+b0))*r+1);
    }
    else{
    	r = (y>0.0)? 1-u:u;
    	r = log(-log(r));
    	x = c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))));
    	x = (y<0.0)? -x:x;
    }
    return x;
}
