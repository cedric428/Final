// BSpricer.cpp
//dhan, 20170917

#include <cmath>

#include "BSpricer.h"
#include "NormalDist.h"


double BSpricer::GetPrice(double S, double K, double T, double vol, double q, double r, char call) {
	// compute BS model price
	NormalDist* norm = new NormalDist;
	double disc = exp(-r*T);
	double PVF = S * exp((-q)*T);
	double d1 = log(PVF / K / disc) / vol / sqrt(T) + vol*sqrt(T) / 2;
	double d2 = log(PVF / K / disc) / vol / sqrt(T) - vol*sqrt(T) / 2;
	if (call == 'c') {
		// call
		BSprice = PVF*norm->CDF(d1) - K * disc * norm->CDF(d2);
		BSdelta = exp(-q*T) * norm->CDF(d1);
		BStheta = -exp(-q*T)*S*norm->PDF(d1) *vol/2/sqrt(T) - r*K*exp(-r*T)*norm->CDF(d2) + q*S*exp(-q*T)*norm->CDF(d1);
	}
	else {
		BSprice = K * disc * norm->CDF(-d2) - PVF*norm->CDF(-d1);
		BSdelta = -exp(-q*T) * norm->CDF(-d1);
		BStheta = -exp(-q*T)*S*norm->PDF(d1) *vol/2/sqrt(T) + r*K*exp(-r*T)*norm->CDF(-d2) - q*S*exp(-q*T)*norm->CDF(-d1);
	}
	BSgamma = exp(-q*T) * norm->PDF(d1) / S / vol / sqrt(T);
	BSvega = exp(-q*T) * norm->PDF(d1) * S * sqrt(T);
	delete norm;
	return BSprice;
}