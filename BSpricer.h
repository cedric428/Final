// BlackScholes pricer
// dhan, 20170917

#ifndef BSPRICER_H
#define BSPRICER_H

class BSpricer{
public:
	double BSprice;
    double BSdelta;
    double BSgamma;
    double BStheta;
    double BSvega;

    double GetPrice(double S, double K, double T, 
    	double vol, double q, double r, char call);


};


#endif