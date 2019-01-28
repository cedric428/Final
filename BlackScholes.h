//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_BLACKSCHOLES_H
#define MONTECARLO_BLACKSCHOLES_H

class BlackScholes {
public:
    BlackScholes(double S, double K, double r, double q, double sigma, double T);

    double getCallValue();

    double getCallDelta();

    double getCallVega();

    double getPutValue();

    double getPutDelta();

    double getPutVega();

    double getDAOValue(double B);

protected:
    double S;
    double K;
    double r;
    double q;
    double sigma;
    double T;

    double d1();

    double d2();

    // CDF of Standard normal.
    double N(double x);

    // PDF of Standard normal.
    double N_prime(double x);
};


#endif //MONTECARLO_BLACKSCHOLES_H
