//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_MONTECARLOPATHDEPENDANT_H
#define MONTECARLO_MONTECARLOPATHDEPENDANT_H


#include <vector>
#include "NormalGenerator.h"
#include "BlackScholes.h"

class MonteCarloPathDependant {
private:
    unsigned long long int m;
    NormalGenerator *normalGenerator;
    std::vector<double> z = std::vector<double>(5120000);
    std::vector<double> s = std::vector<double>(5120000);
    std::vector<double> v = std::vector<double>();

protected:
    double mean(std::vector<double> x, unsigned long long int n);
    double S;
    double K;
    double r;
    double q;
    double sigma;
    double T;
    double B;

public:
    MonteCarloPathDependant(NormalGenerator *normalGenerator, double S, double K, double r,
                                double q, double sigma, double T, double B);

    void runSimulation(unsigned long long int m);

    double getDAOPrice(unsigned long long int N);

    virtual ~MonteCarloPathDependant();
};


#endif //MONTECARLO_MONTECARLOPATHDEPENDANT_H
