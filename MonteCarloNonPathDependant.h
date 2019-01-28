//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_MONTECARLONONPATHDEPENDANT_H
#define MONTECARLO_MONTECARLONONPATHDEPENDANT_H


#include <vector>
#include "InverseTransform.h"
#include "BlackScholes.h"

class MonteCarloNonPathDependant {
private:
    NormalGenerator *normalGenerator;
    std::vector<double> z = std::vector<double>(5120000);
    std::vector<double> s = std::vector<double>(5120000);
    std::vector<double> c = std::vector<double>(5120000);
    std::vector<double> deltaC = std::vector<double>(5120000);
    std::vector<double> vegaC = std::vector<double>(5120000);
    std::vector<double> p = std::vector<double>(5120000);
    std::vector<double> deltaP = std::vector<double>(5120000);
    std::vector<double> vegaP = std::vector<double>(5120000);

protected:
    double mean(std::vector<double> x, unsigned long long int n);
    double S;
    double K;
    double r;
    double q;
    double sigma;
    double T;

public:
    MonteCarloNonPathDependant(NormalGenerator *normalGenerator, double S, double K, double r, double q, double sigma,
                               double T);

    void runSimulation();

    double getCallValue(unsigned long long int N);

    double getCallDelta(unsigned long long int N);

    double getCallVega(unsigned long long int N);

    double getPutValue(unsigned long long int N);

    double getPutDelta(unsigned long long int N);

    double getPutVega(unsigned long long int N);

    virtual ~MonteCarloNonPathDependant();
};


#endif //MONTECARLO_MONTECARLONONPATHDEPENDANT_H
