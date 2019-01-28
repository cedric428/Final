//
// Created by Harish Reddy on 9/17/2017.
//

#include <cmath>
#include <iostream>
#include "MonteCarloNonPathDependant.h"

void MonteCarloNonPathDependant::runSimulation() {
    unsigned long long int i;
    z = normalGenerator->getNormals(z.size());
    for (i = 0; i < z.size(); i++) {
        s[i] = S * exp((r - q - sigma * sigma / 2) * T + sigma * sqrt(T) * z[i]);
        c[i] = exp(-r * T) * std::max(s[i] - K, 0.0);
        deltaC[i] = (s[i] > K) ? (exp(-r * T) * s[i] / S) : 0;
        vegaC[i] = (s[i] > K) ? (exp(-r * T) * s[i] * (-sigma * T + sqrt(T) * z[i])) : 0;

        p[i] = exp(-r * T) * std::max(-s[i] + K, 0.0);;
        deltaP[i] = (s[i] < K) ? -(exp(-r * T) * s[i] / S) : 0;
        vegaP[i] = (s[i] < K) ? -(exp(-r * T) * s[i] * (-sigma * T + sqrt(T) * z[i])) : 0;
    }
}

double MonteCarloNonPathDependant::getCallValue(unsigned long long int N) {
    return mean(c, N);
}

double MonteCarloNonPathDependant::getCallDelta(unsigned long long int N) {
    return mean(deltaC, N);
}

double MonteCarloNonPathDependant::getCallVega(unsigned long long int N) {
    return mean(vegaC, N);
}

double MonteCarloNonPathDependant::getPutValue(unsigned long long int N) {
    return mean(p, N);
}

double MonteCarloNonPathDependant::getPutDelta(unsigned long long int N) {
    return mean(deltaP, N);
}

double MonteCarloNonPathDependant::getPutVega(unsigned long long int N) {
    return mean(vegaP, N);
}

double MonteCarloNonPathDependant::mean(std::vector<double> x, unsigned long long int n) {
    unsigned long long int i;
    double sum = 0;
    for (i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum / n;
}

MonteCarloNonPathDependant::MonteCarloNonPathDependant(NormalGenerator *normalGenerator, double S, double K, double r,
                                                       double q, double sigma, double T) : normalGenerator(
        normalGenerator), S(S), K(K), r(r), q(q), sigma(sigma), T(T) {}


MonteCarloNonPathDependant::~MonteCarloNonPathDependant() {
    delete normalGenerator;
}
