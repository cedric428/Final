//
// Created by Harish Reddy on 9/17/2017.
//

#include <cmath>
#include <iostream>
#include "MonteCarloPathDependant.h"

// For down-and-out call options.
void MonteCarloPathDependant::runSimulation(unsigned long long int m) {
    this->m = m;
    unsigned int i, j;
    const unsigned long long int n = 5120000 / m;
    const double deltaT = T / m;
    v.clear();
    z = normalGenerator->getNormals(z.size());
    for (i = 0; i < n; i++) {
        bool barrierCrossed = false;
        for (j = 0; j < m; j++) {
            s[i * m + j] = (j == 0 ? S : s[i * m + j - 1]) *
                           exp((r - q - sigma * sigma / 2) * deltaT + sigma * sqrt(deltaT) * z[i * m + j]);
            if (s[i * m + j] <= B) {
                barrierCrossed = true;  //
            }
        }
        double p = exp(-r * T) * std::max(s[i * m + m - 1] - K, 0.0);
        v.push_back(barrierCrossed ? 0 : p);
    }
}

double MonteCarloPathDependant::getDAOPrice(unsigned long long int N) {
    unsigned long long int n = N / m;
    return mean(v, n);
}

double MonteCarloPathDependant::mean(std::vector<double> x, unsigned long long int n) {
    unsigned int i;
    double sum = 0;
    for (i = 0; i < n; i++) {
        sum += x.at(i);
    }
    return sum / n;
}

MonteCarloPathDependant::MonteCarloPathDependant(NormalGenerator *normalGenerator, double S, double K, double r,
                                                 double q, double sigma, double T, double B)
        : normalGenerator(normalGenerator),
          S(S), K(K), r(r), q(q),
          sigma(sigma), T(T), B(B) {
}

MonteCarloPathDependant::~MonteCarloPathDependant() {
    delete normalGenerator;
}
