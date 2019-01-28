//
// Created by Harish Reddy on 9/17/2017.
//

#include <cmath>
#include "BlackScholes.h"


BlackScholes::BlackScholes(double S, double K, double r, double q, double sigma, double T) : S(S), K(K), r(r), q(q),
                                                                                             sigma(sigma), T(T) {}

double BlackScholes::getCallValue() {
    return S * exp(-q * T) * N(d1()) - K * exp(-r * T) * N(d2());
}

double BlackScholes::getCallDelta() {
    return exp(-q * T) * N(d1());
}

double BlackScholes::getCallVega() {
    return S * exp(-q * T) * N_prime(d1()) * sqrt(T);
}

double BlackScholes::getPutValue() {
    return -S * exp(-q * T) * N(-d1()) + K * exp(-r * T) * N(-d2());
}

double BlackScholes::getPutDelta() {
    return -exp(-q * T) * N(-d1());
}

double BlackScholes::getPutVega() {
    return S * exp(-q * T) * N_prime(d1()) * sqrt(T);
}

double BlackScholes::d1() {
    return (log(S / K) + (r - q + sigma * sigma / 2) * T) / (sigma * sqrt(T));
}

double BlackScholes::d2() {
    return (log(S / K) + (r - q - sigma * sigma / 2) * T) / (sigma * sqrt(T));
}

double BlackScholes::getDAOValue(double B) {
    double a = (r - q) / (sigma * sigma) - 0.5;
    BlackScholes blackScholes = BlackScholes(B * B / S, K, r, q, sigma, T);
    return this->getCallValue() - pow(B / S, 2 * a) * blackScholes.getCallValue();
}

//cdf for normal distribution, copied from online sources
double BlackScholes::N(double x) {
    // calculate \sqrt{2\pi} upfront once
    static const double RT2PI = sqrt(4.0 * acos(0.0));
    // calculate 10/\sqrt{2} upfront once
    static const double SPLIT = 10. / sqrt(2);
    static const double a[] = {220.206867912376, 221.213596169931, 112.079291497871, 33.912866078383, 6.37396220353165,
                               0.700383064443688, 3.52624965998911e-02};
    static const double b[] = {440.413735824752, 793.826512519948, 637.333633378831, 296.564248779674, 86.7807322029461,
                               16.064177579207, 1.75566716318264, 8.83883476483184e-02};

    const double z = fabs(x);
    // Now N(x) = 1 - N(-x) = 1-\sqrt{2\pi}N'(x)\frac{P(x)}{Q(x)}
    //  so N(-x) = \sqrt{2\pi}N'(x)\frac{P(x)}{Q(x)}
    // now let \sqrt{2\pi}N'(z)\frac{P(x)}{Q(z)} = Nz
    // Therefore we have
    //     Nxm = N(x) = \sqrt{2\pi}N'(z)\frac{P(x)}{Q(z)} = Nz if x<0
    //     Nxp = N(x) = 1 - \sqrt{2\pi}N'(z)\frac{P(x)}{Q(z)} = 1-Nz if x>=0
    double Nz = 0.0;

    // if z outside these limits then value effectively 0 or 1 for machine precision
    if (z <= 37.0) {
        // NDash = N'(z) * sqrt{2\pi}
        const double NDash = exp(-z * z / 2.0) / RT2PI;
        if (z < SPLIT) {
            // here Pz = P(z) is a polynomial
            const double Pz = (((((a[6] * z + a[5]) * z + a[4]) * z + a[3]) * z + a[2]) * z + a[1]) * z + a[0];
            // and Qz = Q(z) is a polynomial
            const double Qz =
                    ((((((b[7] * z + b[6]) * z + b[5]) * z + b[4]) * z + b[3]) * z + b[2]) * z + b[1]) * z + b[0];
            // use polynomials to calculate N(z)  = \sqrt{2\pi}N'(x)\frac{P(x)}{Q(x)}
            Nz = RT2PI * NDash * Pz / Qz;
        } else {
            // implement recurrence relation on F_4(z)
            const double F4z = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
            // use polynomials to calculate N(z), note here that Nz = N' / F
            Nz = NDash / F4z;
        }
    }

    //
    return x >= 0.0 ? 1 - Nz : Nz;
}

//pdf for normal distribution
double BlackScholes::N_prime(double x) {
    return (1 / sqrt(2 * M_PI)) * exp(-0.5 * x * x);
}