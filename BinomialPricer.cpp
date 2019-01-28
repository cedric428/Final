#include "BinomialPricer.h"
// dhan 20170902

#include<cmath>
#include<iostream>
#include<vector>

using namespace std;

BinomialPricer::BinomialPricer() {
	//ctor
}

BinomialPricer::~BinomialPricer() {
	//dtor
}

double
BinomialPricer::BinomialTree(int N, double S, double K, double T, double vol, double r, double q, int call, int amer) {
	// implement the plain Binomial tree
	// output the price of option
	// variable call = 1 means call, -1 means put
	// variable amer = 0 means it is euro put, 1 means amer put, 2 means amer put with variance reduction
	double dt = T / N;    // delta_t
	double u = exp(vol * sqrt(dt));   // u
	double d = 1 / u; // d
	double prn = (exp((r - q) * dt) - d) / (u - d);   // p_RN risk neutral probability
	vector<double> v(N + 1, 1.0);

	for (int i = 0; i <= N; ++i) {
		v[i] = max(0.0, call * (S * pow(u, N - i) * pow(d, i) - K)); // payoff function for all possible closing prices
	}

	double v20, v21, v22, v10, v11; // for greeks computation
	// Counting backward
	for (int k = N - 1; k >= 0; --k) {
		for (int i = 0; i <= k; ++i) {
			v[i] = exp(-r * dt) *
			       (prn * v[i] + (1 - prn) * v[i + 1]); // discounted expected value under risk-neutral prob measure
			if (amer) {  // if amer = 1 : american option
				v[i] = max(max(K - S * pow(u, k - i) * pow(d, i), 0.0), v[i]); // if this is American put
				// Amer call: v[i] = max(max(S*pow(u,k-i)*pow(d,i)-K,0.0),v[i]);
			}
		}
		if (k == 2) {
			v20 = v[0];
			v21 = v[1];
			v22 = v[2];
		}
		if (k == 1) {
			v10 = v[0];
			v11 = v[1];
		}
	}
	delta = (v10 - v11) / (S * u - S * d);
	gamma = (((v20 - v21) / (S * u * u - S * u * d)) - ((v21 - v22) / (S * u * d - S * d * d))) /
	        ((S * u * u - S * d * d) / 2);
	theta = (v21 - v[0]) / (2 * dt);
	price = v[0];

	if (amer == 2) {
		//variance reduction for American put
		double raw_price = price, raw_delta = delta, raw_gamma = gamma, raw_theta = theta; // price without reduction
		BinomialTree(N, S, K, T, vol, r, q, -1, 0); // update greeks with euro put
		BSvalue(S, K, T, vol, r, q, -1); // update BS euro put
		price = raw_price - (price - BSprice);
		delta = raw_delta - (delta - BSdelta);
		gamma = raw_gamma - (gamma - BSgamma);
		theta = raw_theta - (theta - BStheta);
	}

	return price;
}

double BinomialPricer::AvgBinomialTree(int N, double S, double K, double T, double vol, double r, double q, int call,
                                       int amer) {
	// implement the avg Binomial tree
	// variable amer = 0 means it is euro put, 1 means amer put, 2 means amer put with variance reduction

	double price1 = BinomialTree(N, S, K, T, vol, r, q, call, amer);
	double delta1 = delta, gamma1 = gamma, theta1 = theta;
	double price2 = BinomialTree(N + 1, S, K, T, vol, r, q, call, amer);

	delta = (delta + delta1) / 2;
	gamma = (gamma + gamma1) / 2;
	theta = (theta + theta1) / 2;
	return 0.5 * (price1 + price2);
}

double
BinomialPricer::BinomialBlackScholes(int N, double S, double K, double T, double vol, double r, double q, int call,
                                     int amer) {
	// implement the Binomial tree with Black Scholes initialization
	// output the price of option
	// variable call = 1 means call, -1 means put
	double dt = T / N;
	double u = exp(vol * sqrt(dt));
	double d = 1 / u;
	double prn = (exp((r - q) * dt) - d) / (u - d);
	vector<double> v(N + 1, 1.0);

	for (int i = 0; i <= N - 1; ++i) {
		double tmp = 0;
		if ((call == -1) && (amer > 0)) tmp = K - S * pow(u, N - 1 - i) * pow(d, i); // if this is amer put
		v[i] = max(tmp, BSvalue(S * pow(u, N - 1 - i) * pow(d, i), K, dt, vol, r, q,
		                        call)); // init N-1 time value with the BS value, THIS WILL CHANGE BS VALUE MEMBER,CAUTION
	}

	double v20, v21, v22, v10, v11; // for greeks computation
	for (int k = N - 2; k >= 0; --k) {
		for (int i = 0; i <= k; ++i) {
			v[i] = exp(-r * dt) *
			       (prn * v[i] + (1 - prn) * v[i + 1]); // compute discounted value under risk-neutral prob measure
			if (amer) {
				v[i] = max(max(K - S * pow(u, k - i) * pow(d, i), 0.0), v[i]); // if this is American put
			}
		}
		if (k == 2) {
			v20 = v[0];
			v21 = v[1];
			v22 = v[2];
		}
		if (k == 1) {
			v10 = v[0];
			v11 = v[1];
		}
	}

	delta = (v10 - v11) / (S * u - S * d);
	gamma = (((v20 - v21) / (S * u * u - S * u * d)) - ((v21 - v22) / (S * u * d - S * d * d))) /
	        ((S * u * u - S * d * d) / 2);
	theta = (v21 - v[0]) / (2 * dt);
	price = v[0];
	BSvalue(S, K, T, vol, r, q, call); // set BSVALUE back

	if (amer == 2) {
		//variance reduction for American put
		double raw_price = price, raw_delta = delta, raw_gamma = gamma, raw_theta = theta; // price without reduction
		BinomialTree(N, S, K, T, vol, r, q, -1, 0); // update greeks with euro put
		BSvalue(S, K, T, vol, r, q, -1); // update BS euro put
		price = raw_price - (price - BSprice);
		delta = raw_delta - (delta - BSdelta);
		gamma = raw_gamma - (gamma - BSgamma);
		theta = raw_theta - (theta - BStheta);
	}

	return price;
}

double BinomialPricer::BinomialBlackScholesWithRE(int N, double S, double K, double T, double vol, double r, double q,
                                                  int call, int amer) {
	double price1 = BinomialBlackScholes(N, S, K, T, vol, r, q, call, amer);
	double delta1 = delta, gamma1 = gamma, theta1 = theta;
	double price2 = BinomialBlackScholes(N / 2, S, K, T, vol, r, q, call, amer);

	delta = 2 * delta1 - delta;
	gamma = 2 * gamma1 - gamma;
	theta = 2 * theta1 - theta;
	return 2 * price1 - price2;
}

//double BinomialPricer::CDFnormal(double x) {
//    return 1 / (1 + exp(-0.07056 * pow(x, 3) - 1.5976 * x));
//}

//cdf for normal distribution, copied from online sources
double BinomialPricer::CDFnormal(double x) {
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

double BinomialPricer::PDFnormal(double x) {
	double pi = atan(1) * 4;
	return (1 / sqrt(2 * pi)) * exp(-0.5 * x * x);
}

double BinomialPricer::BSvalue(double S, double K, double T, double vol, double r, double q, int call) {
	// compute BS model price
	double disc = exp(-r * T);
	double PVF = S * exp((-q) * T);
	double d1 = log(PVF / K / disc) / vol / sqrt(T) + vol * sqrt(T) / 2;
	double d2 = log(PVF / K / disc) / vol / sqrt(T) - vol * sqrt(T) / 2;
	if (call == 1) {
		// call
		BSprice = PVF * CDFnormal(d1) - K * disc * CDFnormal(d2);
		BSdelta = exp(-q * T) * CDFnormal(d1);
		BStheta = -exp(-q * T) * S * PDFnormal(d1) * vol / 2 / sqrt(T) - r * K * exp(-r * T) * CDFnormal(d2) +
		          q * S * exp(-q * T) * CDFnormal(d1);
	} else {
		BSprice = K * disc * CDFnormal(-d2) - PVF * CDFnormal(-d1);
		BSdelta = -exp(-q * T) * CDFnormal(-d1);
		BStheta = -exp(-q * T) * S * PDFnormal(d1) * vol / 2 / sqrt(T) + r * K * exp(-r * T) * CDFnormal(-d2) -
		          q * S * exp(-q * T) * CDFnormal(-d1);
	}
	BSgamma = exp(-q * T) * PDFnormal(d1) / S / vol / sqrt(T);

	return BSprice;
}

double
BinomialPricer::DOCBinomail(int N, double S, double K, double T, double vol, double r, double q, int call, int amer,
                            double L) {
	// implement the plain Binomial tree
	// output the price of option
	// variable call = 1 means call, -1 means put
	// variable amer = 0 means it is euro put, 1 means amer put, 2 means amer put with variance reduction
	double dt = T / N;    // delta_t
	double u = exp(vol * sqrt(dt));   // u
	double d = 1 / u; // d
	double prn = (exp((r - q) * dt) - d) / (u - d);   // p_RN risk neutral probability
	vector<double> v(N + 1, 1.0);
	double S_final;
	for (int i = 0; i <= N; ++i) {
		S_final = S * pow(u, N - i) * pow(d, i);
		if (S_final > L){
			v[i] = max(0.0, call * (S * pow(u, N - i) * pow(d, i) - K)); // payoff function for all possible closing prices
		}else{
			v[i] =0;
		}

	}

	double v20, v21, v22, v10, v11; // for greeks computation
	// Counting backward
	double S_i =0;
	for (int k = N - 1; k >= 0; --k) {
		for (int i = 0; i <= k; ++i) {
			S_i = S * pow(u, k - i) * pow(d, i);
			if(S_i > L){
				v[i] = exp(-r * dt) *
				       (prn * v[i] + (1 - prn) * v[i + 1]); // discounted expected value under risk-neutral prob measure

			}else{
				v[i] = 0;
			}


			if (amer) {  // if amer = 1 : american option
				v[i] = max(max(K - S * pow(u, k - i) * pow(d, i), 0.0), v[i]); // if this is American put
				// Amer call: v[i] = max(max(S*pow(u,k-i)*pow(d,i)-K,0.0),v[i]);
			}
		}
		if (k == 2) {
			v20 = v[0];
			v21 = v[1];
			v22 = v[2];
		}
		if (k == 1) {
			v10 = v[0];
			v11 = v[1];
		}
	}
	delta = (v10 - v11) / (S * u - S * d);
	gamma = (((v20 - v21) / (S * u * u - S * u * d)) - ((v21 - v22) / (S * u * d - S * d * d))) /
	        ((S * u * u - S * d * d) / 2);
	theta = (v21 - v[0]) / (2 * dt);
	price = v[0];

	if (amer == 2) {
		//variance reduction for American put
		double raw_price = price, raw_delta = delta, raw_gamma = gamma, raw_theta = theta; // price without reduction
		BinomialTree(N, S, K, T, vol, r, q, -1, 0); // update greeks with euro put
		BSvalue(S, K, T, vol, r, q, -1); // update BS euro put
		price = raw_price - (price - BSprice);
		delta = raw_delta - (delta - BSdelta);
		gamma = raw_gamma - (gamma - BSgamma);
		theta = raw_theta - (theta - BStheta);
	}

	return price;
}
