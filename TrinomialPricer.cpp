#include "TrinomialPricer.h"
// dhan 20170909

#include<cmath>
#include<iostream>
#include<vector>

using namespace std;

TrinomialPricer::TrinomialPricer() {
	//ctor
}

TrinomialPricer::~TrinomialPricer() {
	//dtor
}


double TrinomialPricer::TrinomialTree(int N, double S, double K, double T, double vol, double r, double q, int call,
                                      int amer) {
	// implement the plain Trinomial tree
	// output the price of option
	// variable call = 1 means call, -1 means put
	// variable amer = 0 means it is euro put, 1 means amer put, 2 means amer put with variance reduction
	// amer = 3 means Up–and–Out–Down–and–Out–Call
	double dt = T / N;
	double u = exp(vol * sqrt(3 * dt));
	double m = 1.0;
	double d = 1 / u;
	double pu = 1.0 / 6.0 + (r - q - vol * vol * 0.5) * sqrt(dt / 12. / vol / vol);
	double pm = 2.0 / 3.0;
	double pd = 1.0 / 6.0 - (r - q - vol * vol * 0.5) * sqrt(dt / 12. / vol / vol);
	vector<double> v(2 * N + 1, 1.0);
	// 上面不用动

	for (int i = 0; i <= 2 * N; ++i) {
		v[i] = max(0.0, call * (S * pow(u, N - i) - K)); // payoff function for all possible closing prices
	}
	// underlying asset evolution
	double v20, v22, v24, v10, v11, v12; // for greeks computation

	for (int k = N - 1; k >= 0; --k) {
		for (int i = 0; i <= 2 * k; ++i) {
			v[i] = exp(-r * dt) * (pu * v[i] + pm * v[i + 1] +
			                       pd * v[i + 2]); // compute discounted value under risk-neutral prob measure

			if (amer ==1) {
				v[i] = max(max(K - S * pow(u, k - i), 0.0), v[i]); // if this is American put
			}
		}
		if (k == 2) {
			v20 = v[0];
			v22 = v[2];
			v24 = v[4];
		}
		if (k == 1) {
			v10 = v[0];
			v11 = v[1];
			v12 = v[2];
		}
	}
	delta = (v10 - v12) / (S * u - S * d);
	gamma = (((v20 - v22) / (S * u * u - S)) - ((v22 - v24) / (S - S * d * d))) / (S * u - S * d);
	theta = (v11 - v[0]) / (dt);
	price = v[0];

	if (amer == 2) {
		//variance reduction for American put
		double raw_price = price, raw_delta = delta, raw_gamma = gamma, raw_theta = theta; // price without reduction
		TrinomialTree(N, S, K, T, vol, r, q, -1, 0); // update greeks with euro put
		BSvalue(S, K, T, vol, r, q, -1); // update BS euro put
		price = raw_price - (price - BSprice);
		delta = raw_delta - (delta - BSdelta);
		gamma = raw_gamma - (gamma - BSgamma);
		theta = raw_theta - (theta - BStheta);
	}

	return price;
}


double
TrinomialPricer::TrinomialBlackScholes(int N, double S, double K, double T, double vol, double r, double q, int call,
                                       int amer) {
	// implement the Binomial tree with Black Scholes initialization
	// output the price of option
	// variable call = 1 means call, -1 means put
	double dt = T / N;
	double u = exp(vol * sqrt(3 * dt));
	double m = 1.0;
	double d = 1 / u;
	double pu = 1.0 / 6.0 + (r - q - vol * vol * 0.5) * sqrt(dt / 12 / vol / vol);
	double pm = 2.0 / 3.0;
	double pd = 1.0 / 6.0 - (r - q - vol * vol * 0.5) * sqrt(dt / 12 / vol / vol);
	vector<double> v(2 * N + 1, 1.0);

	for (int i = 0; i <= 2 * N - 2; ++i) {
		double tmp = 0;
		if ((call == -1) && (amer > 0)) tmp = K - S * pow(u, N - 1 - i); // if this is amer put
		v[i] = max(tmp, BSvalue(S * pow(u, N - 1 - i), K, dt, vol, r, q,
		                        call)); // init N-1 time value with the BS value, THIS WILL CHANGE BS VALUE MEMBER,CAUTION
	}

	double v20, v22, v24, v10, v11, v12; // for greeks computation
	for (int k = N - 2; k >= 0; --k) {
		for (int i = 0; i <= 2 * k; ++i) {
			v[i] = exp(-r * dt) * (pu * v[i] + pm * v[i + 1] +
			                       pd * v[i + 2]); // compute discounted value under risk-neutral prob measure
			if (amer) {
				v[i] = max(max(K - S * pow(u, k - i), 0.0), v[i]); // if this is American put
			}
		}
		if (k == 2) {
			v20 = v[0];
			v22 = v[2];
			v24 = v[4];
		}
		if (k == 1) {
			v10 = v[0];
			v11 = v[1];
			v12 = v[2];
		}
	}
	delta = (v10 - v12) / (S * u - S * d);
	gamma = (((v20 - v22) / (S * u * u - S)) - ((v22 - v24) / (S - S * d * d))) / (S * u - S * d);
	theta = (v11 - v[0]) / (dt);
	price = v[0];
	BSvalue(S, K, T, vol, r, q, call); // set BSVALUE back

	if (amer == 2) {
		//variance reduction for American put
		double raw_price = price, raw_delta = delta, raw_gamma = gamma, raw_theta = theta; // price without reduction
		TrinomialBlackScholes(N, S, K, T, vol, r, q, -1, 0); // update greeks with euro put
		BSvalue(S, K, T, vol, r, q, -1); // update BS euro put
		price = raw_price - (price - BSprice);
		delta = raw_delta - (delta - BSdelta);
		gamma = raw_gamma - (gamma - BSgamma);
		theta = raw_theta - (theta - BStheta);
	}

	return price;
}

double TrinomialPricer::TrinomialBlackScholesWithRE(int N, double S, double K, double T, double vol, double r, double q,
                                                    int call, int amer) {
	double price1 = TrinomialBlackScholes(N, S, K, T, vol, r, q, call, amer);
	double delta1 = delta, gamma1 = gamma, theta1 = theta;
	double price2 = TrinomialBlackScholes(N / 2, S, K, T, vol, r, q, call, amer);

	delta = 2 * delta1 - delta;
	gamma = 2 * gamma1 - gamma;
	theta = 2 * theta1 - theta;
	return 2 * price1 - price2;
}

//double TrinomialPricer::CDFnormal(double x) {
//    return 1 / (1 + exp(-0.07056 * pow(x, 3) - 1.5976 * x));
//}

//cdf for normal distribution, copied from online sources
double TrinomialPricer::CDFnormal(double x) {
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

double TrinomialPricer::PDFnormal(double x) {
	double pi = atan(1) * 4;
	return (1 / sqrt(2 * pi)) * exp(-0.5 * x * x);
}

double TrinomialPricer::BSvalue(double S, double K, double T, double vol, double r, double q, int call) {
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
TrinomialPricer::TrinomialDoubleBarrier(int N, double S, double K, double T, double vol, double r, double q, double L, double U,
                                        double rebate, int amer, int call) {
	// implement the plain Trinomial tree
	// output the price of option
	// variable call = 1 means call, -1 means put
	// variable amer = 0 means it is euro put, 1 means amer put, 2 means amer put with variance reduction
	// amer = 3 means Up–and–Out–Down–and–Out–Call
	double dt = T / N;
	double u = exp(vol * sqrt(3 * dt));
	double m = 1.0;
	double d = 1 / u;
	double pu = 1.0 / 6.0 + (r - q - vol * vol * 0.5) * sqrt(dt / 12. / vol / vol);
	double pm = 2.0 / 3.0;
	double pd = 1.0 / 6.0 - (r - q - vol * vol * 0.5) * sqrt(dt / 12. / vol / vol);
	vector<double> v(2 * N + 1, 1.0);
	// 上面不用动


	//Calculate the payoff of the option at maturity
	double S_final;
	for (int i = 0; i <= 2 * N; ++i) {
		// payoff function for all possible closing prices
		S_final = S * pow(u, N - i);
		if (S_final<U && S_final>L){
			v[i] = max(0.0, call * (S * pow(u, N - i) - K));
		}else{
			v[i] = rebate;
		}
	}


	double v20, v22, v24, v10, v11, v12; // for greeks computation

	double S_i = 0;
	for (int k = N - 1; k >= 0; --k) {
		for (int i = 0; i <= 2 * k; ++i) {
			S_i = S * pow(u, k - i);
			if (S_i >L&& S_i <U){
				v[i] = exp(-r * dt) * (pu * v[i] + pm * v[i + 1] +
				                       pd * v[i + 2]); // compute discounted value under risk-neutral prob measure
			}else{
				v[i] = exp(-r * dt*(N - k +i ))*rebate;   // N-k+i
			}


			if (amer ==1) {
				v[i] = max(max(K - S * pow(u, k - i), 0.0), v[i]); // if this is American put
			}
		}
		if (k == 2) {
			v20 = v[0];
			v22 = v[2];
			v24 = v[4];
		}
		if (k == 1) {
			v10 = v[0];
			v11 = v[1];
			v12 = v[2];
		}
	}
	delta = (v10 - v12) / (S * u - S * d);
	gamma = (((v20 - v22) / (S * u * u - S)) - ((v22 - v24) / (S - S * d * d))) / (S * u - S * d);
	theta = (v11 - v[0]) / (dt);
	price = v[0];

	if (amer == 2) {
		//variance reduction for American put
		double raw_price = price, raw_delta = delta, raw_gamma = gamma, raw_theta = theta; // price without reduction
		TrinomialTree(N, S, K, T, vol, r, q, -1, 0); // update greeks with euro put
		BSvalue(S, K, T, vol, r, q, -1); // update BS euro put
		price = raw_price - (price - BSprice);
		delta = raw_delta - (delta - BSdelta);
		gamma = raw_gamma - (gamma - BSgamma);
		theta = raw_theta - (theta - BStheta);
	}

	return price;
}
