#ifndef TRINOMIALPRICER_H
#define TRINOMIALPRICER_H
// Trinomial tree method for pricing options
// Trinomial tree
// Trinomial BS
// Trinomial BS with Richardson Extrapolation

// dhan, 20170908

class TrinomialPricer
{
    public:
        // greeks
        double price;
        double delta;
        double gamma;
        double theta;
        double vega;

        double BSprice;
        double BSdelta;
        double BSgamma;
        double BStheta;
        double BSvega;

        // different pricing functions
        TrinomialPricer();
        virtual ~TrinomialPricer();
        double TrinomialTree(int N, double S, double K, double T, double vol, double r, double q, int call, int amer);
        double TrinomialBlackScholes(int N, double S, double K, double T, double vol, double r, double q, int call, int amer);
        double TrinomialBlackScholesWithRE(int N, double S, double K, double T, double vol, double r, double q, int call, int amer);
        double BSvalue(double S, double K, double T, double vol, double r, double q, int call);
		double
		TrinomialDoubleBarrier(int N, double S, double K, double T, double vol, double r, double q, double L, double U,
				                       double rebate, int amer, int call);

    protected:

    private:
        double CDFnormal(double x);
        double PDFnormal(double x);

};


#endif // TRINOMIALPRICER_H
