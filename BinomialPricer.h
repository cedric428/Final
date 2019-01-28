#ifndef BINOMIALPRICER_H
#define BINOMIALPRICER_H
// Binomial tree method for pricing options
// binomial tree
// average binomial tree
// binomial BS
// binomial BS with Richardson Extrapolation

// dhan, 20170902

class BinomialPricer
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
        BinomialPricer();
        virtual ~BinomialPricer();
        double BinomialTree(int N, double S, double K, double T, double vol, double r, double q, int call, int amer);
        double AvgBinomialTree(int N, double S, double K, double T, double vol, double r, double q, int call, int amer);
        double BinomialBlackScholes(int N, double S, double K, double T, double vol, double r, double q, int call, int amer);
        double BinomialBlackScholesWithRE(int N, double S, double K, double T, double vol, double r, double q, int call, int amer);
        double BSvalue(double S, double K, double T, double vol, double r, double q, int call);
		double DOCBinomail(int N, double S, double K, double T, double vol, double r, double q, int call, int amer, double L);

    protected:

    private:
        double CDFnormal(double x);
        double PDFnormal(double x);

};

#endif // BINOMIALPRICER_H
