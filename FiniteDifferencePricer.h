// Finite difference pricing on a fixed computational domain
// dhan 20171013

#ifndef FINITE_DIFF_H
#define FINITE_DIFF_H


#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>
using namespace std;

typedef Eigen :: VectorXd vec ;
typedef Eigen :: MatrixXd mat ;
enum StoppingCriterion { consecutive , residual };

class FiniteDifferencePricer{
public:
	// parameters
	double S,K,B,vol,r,q,T;
	double x_compute;
	double tau_final;
	double x_left;
	double x_right;
	int N_left;
	int N_right;
	unsigned int M; // number of time divisions
	unsigned int N; // number of space division, computed from M
	double dt;
	double dx;
	double alpha; // ratio of dt/dx^2
	double a;
	double b;
	vector<vector<double>> early_ex_premium;

	// grid
	vector<vector<double>> grid;
	vector<double> taxis;
	vector<double> xaxis;

	// system matrix for back euler and crank
	vector<vector<double>> A;

public:
	double Price;
	double Price2; // using interpolation of the heat equation value
	double Delta;
	double Gamma;
	double Theta;

public:
	FiniteDifferencePricer(double S,double K, double B, double vol, double r, double q, double T, unsigned int M, double alpha);
	virtual ~FiniteDifferencePricer();

	void SetBoundaryCondition(char euro);

	vector<vector<double>> FillGridUsingForwardEuler(char euro);
	vector<vector<double>> FillGridUsingBackwardEuler();
	vector<vector<double>> FillGridUsingCrankNicolson();
	vector<vector<double>> FillGridUsingCrankNicolsonUsingSOR(char euro);


	void ComputeResult();

	vector<vector<double>> FindEarlyExerciseDomain() const;

	void run(char euro, string fdm);

private:
	void printMatrix(const vector<vector<double>>& g) const;
	void printvec(const vector<double>& v) const;

	vector<double> addvector(const vector<double>& a, const vector<double> b, double factor = 1.0);
	vector<double> SOR( const mat & A , const vec & b , const vec & x_0 ,
const double tolerance , const StoppingCriterion criterion , const double w);
	
};


#endif