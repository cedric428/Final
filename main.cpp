// MTH9821 HW 1
// dhan 20170902
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <memory>
#include <iomanip>
#include "FiniteDifferencePricer.h"
#include "BSpricer.h"
#include "LinearSystemSolver.h"
#include "BinomialPricer.h"
#include "TrinomialPricer.h"

int precision = 18;
using namespace std;

/**
 * TrinomialPricer.cpp TrinomialPricer.h
        AcceptanceRejection.cpp AcceptanceRejection.h
        AcceptanceRejectionMethod.cpp AcceptanceRejectionMethod.h
        BlackScholes.h BlackScholes.cpp
        BoxMuller.cpp BoxMuller.h
        InverseTransform.cpp InverseTransform.h
        InverseTransformMethod.cpp InverseTransformMethod.h
        LinearCongruentialGenerator.cpp LinearCongruentialGenerator.h
        MonteCarloNonPathDependant.cpp MonteCarloNonPathDependant.h
        MonteCarloPathDependant.cpp MonteCarloPathDependant.h
        NormalGenerator.cpp NormalGenerator.h
        FiniteDifferencePricer.cpp FiniteDifferencePricer.h
        LinearCongruentialGenerator.h LinearCongruentialGenerator.cpp
 */
void Problem1(){
	// Finite difference method for down and out call price
	double S = 42;
	double K = 40;
	double B = 35;
	double T = 7.0/12.0;
	double r = 0.05;
	double q = 0.03;
	double vol = 0.30;
	double alpha = 0.40;
	unsigned int M=4;

	// I/O setting
	ofstream myfile;
	string filename = "HW10_PROBLEM.csv";
	myfile.open(filename.c_str());


	// BS pricer
	shared_ptr<BSpricer> bsp(new BSpricer());
	auto bsprice_1 = bsp->GetPrice(S, K, T, vol, q, r, 'c');
	auto bsprice_2 = bsp->GetPrice(B*B/S, K, T, vol, q, r, 'c');
	auto tmp_a = (r-q)/vol/vol - 0.5;
	auto benchmark = bsprice_1 - pow((B/S),2*tmp_a) * bsprice_2;
	cout << "Closed-form Down–and–Out Call price: " << benchmark <<endl;

	myfile << "Pricing Down–and–Out European Barrier Options using Finite Differences\n\n";

	{
		// first line tables
		// first table on the left
		myfile << "Domain discretization with alpha_temp = 0.4 \n\n";
		myfile << "M,alpha,x_left,x_right,N,dx,dt \n";
		alpha = 0.4;
		M = 1;
		for (int i=1; i<=4; ++i){
			M *= 4;
			shared_ptr<FiniteDifferencePricer> fdp2(new FiniteDifferencePricer(S,K,B,vol,r,q,T,M,alpha));
			myfile << std::setprecision(precision)<< M << "," << fdp2->alpha << ","<<fdp2->x_left << "," << fdp2->x_right << "," << fdp2->N << ","<< fdp2->dx << "," << fdp2->dt << "\n";
		}

		myfile << "\n\n\n";

		// second table on the right
		myfile << "Domain discretization with alpha_temp = 4 \n\n";
		myfile << "M,alpha,x_left,x_right,N,dx,dt \n";
		alpha = 4;
		M =1;
		for (int i=1; i<=4; ++i){
			M *= 4;
			shared_ptr<FiniteDifferencePricer> fdp2(new FiniteDifferencePricer(S,K,B,vol,r,q,T,M,alpha));
			myfile << std::setprecision(precision)<< M << "," << fdp2->alpha << "," <<fdp2->x_left << "," << fdp2->x_right << "," << fdp2->N << "," << fdp2->dx << "," << fdp2->dt << "\n";
		}

		myfile << "\n\n\n";
	}

	vector<vector<double>> grid_forward;
	{
		// 2nd line table
		// left table
		myfile << "Forward Euler with alpha_temp = 0.4\n\n";
		myfile << "\n M,u_value,Option Value, Pointwise Error, Delta_central, Gamma_central, Theta_forward\n";
		alpha = 0.4;
		M = 1;
		vector<vector<double>> table2;
		table2.resize(4,vector<double>(7));
		for (int i=1; i<=4; ++i){
			M *= 4;
			shared_ptr<FiniteDifferencePricer> fdp2(new FiniteDifferencePricer(S,K,B,vol,r,q,T,M,alpha));
			fdp2->run('e',"forward");
			auto g2 = fdp2->grid;
			if (M==4) grid_forward = g2;
			table2[i-1][0] = M;
			table2[i-1][1] = g2[M][fdp2->N_left];
			table2[i-1][2] = fdp2->Price;
			table2[i-1][3] = fabs(benchmark-fdp2->Price);
			table2[i-1][4] = fdp2->Delta;
			table2[i-1][5] = fdp2->Gamma;
			table2[i-1][6] = fdp2->Theta;
		}
		for (auto row: table2){
			for (int i=0; i<row.size(); ++i){
				myfile << std::setprecision(precision)<< row[i] << ", ";
			}
			myfile << "\n";
		}
		myfile << "\n\n\n";

	}

	vector<vector<double>> grid_backward;
	{
		// 3rd line table, 1st table
		myfile << "Backward Euler with LU and alpha_temp = 0.4\n\n";
		myfile << "\n M,u_value,Option Value, Pointwise Error, Delta_central, Gamma_central, Theta_forward\n";
		alpha = 0.4;
		M = 1;
		vector<vector<double>> table2;
		table2.resize(4,vector<double>(7));
		for (int i=1; i<=4; ++i){
			M *= 4;
			shared_ptr<FiniteDifferencePricer> fdp2(new FiniteDifferencePricer(S,K,B,vol,r,q,T,M,alpha));
			fdp2->run('e',"backward");
			auto g2 = fdp2->grid;
			if (M==4) grid_backward = g2;
			table2[i-1][0] = M;
			table2[i-1][1] = g2[M][fdp2->N_left];
			table2[i-1][2] = fdp2->Price;
			table2[i-1][3] = fabs(benchmark-fdp2->Price);
			table2[i-1][4] = fdp2->Delta;
			table2[i-1][5] = fdp2->Gamma;
			table2[i-1][6] = fdp2->Theta;
		}
		for (auto row: table2){
			for (int i=0; i<row.size(); ++i){
				myfile<< std::setprecision(precision) << row[i] << ", ";
			}
			myfile << "\n";
		}
		myfile << "\n\n\n";

	}

	{
		// 3rd line table, 2nd table
		myfile << "Backward Euler with LU and alpha_temp = 4\n\n";
		myfile << "\n M,u_value,Option Value, Pointwise Error, Delta_central, Gamma_central, Theta_forward\n";
		alpha = 4;
		M = 1;
		vector<vector<double>> table2;
		table2.resize(4,vector<double>(7));
		for (int i=1; i<=4; ++i){
			M *= 4;
			shared_ptr<FiniteDifferencePricer> fdp2(new FiniteDifferencePricer(S,K,B,vol,r,q,T,M,alpha));
			fdp2->run('e',"backward");
			auto g2 = fdp2->grid;
			table2[i-1][0] = M;
			table2[i-1][1] = g2[M][fdp2->N_left];
			table2[i-1][2] = fdp2->Price;
			table2[i-1][3] = fabs(benchmark-fdp2->Price);
			table2[i-1][4] = fdp2->Delta;
			table2[i-1][5] = fdp2->Gamma;
			table2[i-1][6] = fdp2->Theta;
		}
		for (auto row: table2){
			for (int i=0; i<row.size(); ++i){
				myfile << std::setprecision(precision)<< row[i] << ", ";
			}
			myfile << "\n";
		}
		myfile << "\n\n\n";
	}

	{
		// 4th line table, 1st table
		myfile << "Crank Nicolson with SOR and alpha_temp = 0.4\n\n";
		myfile << "\n M,u_value,Option Value, Pointwise Error, Delta_central, Gamma_central, Theta_forward\n";
		alpha = 0.4;
		M = 1;
		vector<vector<double>> table2;
		table2.resize(4,vector<double>(7));
		for (int i=1; i<=4; ++i){
			M *= 4;
			shared_ptr<FiniteDifferencePricer> fdp2(new FiniteDifferencePricer(S,K,B,vol,r,q,T,M,alpha));
			fdp2->run('e',"cranknicolsonUsingSOR");
			auto g2 = fdp2->grid;
			table2[i-1][0] = M;
			table2[i-1][1] = g2[M][fdp2->N_left];
			table2[i-1][2] = fdp2->Price;
			table2[i-1][3] = fabs(benchmark-fdp2->Price);
			table2[i-1][4] = fdp2->Delta;
			table2[i-1][5] = fdp2->Gamma;
			table2[i-1][6] = fdp2->Theta;
		}
		for (auto row: table2){
			for (int i=0; i<row.size(); ++i){
				myfile<< std::setprecision(precision) << row[i] << ", ";
			}
			myfile << "\n";
		}
		myfile << "\n\n\n";
	}

	{
		// 4th line table, 1st table
		myfile << "Crank Nicolson with SOR and alpha_temp = 4\n\n";
		myfile << "\n M,u_value,Option Value, Pointwise Error, Delta_central, Gamma_central, Theta_forward\n";
		alpha = 4;
		M = 1;
		vector<vector<double>> table2;
		table2.resize(4,vector<double>(7));
		for (int i=1; i<=4; ++i){
			M *= 4;
			shared_ptr<FiniteDifferencePricer> fdp2(new FiniteDifferencePricer(S,K,B,vol,r,q,T,M,alpha));
			fdp2->run('e',"cranknicolsonUsingSOR");
			auto g2 = fdp2->grid;
			table2[i-1][0] = M;
			table2[i-1][1] = g2[M][fdp2->N_left];
			table2[i-1][2] = fdp2->Price;
			table2[i-1][3] = fabs(benchmark-fdp2->Price);
			table2[i-1][4] = fdp2->Delta;
			table2[i-1][5] = fdp2->Gamma;
			table2[i-1][6] = fdp2->Theta;
		}
		for (auto row: table2){
			for (int i=0; i<row.size(); ++i){
				myfile << std::setprecision(precision)<< row[i] << ", ";
			}
			myfile << "\n";
		}
		myfile << "\n\n\n";
	}

	{
		// 5th line table
		M = 4;
		myfile << "Forward Euler with alpha_temp = 0.4\n\n";
		myfile << "m,x_left,x_compute,x_1,x_2,x_3,x_right\n";
		for (unsigned int i=0; i<=M; ++i){
			myfile << i;
			for (auto row: grid_forward[i]){
				myfile << "," << row;
			}
			myfile << "\n";
		}
	}

	{
		// 6th line table
		M = 4;
		myfile << "Backward Euler with LU and alpha_temp = 0.4\n\n";
		myfile << "m,x_left,x_compute,x_1,x_2,x_3,x_right\n";
		for (unsigned int i=0; i<=M; ++i){
			myfile<< std::setprecision(precision) << i;
			for (auto row: grid_backward[i]){
				myfile << std::setprecision(precision)<< "," << row;
			}
			myfile << "\n";
		}
	}

	// close the output file
	myfile.close();
}

void Problem1() {
	// Binomial tree methods for derivate valuation and hedging parameter computation
	ofstream fs;
	string filename = "problem1_data.csv";
	// create and open the .csv file
	fs.open(filename.c_str());
	// generate date and write data to the file
	BinomialPricer *bp = new BinomialPricer;
	double T = 10. / 12.;
	double S = 50;
	double K = 48;
	double r = 0.02;
	double vol = 0.30;
	double q = 0.01;
	double B = 45;

	double HW10;

	// write the file headers
	fs << "European\n\n";
	// HW10 values
	fs << "HW10:V," <<endl;
	fs << 1<<endl;
	fs << "N" << "," << "BT for Barrier" << "," << "ABT" << "," << "BBS" << "," << "BBSR" << std::endl;

	for (int N = 10; N <= 1000; N++) {
		// DOwn and out call euro
		fs << std::setprecision(precision) << N << "," << bp->DOCBinomail(N, S, K, T, vol, r, q, 1, 0, B) << endl;

		/*"," << bp->AvgBinomialTree(N,S,K,T,vol,r,q,1,0)<< ","
<< bp->BinomialBlackScholes(N,S,K,T,vol,r,q,1,0)<< "," << bp->BinomialBlackScholesWithRE(N,S,K,T,vol,r,q,1,0)<<endl;*/
		// compute prices of Euro put using different pricers
		/*
		fs << std::setprecision(precision) << N << "," << bp->BinomialTree(N, S, K, T, vol, r, q, -1, 0) << ","
		   << bp->AvgBinomialTree(N, S, K, T, vol, r, q, -1, 0) << ","
		   << bp->BinomialBlackScholes(N, S, K, T, vol, r, q, -1, 0) << ","
		   << bp->BinomialBlackScholesWithRE(N, S, K, T, vol, r, q, -1, 0) << endl;
		   */
	}
	/*
	fs << "\n\n\n\n\n";
	fs << "American\n\n";
	fs << "N" <<  "," << "BT" << "," << "ABT"<< "," << "BBS" << ","<< "BBSR" << std::endl;
	for (int N = 10; N <= 100; N++)
	{
		// compute prices of Amer put using different pricers
		fs << std::setprecision(precision) << N << "," << bp->BinomialTree(N,S,K,T,vol,r,q,-1,1) << "," << bp->AvgBinomialTree(N,S,K,T,vol,r,q,-1,1)<< ","
		   << bp->BinomialBlackScholes(N,S,K,T,vol,r,q,-1,1)<< "," << bp->BinomialBlackScholesWithRE(N,S,K,T,vol,r,q,-1,1)<<endl;
	}
*/
	// close the output file
	fs.close();

	delete bp;
}
/*
// implied vol using binomial tree
void BinomialTree(){
	// Binomial tree methods for derivate valuation and hedging parameter computation
	ofstream fs;
	string filename = "impvol.csv";
	// create and open the .csv file
	fs.open(filename.c_str());
	// generate date and write data to the file
	shared_ptr<BinomialPricer> bp(new BinomialPricer);
	unsigned int N = 2500;
	double T = 9./12.;
	double S = 42;
	double K = 45;
	double r = 0.04;
	double q = 0.02;
	double marketPrice = 4.09;

	// using binomial tree with variance reduction
	auto f = [&](double v){return bp->BinomialTree(N, S, K, T, v, r, q, -1, 2) - marketPrice;};

	//secant method
	double old_vol = 0.1, vol = 0.5;
	unsigned int counter = 0;
	while ((fabs(vol-old_vol) >0.0001) && (counter<10000)){
		cout << "Iterations: " << counter++ << ", old_vol: " << old_vol << ", vol: " << vol<<endl;
		double tmp = vol;
		vol = vol - f(vol) * (vol- old_vol) / (f(vol) - f(old_vol));
		old_vol = tmp;
	}

	cout << "Implied Volatility: " << vol <<endl;

	// write the file headers
	fs << "Implied volatility of an American put option\n\n";
	fs << "implied volatility," << vol << "\n";
	fs << "Iterations," << counter << "\n";

	fs << "\n\n";
	fs << "Remarks:, Binomial tree for American option with variance reduction is used here.";
	// close the output file
	fs.close();
}
*/
/*
void TrinomialTree(){
	// trinomial tree method for european options
	ofstream fs;
	string filename = "Trinomial.csv";
	// create and open the .csv file
	fs.open(filename.c_str());
	// generate date and write data to the file
	TrinomialPricer* bp = new TrinomialPricer;


	double T = 1;
	double S = 41;
	double K = 39;
	double r = 0.03;
	double vol = 0.25;
	double q = 0.005;
	auto tmp = bp->BSvalue(S,K,T,vol,r,q,-1);   // Black schole value for reference
	//cout << "BS value: " << tmp << endl;
	// BS values
	fs << "V_BS," << bp->BSprice << ",,Delta_BS," << bp->BSdelta << ",,Gamma_BS," << bp->BSgamma << ",,Theta_BS," << bp->BStheta << "\n" <<endl;

	fs << "Trinomial Tree\n" << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
	for (int N = 10; N<=1280; N*=2){
		double VN = bp->TrinomialTree(N,S,K,T,vol,r,q,-1,0);    // trinomial tree value of euro put
		// variable call = 1 means call, -1 means put
		// variable amer = 0 means it is euro put, 1 means amer put, 2 means amer put with variance reduction
		fs << N << ","
		   << VN << ","
		   << fabs(VN-bp->BSprice) <<","
		   << N*fabs(VN-bp->BSprice) <<","
		   << N*N*fabs(VN-bp->BSprice) <<","
		   << bp->delta << ","
		   << fabs(bp->delta - bp->BSdelta) << ","
		   << bp->gamma << ","
		   << fabs(bp->gamma - bp->BSgamma) << ","
		   << bp->theta << ","
		   << fabs(bp->theta - bp->BStheta) <<endl;
	}
	fs << "\n\n";
	fs << "Trinomial Black Scholes\n" <<
	   "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,"<<
	   "|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
	for (int N = 10; N<=1280; N*=2){
		double VN = bp->TrinomialBlackScholes(N,S,K,T,vol,r,q,-1,0);
		fs << N << ","
		   << VN << ","
		   << fabs(VN-bp->BSprice) <<","
		   << N*fabs(VN-bp->BSprice) <<","
		   << N*N*fabs(VN-bp->BSprice) <<","
		   << bp->delta << ","
		   << fabs(bp->delta - bp->BSdelta) << ","
		   << bp->gamma << ","
		   << fabs(bp->gamma - bp->BSgamma) << ","
		   << bp->theta << ","
		   << fabs(bp->theta - bp->BStheta) <<endl;
	}

	fs << "\n\n";
	fs << "Trinomial Black-Scholes with Richardson Extrapolation\n" <<
	   "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,"
	   << "Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
	for (int N = 10; N<=1280; N*=2){
		double VN = bp->TrinomialBlackScholesWithRE(N,S,K,T,vol,r,q,-1,0);
		fs << N << ","
		   << VN << ","
		   << fabs(VN-bp->BSprice) <<","
		   << N*fabs(VN-bp->BSprice) <<","
		   << N*N*fabs(VN-bp->BSprice) <<","
		   << bp->delta << ","
		   << fabs(bp->delta - bp->BSdelta) << ","
		   << bp->gamma << ","
		   << fabs(bp->gamma - bp->BSgamma) << ","
		   << bp->theta << ","
		   << fabs(bp->theta - bp->BStheta) << "\n";

	}

	fs << "\n\n";
	fs << "Rank the methods in terms of convergence speed (from fastest to slowest)\n";
	fs << "answers....\n";
	fs << "Comment:\n";
	fs << "answers....\n";


	// Double Barrier Option {} 去掉
	{
		double L = 40;
		double U =60;

		//int N =
		double T = 1;
		double S = 50;
		double K = 48;
		double r = 0.02;
		double vol = 0.3;
		double q = 0.005;

		TrinomialPricer* bp = new TrinomialPricer;
		fs << "Trinomial Tree\n" << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
		for (int N = 10; N<=1280; N*=2){ // k = 2^7
			double VN = bp->TrinomialTree(N,S,K,T,vol,r,q,1,0);    // trinomial tree value of euro put
			// variable call = 1 means call, -1 means put
			// variable amer = 0 means it is euro put, 1 means amer put, 2 means amer put with variance reduction
			fs << N << ","
			   << VN << ","
			   << fabs(VN-bp->BSprice) <<","
			   << N*fabs(VN-bp->BSprice) <<","
			   << N*N*fabs(VN-bp->BSprice) <<","
			   << bp->delta << ","
			   << fabs(bp->delta - bp->BSdelta) << ","
			   << bp->gamma << ","
			   << fabs(bp->gamma - bp->BSgamma) << ","
			   << bp->theta << ","
			   << fabs(bp->theta - bp->BStheta) <<endl;
		}
	}


	// close the output file
	fs.close();

	delete bp;
}
 */

int main() {
	//BinomialTree();
	//TrinomialTree();
	Problem1();

	return 0;
}
