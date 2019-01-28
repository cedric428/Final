// MTH9821 HW 2
// dhan 20170907
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;
#include "TrinomialPricer.h"
#include "BinomialPricer.h"
#include <iomanip>

string LOCALPATH = "/home/dhan/Downloads/";

void Problem1(){
    // trinomial tree method for european options
    ofstream fs;
    string filename = LOCALPATH + "problem1_data.csv";
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
    auto tmp = bp->BSvalue(S,K,T,vol,r,q,-1);
    //cout << "BS value: " << tmp << endl;
    // BS values
    fs << "V_BS," << bp->BSprice << ",,Delta_BS," << bp->BSdelta << ",,Gamma_BS," << bp->BSgamma << ",,Theta_BS," << bp->BStheta << "\n" <<endl;
    fs << "Trinomial Tree\n" << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
    for (int N = 10; N<=1280; N*=2){
        double VN = bp->TrinomialTree(N,S,K,T,vol,r,q,-1,0);
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
    // close the output file
    fs.close();

    delete bp;
}


void Problem2(){
    // trinomial tree method for American options
    ofstream fs;
    string filename = LOCALPATH + "problem2_data.csv";
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

    // Exact value using avg binomial tree
    // bp->BSvalue(S,K,T,vol,r,q,-1);
    BinomialPricer* abt = new BinomialPricer;
    auto priceWithoutVR = abt ->AvgBinomialTree(10000,S,K,T,vol,r,q,-1,1);
    auto deltaWithoutVR = abt ->delta;
    auto gammaWithoutVR = abt ->gamma;
    auto thetaWithoutVR = abt ->theta;
    // BS values
    fs << "V_BS," << priceWithoutVR << ",,Delta_BS," << deltaWithoutVR << ",,Gamma_BS," << gammaWithoutVR << ",,Theta_BS," << thetaWithoutVR << "\n" <<endl;
    delete abt;

    fs << "Trinomial Tree,,,,,,,,,,,, Variance Reduction\n"
    << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS|,,"
    << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
    for (int N = 10; N<=1280; N*=2){
        double VN = bp->TrinomialTree(N,S,K,T,vol,r,q,-1,1); // amer put without var reduction
        fs << N << ","
        << VN << ","
        << fabs(VN - priceWithoutVR) <<","
        << N*fabs(VN - priceWithoutVR) <<","
        << N*N*fabs(VN - priceWithoutVR) <<","
        << bp->delta << ","
        << fabs(bp->delta - deltaWithoutVR) << ","
        << bp->gamma << ","
        << fabs(bp->gamma - gammaWithoutVR) << ","
        << bp->theta << ","
        << fabs(bp->theta - thetaWithoutVR) << ",,";

        VN = bp->TrinomialTree(N,S,K,T,vol,r,q,-1,2); // amer put with var reduction
        fs << N << ","
        << VN << ","
        << fabs(VN-priceWithoutVR) <<","
        << N*fabs(VN-priceWithoutVR) <<","
        << N*N*fabs(VN-priceWithoutVR) <<","
        << bp->delta << ","
        << fabs(bp->delta - deltaWithoutVR) << ","
        << bp->gamma << ","
        << fabs(bp->gamma - gammaWithoutVR) << ","
        << bp->theta << ","
        << fabs(bp->theta - thetaWithoutVR) << endl;
    }
    fs << "\n\n";
    fs << "Trinomial Black Scholes,,,,,,,,,,,, Variance Reduction\n"
    << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS|,,"
    << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
    for (int N = 10; N<=1280; N*=2){
        double VN = bp->TrinomialBlackScholes(N,S,K,T,vol,r,q,-1,1);
        fs << N << ","
        << VN << ","
        << fabs(VN-priceWithoutVR) <<","
        << N*fabs(VN-priceWithoutVR) <<","
        << N*N*fabs(VN-priceWithoutVR) <<","
        << bp->delta << ","
        << fabs(bp->delta - deltaWithoutVR) << ","
        << bp->gamma << ","
        << fabs(bp->gamma - gammaWithoutVR) << ","
        << bp->theta << ","
        << fabs(bp->theta - thetaWithoutVR) << ",,";

        VN = bp->TrinomialBlackScholes(N,S,K,T,vol,r,q,-1,2); // amer put with var reduction
        fs << N << ","
        << VN << ","
        << fabs(VN-priceWithoutVR) <<","
        << N*fabs(VN-priceWithoutVR) <<","
        << N*N*fabs(VN-priceWithoutVR) <<","
        << bp->delta << ","
        << fabs(bp->delta - deltaWithoutVR) << ","
        << bp->gamma << ","
        << fabs(bp->gamma - gammaWithoutVR) << ","
        << bp->theta << ","
        << fabs(bp->theta - thetaWithoutVR) << endl;
    }

    fs << "\n\n";
    fs << "Trinomial Black-Scholes with Richardson Extrapolation,,,,,,,,,,,, Variance Reduction\n"
    << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS|,,"
    << "N,V(N),|V(N)-V_BS|,N*|V(N) - V_BS|,N^2*|V(N) - V_BS|,Delta_approx,|Delta_approx-Delta_BS|,Gamma_approx,|Gamma_approx-Gamma_BS|,Theta_approx,|Theta_approx-Theta_BS| \n";
    for (int N = 10; N<=1280; N*=2){
        double VN = bp->TrinomialBlackScholesWithRE(N,S,K,T,vol,r,q,-1,1);
        fs << N << ","
        << VN << ","
        << fabs(VN-priceWithoutVR) <<","
        << N*fabs(VN-priceWithoutVR) <<","
        << N*N*fabs(VN-priceWithoutVR) <<","
        << bp->delta << ","
        << fabs(bp->delta - deltaWithoutVR) << ","
        << bp->gamma << ","
        << fabs(bp->gamma - gammaWithoutVR) << ","
        << bp->theta << ","
        << fabs(bp->theta - thetaWithoutVR) << ",,";

        VN = bp->TrinomialBlackScholesWithRE(N,S,K,T,vol,r,q,-1,2); // amer put with var reduction
        fs << N << ","
        << VN << ","
        << fabs(VN-priceWithoutVR) <<","
        << N*fabs(VN-priceWithoutVR) <<","
        << N*N*fabs(VN-priceWithoutVR) <<","
        << bp->delta << ","
        << fabs(bp->delta - deltaWithoutVR) << ","
        << bp->gamma << ","
        << fabs(bp->gamma - gammaWithoutVR) << ","
        << bp->theta << ","
        << fabs(bp->theta - thetaWithoutVR) << "\n";

    }

    fs << "\n\n";
    fs << "Rank the methods in terms of convergence speed (from fastest to slowest)\n";
    fs << "answers....\n";
    fs << "Comment:\n";
    fs << "answers....\n";

    // close the output file
    fs.close();

    delete bp;
}


int main()
{

    cout << "Solving problems... Please wait...\n";
    Problem1();
    Problem2();
    cout << "Complete! Please check data.";

    return 0;
}
