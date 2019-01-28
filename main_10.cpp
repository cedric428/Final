// MTH9821 HW 10 
//Pricing Down–and–Out European Barrier Options using Finite Differences
// dhan 20171115
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <memory>
#include <iomanip>  
#include "FiniteDifferencePricer.h"
#include "BSpricer.h"
#include "LinearSystemSolver.h"


using namespace std;
int precision = 18;

//string LOCALPATH = "/media/sf_Dropbox/1 Baruch Studies/MTH9821 Numerical Methods for Finance 2017 Fall Sem1/hw 10/";

string LOCALPATH = "";

void printGrid(vector<vector<double>> g){
    cout << "[ \n";
    for (auto row: g){
        for (auto item:row){
            cout << item << " ";
        }
        cout << "\n";
    }
    cout << " ]" <<endl;
}

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
    string filename = LOCALPATH + "HW10_PROBLEM.csv";
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



int main()
{
    Problem1();


    return 0;
}
