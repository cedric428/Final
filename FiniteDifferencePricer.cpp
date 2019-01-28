// FiniteDifferencePricer.cpp
// this file is for down and out call 

#include "FiniteDifferencePricer.h"
#include <iostream>
#include <cmath>
#include <memory>
#include "LinearSystemSolver.h"
using namespace std;

FiniteDifferencePricer::FiniteDifferencePricer(double S,double K, double B,
	double vol, double r, double q, double T, unsigned int M, double alpha_temp):
	S(S),K(K),B(B),vol(vol),r(r),q(q),T(T),M(M),alpha(alpha_temp),Price(-1),Delta(-1),Gamma(-1),Theta(-1){
	x_compute = log(S/K);
	tau_final = T*vol*vol/2.;

	dt = tau_final/M;
	dx = sqrt(dt/alpha_temp);

	x_left = log(B/K);
	N_left = floor((x_compute - x_left)/dx);
	dx = (x_compute-x_left)/N_left;
	alpha = dt/dx/dx;

	auto x_right_temp = log(S/K) + (r-q-vol*vol/2)*T + 3*vol*sqrt(T);
	N_right = ceil((x_right_temp-x_compute)/dx);

	N = N_left+N_right;
	x_right = x_compute + N_right * dx;

	a = (r-q)/vol/vol - 0.5;
	b = pow(((r-q)/vol/vol+ 0.5),2) + 2*q/vol/vol;

	grid.resize(M+1,vector<double>(N+1));

	// init time axis
	taxis.push_back(0.0);
	for (unsigned int i=1; i <= M; ++i){
		taxis.push_back(taxis[i-1]+dt);
	}

	// init x axis
	xaxis.push_back(x_left);
	for (unsigned int i=1; i <= N; ++i){
		xaxis.push_back(xaxis[i-1]+dx);
	}

	// early exercise premium
	early_ex_premium.resize(M+1,vector<double>(N+1));
	for (unsigned int m=0; m<M+1; ++m){
		for (unsigned int n=0; n<N; ++n){
			early_ex_premium[m][n] = K*exp(a*xaxis[n]+b*taxis[m])*max(1-exp(xaxis[n]),0.0);
		}
	}

	// system matrix
	// by defalut is backward
	if (N>=3){
		A.resize(N-1,vector<double>(N-1));
	}
}

FiniteDifferencePricer::~FiniteDifferencePricer(){
}

void FiniteDifferencePricer::SetBoundaryCondition(char euro){
	if (euro == 'e'){
		// european
		// fill the grid with boundary condition
		// u(x,0)
		for (unsigned int x=0; x<=N; ++x){
			grid[0][x] = K*exp(a*xaxis[x])*max(exp(xaxis[x])-1,0.0);
		}
		// u(x_left, tau)
		for (unsigned int t=0; t<=M; ++t){
			grid[t][0] = 0;
		}
		// u(x_right,tau)
		for (unsigned int t=0; t<=M; ++t){
			grid[t][N] = K*exp(a*x_right+b*taxis[t])*(exp(x_right-2*q*taxis[t]/vol/vol)-exp(-2*r*taxis[t]/vol/vol));
		}
	}
	else{
		cout << "American not support here" <<endl;
	}
}

vector<vector<double>> FiniteDifferencePricer::FillGridUsingForwardEuler(char euro){
	//cout << "alpha in fill grid: " <<alpha <<endl;
	if (euro=='e'){
		// use FDM to fill the grid
		for (unsigned int m=0; m<M; ++m){
			for (unsigned int n=1; n<N; ++n){
				grid[m+1][n] = alpha * grid[m][n-1] + (1-2*alpha)*grid[m][n] + alpha*grid[m][n+1];
			}
		}
	}
	else{
		// amer put
		// use FDM to fill the grid

		for (unsigned int m=0; m<M; ++m){
			for (unsigned int n=1; n<N; ++n){
				grid[m+1][n] = alpha * grid[m][n-1] + (1-2*alpha)*grid[m][n] + alpha*grid[m][n+1];
				grid[m+1][n] = max(grid[m+1][n], early_ex_premium[m+1][n]);
			}
		}
		
	}
	return grid;
}

vector<vector<double>> FiniteDifferencePricer::FillGridUsingBackwardEuler(){
	// set A
	A[0][0] = 1+2*alpha; A[0][1]= -alpha;
	for (unsigned int i=1;i<N-2; ++i){
		A[i][i] = 1+2*alpha;
		A[i][i-1] = -alpha;
		A[i][i+1] = -alpha;
	}
	A[N-2][N-2] = 1+2*alpha;
	A[N-2][N-3] = -alpha;
	// set boundary vector
	vector<double> boundaryvector(N-1,0);
	// define a solver
	shared_ptr<LinearSystemSolver> lss(new LinearSystemSolver());
	// compute decomposition first and reuse them
	lss->SetupTridiagSolverUsingLU(A);
	for (unsigned int i=1; i<M+1; ++i){
		boundaryvector[0] = alpha*grid[i][0]; 
		boundaryvector[N-2] = alpha*grid[i][N];
		auto y = std::vector<double>(grid[i-1].begin()+1, grid[i-1].begin()+N);
		auto res = lss->SolveTridiagExpress(addvector(y,boundaryvector,1.0));
		for (unsigned int j=1; j<N; ++j){
			grid[i][j] = res[j-1];
		}
	}
	return grid;
}

vector<vector<double>> FiniteDifferencePricer::FillGridUsingCrankNicolson(){
	// set A
	A[0][0] = 1+alpha; A[0][1]= -alpha/2;
	for (unsigned int i=1;i<N-2; ++i){
		A[i][i] = 1+alpha;
		A[i][i-1] = -alpha/2;
		A[i][i+1] = -alpha/2;
	}
	A[N-2][N-2] = 1+alpha;
	A[N-2][N-3] = -alpha/2;
	// set boundary vector
	vector<double> boundaryvector(N-1,0);
	// define a solver
	shared_ptr<LinearSystemSolver> lss(new LinearSystemSolver());
	// compute decomposition first and reuse them
	lss->SetupTridiagSolverUsingLU(A);
	for (unsigned int i=1; i<M+1; ++i){
		boundaryvector[0] = alpha*grid[i][0]/2; 
		boundaryvector[N-2] = alpha*grid[i][N]/2;
		//cout << "this row : \n";
		//if (M<10) { cout << i << ","; printvec(grid[i-1]);}
		auto y1 = std::vector<double>(grid[i-1].begin()+2, grid[i-1].end());
		auto y2 = std::vector<double>(grid[i-1].begin()+1, grid[i-1].end()-1);
		auto y3 = std::vector<double>(grid[i-1].begin(), grid[i-1].end()-2);
		auto bvec = addvector(boundaryvector, y1, alpha/2);
		//if (M<10) printvec(y1);
		//if (M<10) printvec(y2);
		//if (M<10) printvec(y3);
		//if (M<10) printvec(boundaryvector);


		bvec = addvector(bvec, y2, (1-alpha));
		bvec = addvector(bvec, y3, alpha/2);

		auto res = lss->SolveTridiagExpress(bvec);
		for (unsigned int j=1; j<N; ++j){
			grid[i][j] = res[j-1];
		}
	}
	return grid;
}

void FiniteDifferencePricer::ComputeResult(){
	Price = exp(-a*xaxis[N_left] - b*tau_final) * grid[M][N_left];

	auto S_minus1 = K*exp(xaxis[N_left-1]);
	auto S_0 = K*exp(xaxis[N_left]);
	auto S_1 = K*exp(xaxis[N_left+1]);
	auto V_minus1 = exp(-a*xaxis[N_left-1] - b*tau_final) * grid[M][N_left-1];
	auto V_0 = exp(-a*xaxis[N_left] - b*tau_final) * grid[M][N_left];
	auto V_1 = exp(-a*xaxis[N_left+1] - b*tau_final) * grid[M][N_left+1];
	auto V_approx_dt = exp(-a*xaxis[N_left] -b*(tau_final-dt)) * grid[M-1][N_left];


	Delta = (V_1 - V_minus1) / (S_1 - S_minus1);
	Gamma = ( (S_0-S_minus1 )*V_1 - (S_1-S_minus1)*V_0 +(S_1-S_0)*V_minus1 ) / ( (S_0-S_minus1)*(S_1-S_0)*((S_1-S_minus1)/2) );
	cout << "debug:" << Price << ", " << V_approx_dt << "," << dt <<endl;
	Theta = (Price - V_approx_dt)/(2*dt/vol/vol);

}


void FiniteDifferencePricer::run(char euro, string fdm){
	this->SetBoundaryCondition(euro);
	if (fdm=="backward" && euro=='e'){
		this->FillGridUsingBackwardEuler();
	}
	else if( fdm == "cranknicolson" && euro=='e'){
		this->FillGridUsingCrankNicolson();
	}
	else if (fdm == "cranknicolsonUsingSOR" && euro=='a'){
		this->FillGridUsingCrankNicolsonUsingSOR(euro);
	}
	else if (fdm == "cranknicolsonUsingSOR" && euro=='e'){
		this->FillGridUsingCrankNicolsonUsingSOR(euro);
	}
	else{
		this->FillGridUsingForwardEuler(euro);
	}
	this->ComputeResult();
}

vector<vector<double>> FiniteDifferencePricer::FindEarlyExerciseDomain() const{
	vector<vector<double>> domain;
	//printMatrix(early_ex_premium);
	//printMatrix(grid);
	domain.push_back(vector<double>{T, K});//at maturity,if S>K, should exercise 
	for (unsigned int m=0; m<=M; ++m){
		for (unsigned int n=0; n<N; ++n){
			if ((grid[m][n]==early_ex_premium[m][n]) && (grid[m][n+1]>early_ex_premium[m][n+1])){
				domain.push_back(vector<double>{T-2*m*dt/vol/vol, (K*exp(xaxis[n])+K*exp(xaxis[n+1]))/2 });
			}
		}
	}
	//printMatrix(domain);
	return domain;

}

void FiniteDifferencePricer::printMatrix(const vector<vector<double>>& g) const{
	cout << "[ \n";
    for (auto row: g){
        for (auto item:row){
            cout << item << " ";
        }
        cout << "\n";
    }
    cout << " ]" <<endl;
}

vector<double> FiniteDifferencePricer::addvector(const vector<double>& a, const vector<double> b, double factor){
	if (a.size()==b.size()){
		std::vector<double> v(a.size(),0.0);
		for (unsigned int i=0; i<a.size(); ++i){
			v[i] = a[i]+ factor * b[i];
		}
		return v;
	}
	return {};
}



void FiniteDifferencePricer::printvec(const vector<double>& v) const{
	cout << "[";
	for (auto& i: v){
		cout << i <<" ";
	}
	cout << "]" << endl;
}


vector<double> FiniteDifferencePricer::SOR ( const mat & A , const vec & b , const vec & x_0 ,
const double tolerance , const StoppingCriterion criterion , const double w){
	int n = b.size();
	vec res = x_0;
	int counter = 0;
	function<double (const vec&v1, const vec& v2)> f;
	switch(criterion)
	{
		case consecutive: {
			f = [&](const vec& v1, const vec& v2)-> double{return (v1-v2).norm();};
		}
		break;
		case residual: {
			f = [&](const vec& v1, const vec& v2)-> double{return (A*v1-b).norm();};
		}
		break;
	}
	auto res_old = b;
	while ( f(res,res_old) > tolerance && counter < 10000){
		counter++;
		res_old = res;
		for (int i=0; i<n; ++i){
			res(i) = 0;
			for (int k=0; k<n; ++k){
				if (k<i) res(i) += w * A(i,k)*res(k);
				if (k>i) res(i) += w * A(i,k)*res_old(k);
			}
			res(i) += (w-1) * A(i,i) * res_old(i);
			res(i) += -w* b(i);
			res(i) /= -(A(i,i));
		}
		//cout << res.transpose() <<endl; 
	}
	vector<double> res_vec;
	for (int i=0; i<n; ++i){
		res_vec.push_back(res(i));
	}
	return res_vec;
}

vector<vector<double>> FiniteDifferencePricer::FillGridUsingCrankNicolsonUsingSOR(char euro){
	// set boundary vector
	vector<double> boundaryvector(N-1,0);
	for (unsigned int i=1; i<M+1; ++i){
		boundaryvector[0] = alpha*grid[i][0]/2; 
		boundaryvector[N-2] = alpha*grid[i][N]/2;
		auto y1 = std::vector<double>(grid[i-1].begin()+2, grid[i-1].end());
		auto y2 = std::vector<double>(grid[i-1].begin()+1, grid[i-1].end()-1);
		auto y3 = std::vector<double>(grid[i-1].begin(), grid[i-1].end()-2);
		auto bvec = addvector(boundaryvector, y1, alpha/2);


		bvec = addvector(bvec, y2, (1-alpha));
		bvec = addvector(bvec, y3, alpha/2);

		// copy bvec into a vec object
		vec b_tmp(N-1);
		vec res_old(N-1);
		vec res(N-1);
		for (int k=0;k<N-1;++k){
			b_tmp(k) = bvec[k];
			res_old(k) = bvec[k];
			res(k) = early_ex_premium[i][k+1];
		}

		//cout << b_tmp <<endl;
		//cout << "|||||||\n" << res<<endl;
		double tolerance= 0.000001;
		int counter = 0;
		double w = 1.2;
		
		while ( (res-res_old).norm() > tolerance && counter < 10000){
			counter++;
			res_old = res;
			for (int j=0; j<N-1; ++j){
				double tmp1;
				if (j==0){
					tmp1 = res_old(j+1); // + grid[i][0];
				}
				else if (j==N-2){
					tmp1 = res(j-1); //+ grid[i][N];
				}
				else{
					tmp1 = res(j-1)+res_old(j+1);
				}
				res(j) = (1-w)*res_old(j)+(w*alpha/2/(1+alpha))*(tmp1) + w/(1+alpha)*b_tmp(j);
				//cout << res(j) << ",,,,,," << early_ex_premium[i][j+1]<<"\n";
				//cout << res(j)-early_ex_premium[i][j+1];
				if (euro == 'a') res(j) = max(res(j), early_ex_premium[i][j+1]);
			}
			// cout << res.transpose() <<endl; 
		}

		//cout << "SOR finished with "<< counter << " iterations"<<endl;
		for (unsigned int j=1; j<N; ++j){
			grid[i][j] = res(j-1);
		}
	}
	return grid;
}
