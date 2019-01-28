// LinearSystemSolver.cpp

#include "LinearSystemSolver.h"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
using namespace std;

void LinearSystemSolver::DecomposeTridiagUsingChol(){
	unsigned int n = MTX.size();
	for (unsigned int i=0; i<n-1; ++i){
		U[i][i] = sqrt(MTX[i][i]);
		L[i][i] = U[i][i];
		U[i][i+1] = MTX[i][i+1]/U[i][i];
		L[i+1][i] = U[i][i+1];
		MTX[i+1][i+1] -=  U[i][i+1]*U[i][i+1];
	}
	U[n-1][n-1] = sqrt(MTX[n-1][n-1]);
	L[n-1][n-1] = U[n-1][n-1];
}

void LinearSystemSolver::DecomposeTridiagUsingLU(){
	unsigned int n = MTX.size();
	for (unsigned int i=0; i<n-1; ++i){
		L[i][i] = 1;
		L[i+1][i] = MTX[i+1][i]/MTX[i][i];
		U[i][i]= MTX[i][i];
		U[i][i+1] = MTX[i][i+1];
		MTX[i+1][i+1] -=  L[i+1][i]*U[i][i+1];
	}
	L[n-1][n-1]=1;
	U[n-1][n-1] = MTX[n-1][n-1];
}

vector<double> LinearSystemSolver::ForwardSubstitutionTridiag(const vector<vector<double>> H, const vector<double>& b){
	vector<double> res(b.size(),1);
	res[0] = b[0]/H[0][0];
	for (unsigned int i=1; i<b.size(); ++i){
		res[i] = (b[i] - H[i][i-1]*res[i-1])/H[i][i];
	}
	return res;
}

vector<double> LinearSystemSolver::BackwardSubstitutionTridiag(const vector<vector<double>> H, const vector<double>& b){
	unsigned int n = b.size();
	vector<double> res(n,0);
	res[n-1] = b[n-1]/H[n-1][n-1];
	for (int i=n-2; i>=0; --i){
		res[i] = (b[i]-H[i][i+1]*res[i+1])/H[i][i];
	}
	return res;
}

vector<double> LinearSystemSolver::SolveTridiagUsingCholesky(vector<vector<double>> A,const vector<double>& b){
	MTX = A;
	L.resize(A.size(),vector<double>(A.size()));
	U.resize(A.size(),vector<double>(A.size()));
	DecomposeTridiagUsingChol();
	auto y = ForwardSubstitutionTridiag(L,b);
	auto x = BackwardSubstitutionTridiag(U,y);
	//printmatrix(L);
	//printmatrix(U);
	//printvec(x);
	return x;
}

vector<double> LinearSystemSolver::SolveTridiagUsingLU(vector<vector<double>> A,const vector<double>& b){
	MTX = A;
	L.resize(A.size(),vector<double>(A.size()));
	U.resize(A.size(),vector<double>(A.size()));
	DecomposeTridiagUsingLU();
	auto y = ForwardSubstitutionTridiag(L,b);
	auto x = BackwardSubstitutionTridiag(U,y);
	//printmatrix(L);
	//printmatrix(U);
	//printvec(x);
	return x;
}

void LinearSystemSolver::SetupTridiagSolverUsingChol(vector<vector<double>> A){
	MTX = A;
	L.resize(A.size(),vector<double>(A.size()));
	U.resize(A.size(),vector<double>(A.size()));
	DecomposeTridiagUsingChol();
}

void LinearSystemSolver::SetupTridiagSolverUsingLU(vector<vector<double>> A){
	MTX = A;
	L.resize(A.size(),vector<double>(A.size()));
	U.resize(A.size(),vector<double>(A.size()));
	DecomposeTridiagUsingLU();
}


vector<double> LinearSystemSolver::SolveTridiagExpress(const vector<double>& b){
	auto y = ForwardSubstitutionTridiag(L,b);
	auto x = BackwardSubstitutionTridiag(U,y);
	return x;
}

void LinearSystemSolver::printvec(vector<double> v) const{
	cout << "[";
	for (auto& i: v){
		cout << i <<" ";
	}
	cout << "]" << endl;
}

void LinearSystemSolver::printmatrix(vector<vector<double>> m) const{
	cout << "[";
	for (auto& v: m){
		printvec(v);
	}
	cout << "]" << endl;
}
