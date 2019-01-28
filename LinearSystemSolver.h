// linear system solver
// dhan 20171020

#ifndef LSS_H
#define LSS_H

#include <vector>
#include <string>
using namespace std;


class LinearSystemSolver{
public:

	// can be used for general A and b
	vector<double> SolveTridiagUsingCholesky(vector<vector<double>> A,const vector<double>& b);
	vector<double> SolveTridiagUsingLU(vector<vector<double>> A,const vector<double>& b);

	// can be repeatedly used for different b given A
	void SetupTridiagSolverUsingChol(vector<vector<double>> A);
	void SetupTridiagSolverUsingLU(vector<vector<double>> A);
	vector<double> SolveTridiagExpress(const vector<double>& b);

private:
	vector<vector<double>> MTX;
	vector<vector<double>> L; // L (lower) matrix for LU/chol
	vector<vector<double>> U; // U (upper) matrix for LU/chol

	// decompose A which is tri-diagonal
	void DecomposeTridiagUsingChol();
	void DecomposeTridiagUsingLU();

	// backward and forward substitution
	vector<double> ForwardSubstitutionTridiag(vector<vector<double>> H, const vector<double>& b);
	vector<double> BackwardSubstitutionTridiag(vector<vector<double>> H, const vector<double>& b);

	void printmatrix(vector<vector<double>> m) const;
	void printvec(vector<double> v) const;
};


#endif