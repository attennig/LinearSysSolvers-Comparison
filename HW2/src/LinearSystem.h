#pragma once

class LinearSystem {
	public:
		int nRows;
		int nCols;
		double *A;
		double *b;
		double *x;

		// Constructors 
		LinearSystem(int n, int m);
		LinearSystem(int n, int m, double *M, double *kt);
		// Destructor
		~LinearSystem();
		// Setters
		void setCoefficientMatrix(double *M);
		void setKnownTermVector(double *kt);
		// Utils
		void print();
		void printSolution();
		// Solvers
		void GEPP();
		void GENP();
		void CHOLESKY();
		void MKL();
		double getError();
	private:
		void forwardSubstitution(double *M, double *y);
		void backwardSubstitution(double *M);

		
};