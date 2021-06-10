#include "LinearSystem.h"
#include "MatrixUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <mkl.h>

LinearSystem::LinearSystem(int n, int m) {
	nRows = n;
	nCols = m;
	A = (double*)malloc(n * m * sizeof(double));
	b = (double*)malloc(n * sizeof(double));
	x = (double*)malloc(n * sizeof(double));
}

LinearSystem::LinearSystem(int n, int m, double* M, double* kt) {
	nRows = n;
	nCols = m;
	A = M;
	b = kt;
	x = (double*)malloc(n * sizeof(double));
}

LinearSystem::~LinearSystem() {
}
void LinearSystem::setCoefficientMatrix(double* M) {
	A = M;
}
void LinearSystem::setKnownTermVector(double* kt) {
	b = kt;
}

void LinearSystem::print() {
	printf("______________\n");
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			//printf("%f ", M[m*i+j]);
			printf("%f ", A[index(i, j,nCols)]);
		}
		printf("| %f \n", b[i]);
	}
	printf("______________\n");
	return;

}

void LinearSystem::printSolution() {
	for (int i = 0; i < nRows; i++) {
		printf("x%d = %f; \n", i + 1, x[i]);
	}
	return;
}

void LinearSystem::GEPP() {
	// GEPP = Gaussian Elimination Partial Pivoting
	const double epsilon = 0.1;
	for (int i = 0; i < nRows; i++) {
		if (fabs(A[index(i, i, nCols)]) <= epsilon) {
			int max = i;
			for (int k = i + 1; k < nRows; k++) {
				if (fabs(A[index(k, i, nCols)]) > fabs(A[index(max, i, nCols)])) max = k;
			}
			swap_rows(A, nCols, i, max);
			swap_rows(b, 1, i, max);
		}
	}
	GENP();
	return;
}

void LinearSystem::GENP() {
	// GENP = Gaussian Elimination No Pivoting
	double* M = (double*)malloc(nRows * nCols * sizeof(double));
	memcpy(M, A, nRows * nCols * sizeof(double));
	double* original_b = (double*)malloc(nRows * sizeof(double));
	memcpy(original_b, b, nRows * sizeof(double));

	for (int k = 0; k < nRows - 1; k++) {
		for (int i = k + 1; i < nRows; i++) {
			double m = M[index(i, k, nCols)] / M[index(k, k, nCols)];
			for (int j = k; j < nCols; j++) {
				M[index(i, j, nCols)] = M[index(i, j, nCols)] - m * M[index(k, j, nCols)];
			}
			b[i] = b[i] - m * b[k];
		}
	}

	backwardSubstitution(M);

	memcpy(b, original_b, nRows * sizeof(double));

	free(M);
	free(original_b);
	return;
}

void LinearSystem::CHOLESKY() {

	assert(nRows == nCols);

	double* M = (double*)malloc(nRows * nCols * sizeof(double));
	memcpy(M, A, nRows * nCols * sizeof(double));
	double* L = (double*)malloc(nRows * nCols * sizeof(double));
	double* LT = (double*)malloc(nRows * nCols * sizeof(double));
	for (int i = 0; i < nRows * nCols; i++) L[i] = 0;
	for (int i = 0; i < nRows; i++) {
		L[index(i, i, nCols)] = sqrt(M[index(i, i, nCols)]);
		if (i < nRows - 1) {
			for (int j = i + 1; j < nRows; j++) {
				L[index(j, i, nCols)] = M[(index(j, i, nCols))] / L[index(i, i, nCols)];
				for (int k = i + 1; k < nCols; k++) {
					M[(index(j, k, nCols))] -= L[index(j, i, nCols)]*L[index(k, i, nCols)];
				}
			}
		}
	}

	double *y = (double*)malloc(nRows * sizeof(double));
	double *original_b = (double*)malloc(nRows * sizeof(double));
	memcpy(original_b, b, nRows * sizeof(double));
	forwardSubstitution(L, y);
	b = y;
	memcpy(b, y, nRows * sizeof(double));
	LT = transpose(L, nRows, nCols);
	backwardSubstitution(LT);
	memcpy(b, original_b, nRows * sizeof(double));
	free(M);
	free(L);
	free(LT);
	free(original_b);
}

void LinearSystem::MKL() {
	int *ipiv = (int*)malloc(nRows*sizeof(int));
	memcpy(x, b, nRows * sizeof(double));
	double* M = (double*)malloc(nRows * nCols * sizeof(double));
	memcpy(M, A, nRows * nCols * sizeof(double));
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, nRows, 1, M, nRows, ipiv, x, 1);
	free(M);
	free(ipiv);
}

void LinearSystem::forwardSubstitution(double *M, double *y) {
	// Forward substitution
	double sum;
	for (int i = 0; i < nRows; i++) {
		sum = 0;
		for (int j = 0; j < i; j++) {
			sum += M[index(i, j, nCols)] * y[j];
		}
		y[i] = (b[i] - sum) / M[index(i, i, nCols)];
	}
	return;
}

void LinearSystem::backwardSubstitution(double* M) {
	// Backward substitution
	double sum;
	for (int i = nRows - 1; i >= 0; i--) {
		sum = 0;
		for (int j = i + 1; j < nCols; j++) {
			sum += M[index(i, j, nCols)] * x[j];
		}
		x[i] = (b[i] - sum) / M[index(i, i, nCols)];
	}
	return;
}

double LinearSystem::getError() {
	// ||Ax-b||
	double error = 0;
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			error += A[index(i, j, nCols)] * x[j];
		}
		error -= b[i];
	}
	return error;
}