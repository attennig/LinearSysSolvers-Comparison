#include "MatrixUtils.h"
#include <mkl.h>
#include <stdio.h>
#include <math.h>
#define MIN -100
#define MAX 100
int index(int i, int j, int m) {
	return i * m + j;
}

void generate_random_matrix(double* M, int n, int m) {
	for (int i = 0; i < m * n; i++)
		M[i] = rand_coeff();
	return;
}

void generate_random_convergent_matrix(double* M, int n, int m) {
	// Simmetric, diagonal elements > 0 and strictly diag dominant
	double min_value, rnd_value;
	for (int i = 0; i < n; i++) {
		for (int j = n-1; j >= i; j--) {
			rnd_value = rand_coeff();
			if (i == j) {
				min_value = 0;
				for (int k = 0; k < n; k++) {			
					if (k != i) min_value += fabs(M[index(i, k, m)]);	
				}		
				// To be strictily diag dominant Mij should be > min_value
				M[index(i, j, m)] = min_value + rnd_value;
			}
			else {
				// Symmetric 
				M[index(i, j, m)] = rnd_value;
				M[index(j, i, m)] = rnd_value;
			}
		}
	}
		
	return;
}


void print_matrix(double* M, int n, int m) {
	for (int i = 0; i < m * n; i++) {
		printf("%f ", M[i]);
		if (i % m == m - 1) printf("\n");
	}
	return;
}

void swap_rows(double* M, int m, int i, int k) {
	double app;
	for (int j = 0; j < m; j++) {
		app = M[index(i, j, m)];
		M[index(i, j, m)] = M[index(k, j, m)];
		M[index(k, j, m)] = app;
	}
	return;
}

double rand_coeff() {
	return MIN + (MAX - MIN) * (rand() / ((double)RAND_MAX));
}

double* transpose(double* M, int n, int m) {
	double* T = (double*)malloc(n * m * sizeof(double));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			T[index(j, i, n)] = M[index(i, j, m)];
		}
	}
	return T;
}