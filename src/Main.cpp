#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "LinearSystem.h"
#include "MatrixUtils.h"

void Gaussian(const char* filename, int nTests, int minDim, int maxDim, int incDim) {
	clock_t start, end;
	double time_GENP, time_GEPP, time_MKL;
	double error_GENP, error_GEPP, error_MKL;
	FILE* out_file = fopen(filename, "w");
	fprintf(out_file, "n;timeGENP; errorGENP; timeGEPP; errorGEPP; timeMKL; errorMKL\n");
	fclose(out_file);
	for (int n = minDim; n <= maxDim; n = n + incDim) {
		double* A = (double*)malloc(n * n * sizeof(double));
		double* b = (double*)malloc(n * sizeof(double));
		generate_random_matrix(A, n, n);
		generate_random_matrix(b, n, 1);
		LinearSystem sys = LinearSystem(n, n, A, b);
		time_GENP = time_GEPP = time_MKL = error_GENP = error_GEPP = error_MKL = 0;
		for (int test = 0; test < nTests; test++) {
			printf("Dim: %d , %d / %d;\n", n, test+1, nTests);
			start = clock();
			sys.GENP();
			end = clock();
			time_GENP += ((double)(end - start) / CLOCKS_PER_SEC) / nTests;
			error_GENP += sys.getError() / nTests;

			start = clock();
			sys.GEPP();
			end = clock();
			time_GEPP += ((double)(end - start) / CLOCKS_PER_SEC) / nTests;
			error_GEPP += sys.getError() / nTests;
			
			start = clock();
			sys.MKL();
			end = clock();
			time_MKL += ((double)(end - start) / CLOCKS_PER_SEC) / nTests;
			error_MKL += sys.getError() / nTests;
		}
		out_file = fopen(filename, "a");
		fprintf(out_file, "%d;%.20f;%.20f;%.20f;%.20f;%.20f;%.20f\n", n, time_GENP, error_GENP, time_GEPP, error_GEPP, time_MKL, error_MKL);
		fclose(out_file);
	}
}


void Cholesky(const char* filename, int nTests, int minDim, int maxDim, int incDim) {
	clock_t start, end;
	double time_CHO, time_MKL, error_CHO, error_MKL;
	FILE* out_file = fopen(filename, "w");
	fprintf(out_file, "n;timeCHO;errorCHO;timeMKL;errorMKL\n");
	fclose(out_file);
	for (int n = minDim; n <= maxDim; n = n + incDim) {	
		double* A = (double*)malloc(n * n * sizeof(double));
		double* b = (double*)malloc(n * sizeof(double));
		generate_random_convergent_matrix(A, n, n);
		generate_random_matrix(b, n, 1);
		LinearSystem sys = LinearSystem(n, n, A, b);
		time_CHO = time_MKL = error_CHO = error_MKL = 0;
		for (int test = 0; test < nTests; test++) {
			printf("Dim: %d , %d / %d;\n", n, test+1, nTests);
			start = clock();
			sys.CHOLESKY();
			end = clock();
			time_CHO += ((double)(end - start) / CLOCKS_PER_SEC) / nTests;
			error_CHO += sys.getError() / nTests;

			start = clock();
			sys.MKL();
			end = clock();
			time_MKL += ((double)(end - start) / CLOCKS_PER_SEC) / nTests;
			error_MKL += sys.getError() / nTests;
		}
		out_file = fopen(filename, "a");
		fprintf(out_file, "%d;%.20f;%.20f;%.20f;%.20f\n", n, time_CHO, error_CHO, time_MKL, error_MKL);
		fclose(out_file);
	}
}
int main() {

	/* EX 1 */
	//printf("Ex 1\n");
	//Gaussian("out_Gaussian.csv", 5, 1000, 10000, 1000);

	/* EX 2 */
	printf("Ex 2\n");
	Cholesky("out_Cholesky.csv", 5, 1000, 10000, 1000);
	
	return 0;
}
