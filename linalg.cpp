#include "linalg.h"

double* vectorAdd(double* v1, double* v2, int n) {
	double *result = (double *)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		result[i] = v1[i] + v2[i];
	}

	return result;
}



double* scalerProduct(double scaler, double* v, int n) {
	double *result = (double *)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		result[i] = v[i] * scaler;
	}


	return result;
}

double* vectorSubstract(double* v1, double* v2, int n) {
	return vectorAdd(v1, scalerProduct(-1.0, v2, n), n);
}

double dotProduct(double* v1, double* v2, int n) {
	double result = 0.0;
	for (int j = 0; j < n; j++) {
		result += (v1[j] * v2[j]);
	}
	return result;
}

/*
 Assume mat is a square matrix of size n*n
 */
double* matrixTimesVector(double* mat, double* v, int n) {
	double *result = (double *)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		result[i] = 0.0;
		for (int j = 0; j < n; ++j) {
			result[i] += (v[j] * mat[i*n + j]);
		}
	}

	return result;

}