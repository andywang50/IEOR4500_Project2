#include "linalg.h"

/*
	Sum of two vectors of the same size
	Parameters
		double *v1: first vector
		double *v2: second vector
		int n: size of the vectors
	Return
		sum of the two vectors as double*
*/
double* vectorAdd(double* v1, double* v2, int n) {
	double *result = (double *)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		result[i] = v1[i] + v2[i];
	}

	return result;
}

/*
	Scaler times a vector (or matrix)
	Parameters
		double scaler: scaler
		double *v: vector
		int n: size of the vectors
	Return
		the scaler product as double*
*/
double* scalerProduct(double scaler, double* v, int n) {
	double *result = (double *)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		result[i] = v[i] * scaler;
	}


	return result;
}

/*
	Difference between two vectors of the same size
	Parameters
		double *v1: first vector
		double *v2: second vector
		int n: size of the vectors
	Return
		difference between the two vectors as double*
*/
double* vectorSubstract(double* v1, double* v2, int n) {
	return vectorAdd(v1, scalerProduct(-1.0, v2, n), n);
}

/*
	Dot product of two vectors of the same size
	Parameters
		double *v1: first vector
		double *v2: second vector
		int n: size of the vectors
	Return
		Dot product of two vectors as double
*/
double dotProduct(double* v1, double* v2, int n) {
	double result = 0.0;
	for (int j = 0; j < n; j++) {
		result += (v1[j] * v2[j]);
	}
	return result;
}

/*
	norm of a vector
	Parameters
		double *v: vector
		int n: size of the vector
	Return
		norm of the vector as double
*/
double norm(double* v, int n) {
	double result = 0.0;
	result = dotProduct(v, v, n);
	result = sqrt(result);
	return result;
}

/*
	Cross product (outer product) of two vectors of the same size
	Parameters
		double *v1: first vector
		double *v2: second vector
		int n: size of the vectors
	Return
		Cross product (outer product) of the two vectors as double*
*/
double* crossProduct(double* v1, double* v2, int n) {
	double *result = (double *)calloc(n*n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result[i*n + j] = v1[i] * v2[j];
		}
	}
	return result;
}


/*
	Matrix times vector  (Assume mat is a square matrix of size n*n)
	Parameters
		double *mat: matrix
		double *v: vector
		int n: size of the vector
	Return
		Product as double*
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

/*
	normalize the vector in place
	Parameters
		double **pv: address of the vector
		int n: size of the vector
	Return
		void
 */
void normalize(double** pv, int n) {
	double normv = norm(*pv, n);
	for (int i = 0; i < n; ++i) {
		(*pv)[i] = (*pv)[i] / normv;
	}
}