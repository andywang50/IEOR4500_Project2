# include "power.h"

void run_power_once(double *matrix, int n, double *evalue, double **evector) {
	double *v = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i) {
		v[i] = rand() / (double)RAND_MAX;
	}
	normalize(&v, n);
	int T = 10000;
	double tol = 1e-6;
	double oldnorm = 0.0;
	for (int t = 0; t < T; ++t) {
		v = matrixTimesVector(matrix, v, n);
		double normv = norm(v, n);
		normalize(&v, n);
		if (fabs(normv - oldnorm) < tol){
			break;
		}	
		oldnorm = normv;
	}

	*evalue = oldnorm;
	*evector = v;
}

void run_power(double *matrix, int n, int num_spectra, double **new_matrix) {
	double *m = (double*)calloc(n*n, sizeof(double));
	for (int i = 0; i < n*n; ++i) {
		m[i] = 0.0;
	}
	double *matrix_copy = (double*)calloc(n*n, sizeof(double));
	for (int i = 0; i < n*n; ++i) {
		matrix_copy[i] = matrix[i];
	}
	for (int count = 0; count < num_spectra; ++count) {
		double evalue;
		double* evector;
		run_power_once(matrix_copy, n, &evalue, &evector);
		double* tmp = crossProduct(evector, evector, n);
		tmp = scalerProduct(evalue, tmp, n*n);
		m = vectorAdd(m, tmp, n*n);
		matrix_copy = vectorSubstract(matrix_copy, tmp, n*n);

	}
	
	free(matrix_copy);
	*new_matrix = m;
}

void run_power(double *matrix, int n, double tol, double **new_matrix) {
	double *m = (double*)calloc(n*n, sizeof(double));
	for (int i = 0; i < n*n; ++i) {
		m[i] = 0.0;
	}
	double *matrix_copy = (double*)calloc(n*n, sizeof(double));
	for (int i = 0; i < n*n; ++i) {
		matrix_copy[i] = matrix[i];
	}

	double largest_evalue = 1.0 / (double)RAND_MAX;
	int run_next = 1;
	int count = 0;
	while (count < n) {
		double evalue;
		double* evector;
		run_power_once(matrix_copy, n, &evalue, &evector);

		if (largest_evalue == 1.0 / (double)RAND_MAX) {
			largest_evalue = evalue;
		}
		if (evalue / largest_evalue < tol) {
			break;
		}
		double* tmp = crossProduct(evector, evector, n);
		tmp = scalerProduct(evalue, tmp, n*n);
		m = vectorAdd(m, tmp, n*n);
		matrix_copy = vectorSubstract(matrix_copy, tmp, n*n);
		count++;
	}
		

	

	free(matrix_copy);
	*new_matrix = m;
}