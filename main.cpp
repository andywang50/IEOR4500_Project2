#define _CRT_SECURE_NO_WARNINGS

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "algo.h"
#include "power.h"

int readit(char *nameoffile, int *addressofn, double **, double **, double **, double **, double *);

void cleanup(double **lb, double **ub, double **mu, double **covariance, double **x, double **new_cov);


int main(int argc, char **argv)
{
	int retcode = 0;
	int n;
	double *lb = NULL, *ub = NULL, *covariance = NULL, *mu = NULL, lambda;
	srand((long)time(NULL));

	if (argc != 2) {
		printf("usage: qp1 filename\n");  retcode = 1;
		return retcode;
	}

	retcode = readit(argv[1], &n, &lb, &ub, &mu, &covariance, &lambda);

	double *x = NULL;
	int is_feasible;
	is_feasible = feasible(lb, ub, mu, &x, n);
	if (is_feasible == 1) {
		retcode = 2;
		printf("Infeasible.");
		return retcode;
	}

	optimize(n, lb, ub, mu, covariance, lambda, &x);
	print_vector(x, n);

	printf("\n\nExtra Credit\n");
	int num_spectra = 1;
	printf("Number of spectra in PCA: %d\n", num_spectra);
	double tol = 0.01;
	double* new_cov = NULL;
	run_power(covariance, n, num_spectra, &new_cov);
	//run_power(covariance, n, tol, &new_cov);


	for (int i = 0; i < n; ++i) {
		new_cov[i*n + i] = covariance[i*n + i];
	}
	print_square_matrix(new_cov, n);
	optimize(n, lb, ub, mu, new_cov, lambda, &x);
	print_vector(x, n);

	cleanup(&lb, &ub, &mu, &covariance, &x, &new_cov);


	return retcode;
	//BACK:
	//    return retcode;
}

/*
	read data from datafile.
	Parameters:
		char* filename: filename
		int* address_of_n: address of n, where n is the size of the QP.
		double **plb address of lower bound array
		double **pub address of upper bound array
		double **pmu address of mu array
		double **pcovariance address of covariance array
		double *lambda: address of lambda
	Returns:
		retcode
*/
int readit(char *filename, int *address_of_n, double **plb, double **pub,
	double **pmu, double **pcovariance, double *lambda)
{
	int readcode = 0, fscancode;
	FILE *datafile = NULL;
	char buffer[100];
	int n, i, j;
	double *lb = NULL, *ub = NULL, *mu = NULL, *covariance = NULL;

	datafile = fopen(filename, "r");
	if (!datafile) {
		printf("cannot open file %s\n", filename);
		readcode = 2;  goto BACK;
	}

	printf("reading data file %s\n", filename);

	fscanf(datafile, "%s", buffer);
	fscancode = fscanf(datafile, "%s", buffer);
	if (fscancode == EOF) {
		printf("problem: premature file end at ...\n");
		readcode = 4; goto BACK;
	}

	n = *address_of_n = atoi(buffer);

	printf("n = %d\n", n);

	lb = (double *)calloc(n, sizeof(double));
	if (lb == NULL) {
		printf("not enough memory for lb\n"); readcode = 3; goto BACK;
	}
	*plb = lb;
	ub = (double *)calloc(n, sizeof(double));
	if (ub == NULL) {
		free(lb);
		printf("not enough memory for ub\n"); readcode = 3; goto BACK;
	}
	*pub = ub;
	mu = (double *)calloc(n, sizeof(double));
	if (mu == NULL) {
		free(lb);
		free(ub);
		printf("not enough memory for mu\n"); readcode = 3; goto BACK;
	}
	*pmu = mu;
	covariance = (double *)calloc(n*n, sizeof(double));
	if (covariance == NULL) {
		free(lb);
		free(ub);
		free(mu);
		printf("not enough memory for covariance\n"); readcode = 3; goto BACK;
	}
	*pcovariance = covariance;

	fscanf(datafile, "%s", buffer);

	for (j = 0; j < n; j++) {
		fscanf(datafile, "%s", buffer);
		fscanf(datafile, "%s", buffer);
		lb[j] = atof(buffer);
		fscanf(datafile, "%s", buffer);
		ub[j] = atof(buffer);
		fscanf(datafile, "%s", buffer);
		mu[j] = atof(buffer);
		printf("j = %d lb = %g ub = %g mu = %g\n", j, lb[j], ub[j], mu[j]);
	}


	fscanf(datafile, "%s", buffer);

	fscanf(datafile, "%s", buffer);
	*lambda = atof(buffer);
	printf("lambda = %g\n", *lambda);

	fscanf(datafile, "%s", buffer); /* reading 'covariance'*/

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fscanf(datafile, "%s", buffer);
			covariance[i*n + j] = atof(buffer);
		}
	}


	fscanf(datafile, "%s", buffer);
	if (strcmp(buffer, "END") != 0) {
		printf("possible error in data file: 'END' missing\n");
	}


	fclose(datafile);

BACK:

	return readcode;
}


void cleanup(double **lb, double **ub, double **mu, double **covariance, double **x, double **new_cov) {
	// free up dynamic allocated memory
	free(*lb);
	free(*ub);
	free(*mu);
	free(*covariance);
	free(*x);
	free(*new_cov);
}