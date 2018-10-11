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
void print_vector(double*, int);
void print_square_matrix(double*, int);
void optimize(int, double*, double*, double*, double*, double, double**, int num_iter = 100000);


int main(int argc, char **argv)
{
    int retcode = 0;
    int n;
    double *lb, *ub, *covariance, *mu, lambda;
	srand((long)time(NULL));
    
    if (argc != 2){
        printf("usage: qp1 filename\n");  retcode = 1;
		return retcode;
    }
    
    retcode = readit(argv[1], &n, &lb, &ub, &mu, &covariance, &lambda);
    
    double *x;
    int is_feasible;
    is_feasible = feasible(lb, ub, mu, &x, n);
    if (is_feasible == 1){
        retcode = 2;
        return retcode;
    }
    
	optimize(n, lb, ub, mu, covariance, lambda, &x);
	print_vector(x, n);

	printf("Extra Credit\n\n");
	int num_spectra = 1;
	double tol = 0.01;
	double* new_cov;
	run_power(covariance, n, num_spectra, &new_cov);
	//run_power(covariance, n, tol, &new_cov);

	for (int i = 0; i < n; ++i) {
		new_cov[i*n + i] = covariance[i*n + i];
	}
	print_square_matrix(new_cov, n);
	optimize(n, lb, ub, mu, new_cov, lambda, &x);
	print_vector(x, n);

	free(lb);
	free(ub);
	free(covariance);
	free(mu);
	free(x);
	return retcode;
//BACK:
//    return retcode;
}

/*
	read data from datafile.
	Parameters:
		char* filename: filename
		int* address_of_n: address of n, where n is the size of the QP.
		double **plb£º address of lower bound array
		double **pub£º address of upper bound array
		double **pmu£º address of mu array
		double **pcovariance£º address of covariance array
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
    if (!datafile){
        printf("cannot open file %s\n", filename);
        readcode = 2;  goto BACK;
    }
    
    printf("reading data file %s\n", filename);
    
    fscanf(datafile, "%s", buffer);
    fscancode = fscanf(datafile, "%s", buffer);
    if (fscancode == EOF){
        printf("problem: premature file end at ...\n");
        readcode = 4; goto BACK;
    }
    
    n = *address_of_n = atoi(buffer);
    
    printf("n = %d\n", n);
    
    lb = (double *)calloc(n, sizeof(double));
    *plb = lb;
    ub = (double *)calloc(n, sizeof(double));
    *pub = ub;
    mu = (double *)calloc(n, sizeof(double));
    *pmu = mu;
    covariance = (double *)calloc(n*n, sizeof(double));
    *pcovariance = covariance;
    
    if (!lb || !ub || !mu || !covariance){
        printf("not enough memory for lb ub mu covariance\n"); readcode = 3; goto BACK;
    }
    
    fscanf(datafile, "%s", buffer);
    
    for (j = 0; j < n; j++){
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
    
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            fscanf(datafile, "%s", buffer);
            covariance[i*n + j] = atof(buffer);
        }
    }
    
    
    fscanf(datafile, "%s", buffer);
    if (strcmp(buffer, "END") != 0){
        printf("possible error in data file: 'END' missing\n");
    }
    
    
    fclose(datafile);
    
BACK:
    
    return readcode;
}

/*
	print a vector in console
	Parameters:
		double *v: array
		int n: size
	Returns:
		void
*/
void print_vector(double* v, int n) {
	for (int i = 0; i < n; ++i) {
		printf("%f, ", v[i]);
	}
	printf("\n");
}

/*
	print a matrix in console
	Parameters:
		double *matrix: the matrix array (assumes its a square matrix with size n)
		int n: size of the square matrix (n*n entries)
	Returns:
		void
*/
void print_square_matrix(double* mat, int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%f, ", mat[i*n + j]);
		}
		printf("\n");
	}
}

/*
	algorithm2 as described in the pdf
	Parameters:
		int n: size of the Quadratice Programming problem
		double *lb£º lower bound array
		double *ub£º upper bound array
		double *mu£º mu array
		double *covariance£º covariance array
		double lambda: lambda in the QP
		double **px: address of the solution array
		int num_iter: (optional) maximum iteration, default = 100000
	Returns:
		void
*/
void optimize(int n, double *lb, double *ub, double *mu, double *covariance, double lambda, double **px, int num_iter) {
	//print_vector(*px, n);
	for (int k = 0; k < num_iter; ++k) {
		if (k % 10000 == 0) {
			printf("%d, %f \n", k, eval_objective(lambda, covariance, mu, *px, n));
		}

		double* gk = scalerProduct(2 * lambda, matrixTimesVector(covariance, *px, n), n);
		gk = vectorSubstract(gk, mu, n);
		double *y = (double *)calloc(n, sizeof(double));

		double best_improvement = DBL_MAX;

		// improvement phase 1
		for (int i = 0; i < n; ++i) {
			double threshold = gk[i];
			double *y_cand = (double *)calloc(n, sizeof(double));
			double y_cand_sum = 0.0;
			for (int j = 0; j < n; ++j) {
				if (i == j) {
					continue;
				}
				if (gk[j] < threshold) {
					double delta = ub[j] - (*px)[j];
					y_cand[j] = delta;
					y_cand_sum += delta;
				}
				else {
					double delta = lb[j] - (*px)[j];
					y_cand[j] = delta;
					y_cand_sum += delta;
				}
			}
			y_cand[i] = -y_cand_sum;
			if (((*px)[i] + y_cand[i] <= ub[i]) && ((*px)[i] + y_cand[i] >= lb[i])) {
				double improvement = dotProduct(gk, y_cand, n);
				if (improvement < best_improvement) {
					best_improvement = improvement;
					y = y_cand;
				}
			}
		}

		free(gk);

		//print_vector(y,n);
		if (fabs(best_improvement) < 1e-6) break;
		// improvement phase 2
		double* Sy = matrixTimesVector(covariance, y, n);
		double xSy = dotProduct(*px, Sy, n);
		double s = dotProduct(mu, y, n) - 2 * lambda * xSy;
		s = s / (2 * lambda);
		s = s / dotProduct(y, Sy, n);
		if (s < 0) {
			s = 0;
		}
		if (s > 1) {
			s = 1;
		}
		*px = vectorAdd(*px, scalerProduct(s, y, n), n);

		free(y);
		free(Sy);
	}
}
