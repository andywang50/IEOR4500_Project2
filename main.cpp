#define _CRT_SECURE_NO_WARNINGS

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "algo.h"

int readit(char *nameoffile, int *addressofn, double **, double **, double **, double **, double *);
void print_vector(double* v, int n){
    for (int i=0; i < n; ++i){
        printf("%f", v[i]);
    }
}

int main(int argc, char **argv)
{
    int retcode = 0;
    int n;
    double *lb, *ub, *covariance, *mu, lambda;
    
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
    

    int num_iter = 100000;
    for (int k = 0; k < num_iter; ++k){
        if (k%10000 == 0){
            printf("%d, %f \n",k, eval_objective(lambda, covariance, mu, x, n));
        }
        
        double* gk = scalerProduct(2 * lambda, matrixTimesVector(covariance, x, n), n);
        gk = vectorSubstract(gk, mu, n);
        double *y = (double *)calloc(n, sizeof(double));
        
        double best_improvement = DBL_MAX;
        // improvement phase 1
        for (int i = 0; i < n; ++i){
            double threshold = gk[i];
            double *y_cand = (double *)calloc(n, sizeof(double));
            double y_cand_sum = 0.0;
            for (int j = 0; j < n; ++j){
                if (i==j){
                    continue;
                }
                if (gk[j] < threshold){
                    double delta =ub[j] - x[j];
                    y_cand[j] = delta;
                    y_cand_sum += delta;
                }
                else{
                    double delta = lb[j] - x[j];
                    y_cand[j] = delta;
                    y_cand_sum += delta;
                }
            }
            y_cand[i] = -y_cand_sum;
            if ((x[i] + y_cand[i] <= ub[i]) && (x[i] + y_cand[i] >= lb[i])){
                double improvement = dotProduct(gk, y_cand, n);
                if (improvement < best_improvement){
                    best_improvement = improvement;
                    y = y_cand;
                }
            }
        }
        //print_vector(y,n);
        if (fabs(best_improvement) < 1e-6) break;
        // improvement phase 2
        double* Sy = matrixTimesVector(covariance, y, n);
        double xSy = dotProduct(x, Sy, n);
        double s = dotProduct(mu ,y, n) - 2 * lambda * xSy;
        s = s / (2 * lambda);
        s = s / dotProduct(y, Sy, n);
        if (s < 0){
            s = 0;
        }
        if (s > 1){
            s = 1;
        }
        x = vectorAdd(x, scalerProduct(s, y, n), n);
 
        
        
    }
	return retcode;
//BACK:
//    return retcode;
}

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
