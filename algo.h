#pragma once
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"

int feasible(double *lb, double *ub, double *mu, double **px, int n);
double eval_objective(double lambdaval, double *cov, double *mu, double *x, int n);
void optimize(int, double*, double*, double*, double*, double, double**, int num_iter = 100000);
