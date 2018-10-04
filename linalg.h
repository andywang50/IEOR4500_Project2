#pragma once
#include <stdio.h>
#include <stdlib.h>

double* vectorAdd(double*, double*, int);
double* scalerProduct(double scaler, double* v, int n);
double* vectorSubstract(double* v1, double* v2, int n);
double dotProduct(double* v1, double* v2, int n);
double* matrixTimesVector(double* mat, double* v, int n);