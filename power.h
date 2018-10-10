#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"

void run_power_once(double *matrix, int n, double *evalue, double **evector);

void run_power(double *matrix, int n, int num_spectra, double **new_matrix);

void run_power(double *matrix, int n, double tol, double **new_matrix);