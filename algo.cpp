//
//  algo.cpp
//  Project2
//
//  Created by WangGuoan on 10/4/18.
//  Copyright © 2018 Guoan Wang. All rights reserved.
//

#include "algo.h"

/*
	Checks whether QP is feasible. If yes, initiaites a feasible solution.
	Parameters:
		double *lb: lower bound array
		double *ub: upper bound array
		double *mu: mu array
		double **px: address of solution array
		int n: size of the QP problem
	Returns:
		retcode
*/
int feasible(double *lb, double *ub, double *mu, double **px, int n) {
	double lb_sum = 0.0;
	double ub_sum = 0.0;
	for (int i = 0; i < n; ++i) {
		lb_sum += lb[i];
		ub_sum += ub[i];
	}

	if ((lb_sum > 1.0) || (ub_sum < 1.0)) {
		return 1;
	}

	double *x = (double*)calloc(n, sizeof(double));
	if (x == NULL) {
		printf("no memory for x\n");
		return 1;
	}
	for (int i = 0; i < n; ++i) {
		x[i] = lb[i];
	}

	double sumx = lb_sum;

	for (int i = 0; i < n; ++i) {
		if (sumx + (ub[i] - lb[i]) >= 1.0) {
			x[i] = 1.0 - sumx + lb[i];
			break;
		}
		else {
			x[i] = ub[i];
			double delta = ub[i] - lb[i];
			sumx += delta;
		}
	}

	*px = x;
	return 0;
}

/*
	Evaluate the objective function
	Parameters:
		double lambdaval: value of lambda
		double *cov: covariance matrix array
		double *mu: mu array
		double *x: current solution array
		int n: size of the QP problem
	Returns:
		the value of the objective function
*/
double eval_objective(double lambdaval, double *cov, double *mu, double *x, int n) {
	double* tmp1 = matrixTimesVector(cov, x, n);
	double tmp2 = dotProduct(tmp1, x, n);
	double tmp3 = dotProduct(mu, x, n);

	free(tmp1);
	return lambdaval * tmp2 - tmp3;

}

/*
	algorithm2 as described in the pdf
	Parameters:
		int n: size of the Quadratice Programming problem
		double *lb lower bound array
		double *ub upper bound array
		double *mu mu array
		double *covariance covariance array
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
	printf("Result: %f \n", eval_objective(lambda, covariance, mu, *px, n));
}