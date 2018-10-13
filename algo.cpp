//
//  algo.cpp
//  Project2
//
//  Created by WangGuoan on 10/4/18.
//  Copyright © 2018 Guoan Wang. All rights reserved.
//

#include "algo.h"


int feasible(double *lb,double *ub, double *mu, double **px, int n){
    double lb_sum = 0.0;
    double ub_sum = 0.0;
    for (int i = 0; i < n; ++i){
        lb_sum += lb[i];
        ub_sum += ub[i];
    }
    
    if ((lb_sum > 1.0) || (ub_sum < 1.0)){
        return 1;
    }
    
    double *x = (double*) calloc(n, sizeof(double));
    if (x == NULL) {
        printf("no memory for x\n");
        return 1;
    }
    for (int i =0; i < n; ++i){
        x[i] = lb[i];
    }
    
    double sumx = lb_sum;
    
    for (int i = 0; i < n; ++i){
        if (sumx + (ub[i] - lb[i]) >= 1.0){
            x[i] = 1.0 - sumx + lb[i];
            break;
        }
        else{
            x[i] = ub[i];
            double delta = ub[i] - lb[i];
            sumx += delta;
        }
    }
    
    *px = x;
    return 0;
}

double eval_objective(double lambdaval, double *cov, double *mu, double *x, int n){
    double* tmp1 = matrixTimesVector(cov, x, n);
    double tmp2 = dotProduct(tmp1, x, n);
    double tmp3 = dotProduct(mu,x,n);
    return lambdaval * tmp2 - tmp3;
    
}
