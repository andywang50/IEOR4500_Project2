# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 22:27:34 2018

@author: andy
"""


import numpy as np
import sys
from cvxopt import matrix, solvers 

from func import readdata, feasible, eval_objective
from power import runpower


def optimize(alldata, num_iter, covariance_type = 'original'):
	"""
	Optimization algorithm, as described in algo2.pdf
	Parameters
	---------
		alldata: a dictionary of the following form:
			'n': size of the QP problem
			'lower': lower bound vector
			'upper': upper bound vector
			'x': solution vector
			'lambda': value of lambda
			'covariance': covariance matrix
		num_iter: int, maximum number of iterations
		covariance_type: string
			'original': original version (question1)
			'modified': extra credit
			otherwise will raise exception
	Returns
	-------
		No returns
	"""
	assert covariance_type in ['original', 'modified']
	# our algo
	x = alldata['x']
	lambdaval = alldata['lambda']
	n = alldata['n']
	mu = alldata['mu']
	if covariance_type == 'original':
		cov = alldata['covariance']
	else:
		cov = alldata['modified_covariance']
	lb = alldata['lower']
	ub = alldata['upper']
	
	for k in range(0,num_iter):
		## improvement phase 1
		if k%10000 == 0:
			print(k, eval_objective(lambdaval, cov, mu, x))
		gk = 2 * lambdaval * np.matmul(cov, x) - mu
		y = np.zeros(n)
		best_improvement = np.inf
		argsort = gk.argsort();
		#ranks = np.empty_like(argsort)
		#ranks[argsort] = np.arange(len(gk))	
		for id_m in range(0,n):
			m = argsort[id_m]
			y_cand = np.zeros(n)
			for id_j in range(0,id_m):
				j = argsort[id_j]
				y_cand[j] = ub[j] - x[j]
			for id_j in range(id_m+1, n):
				j = argsort[id_j]				
				y_cand[j] = lb[j] - x[j]
			y_cand[m] = -np.sum(y_cand)
			if x[m] + y_cand[m] <= ub[m] and x[m] + y_cand[m] >= lb[m]:
				improvement = np.dot(gk, y_cand)
				if improvement < best_improvement:
					best_improvement = improvement
					y = y_cand
		if abs(best_improvement) < 1e-6:
			break
		## improvement phase 2
		s = np.dot(mu, y) - 2*lambdaval*np.matmul(np.matmul(x.T, cov),y)
		s = s / (2 * lambdaval * np.matmul(np.matmul(y.T, cov),y))

		if s < 0:
			s = 0
		if s > 1:
			s = 1
		x = x + s * y
	alldata['x'] = x
	return
	
def apply_factor(alldata, total_components):
	"""
	Apply factor model as described in factor.pdf to modify our QP program
	Spectra are computed using power method
	Parameters
	----------
		alldata: a dictionary of the following form:
			'n': size of the QP problem
			'lower': lower bound vector
			'upper': upper bound vector
			'x': solution vector
			'lambda': value of lambda
			'covariance': covariance matrix
		total_components: number of spectra to use
	Returns
	-------
		No returns
	"""
#==============================================================================
# 	e_values, e_vectors = np.linalg.eigh(alldata['covariance'])
# 	e_values = e_values[::-1]
# 	e_vectors = e_vectors[:, ::-1]
#==============================================================================
	cov = alldata['covariance']
	n = cov.shape[0]
	e_values, e_vectors = runpower(alldata['covariance'],n,tolerance=0, max_num=total_components )
	
	new_cov = np.zeros(alldata['covariance'].shape)
	for nth_component in range(total_components):
		new_cov += e_values[nth_component] \
			* np.outer(e_vectors[:,nth_component],e_vectors[:,nth_component])
	new_cov += np.diag(np.diag(alldata['covariance'] - new_cov))
	alldata['modified_covariance'] = new_cov

	
if __name__ == "__main__":
	if len(sys.argv) != 2:
		sys.exit("datafile")
	num_spectra = 1
	filename = sys.argv[1]
	print("input: ", sys.argv[1])
	
	alldata = readdata(filename)
    
	feasible(alldata)

	# use outside library
	Q = 2 * matrix(alldata['covariance'] * alldata['lambda'])
	p = -matrix(alldata['mu'])
	A = matrix(np.ones((1,alldata['n'])))
	b = matrix(1.0)
	G = matrix(np.vstack((np.eye(alldata['n']),-1*np.eye(alldata['n']))))
	h = matrix(np.hstack((alldata['upper'],-1*alldata['lower'])))
	
	sol=solvers.qp(Q, p, G, h, A, b)

	print(sol['x'], sol['primal objective'])
	
	# our algo
	num_iter = 100000

	optimize(alldata, num_iter)
	print(alldata['x'], eval_objective(alldata['lambda'], alldata['covariance'],
							alldata['mu'], alldata['x']))
		

	# Extra Credit					
	apply_factor(alldata, num_spectra)
	
	# use outside library
	Q = 2 * matrix(alldata['modified_covariance'] * alldata['lambda'])
	sol=solvers.qp(Q, p, G, h, A, b)
	
	print(sol['x'], sol['primal objective'])
	
	# our algo	
	
	optimize(alldata, num_iter, 'modified')
	print(alldata['x'], eval_objective(alldata['lambda'], alldata['modified_covariance'],
							alldata['mu'], alldata['x']))
							
	

	

