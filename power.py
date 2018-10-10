# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 13:24:35 2018
Note: the following code is written in Python 3.5.2 not Python 2.7

"""
import math
import numpy as np
import pandas as pd
import sys


def runpower_one(matrix, n):
	"""
	Calculate the leading eigenvalue and its corresponding eigenvector (normalized),
	using power method.
	Parameters
	----------
	matrix: 2D numpy array. 
		Assumed to be positive semi-definite.
	n: int
		size of the square matrix
	Returns
	-------
	eigenvalue: float
	eigenvector: np array. normalized so that L2 norm = 1.0
	"""
	#get initial vector
	v = np.zeros(n)
	w = np.zeros(n)
	for j in range(n):
		v[j] = np.random.uniform(0,1)
	#print 'matrix', matrix
	#print 'v', v
	T = 10000 #number of iterations
	tol = 1e-06
	oldnormw = 0
	for t in range(T):
		w = matrix.dot(v)
		#print 't', t, 'w',w
		normw = (np.inner(w,w))**.5
		v = w/normw
		#print 't',t,'v',v
		#print 't',t,'normw',normw, 'old', oldnormw
		if np.abs(normw - oldnormw)/normw < tol:
			#print ' breaking'
			break
		oldnormw = normw
	return normw, v
 

def runpower(matrix, n, tolerance, max_num=None, return_vector=True):
	"""
	Returns all the eigenvalues such that they are no smaller than a specific 
	fraction (specified by 'tolerance') than the leading eigenvalue.
	Calculation of eigenvalues is done using power method.
	Parameters
	----------
	matrix: 2D numpy array. 
		Assumed to be positive semi-definite.
	n: int
		size of the square matrix
	tolerance: float
		the tolerance e.g. 0.01
	max_num: int
		maximum number of eigenvalues to return, default=None
	return_vector: boolean
		whether eigenvectors will also be returned
	Returns
	-------
	list of eigenvalues in decreasing order
	"""
	calculate_next = True
	eigenvalue_list = []
	eigenvector_list = []
	leading_eigenvalue = np.nan
	while(calculate_next):	
		new_eigenvalue, v = runpower_one(matrix, n)
		if np.isnan(leading_eigenvalue):
			leading_eigenvalue = new_eigenvalue
		eigenvalue_list.append(new_eigenvalue)
		eigenvector_list.append(v)
		if max_num is not None and len(eigenvalue_list) == max_num:
			break
		if abs(1.0 * new_eigenvalue / leading_eigenvalue) < tolerance:
			calculate_next = False
		else:
			matrix = matrix - new_eigenvalue * np.outer(v,v)
	if return_vector:
		return eigenvalue_list, np.asarray(eigenvector_list).T
	else:
		return eigenvalue_list
