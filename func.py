# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 13:24:35 2018
Note: the following code is written in Python 3.5.2 not Python 2.7

"""
import math
import numpy as np
import pandas as pd
import sys

def eval_objective(lambdaval, cov, mu, x):
	"""
	Evaluate the objective function
	Parameters
	----------
		lambdaval: value of lambda
		cov: covariance matrix
		mu: mu vector
		x: current solution
	Returns
	-------
		value of the objective function at x
	"""
	tmp1 = np.matmul(np.matmul(x.T, cov),x)
	tmp2 = np.dot(mu,x)
	return lambdaval*tmp1 - tmp2


def feasible(alldata): 
	"""
	Checks whether the QP is feasible. If yes, provide a feasible solution,
	and store it in alldata['x'].
	Parameters
	----------
		alldata: a dictionary of the following form:
			'n': size of the QP problem
			'lower': lower bound vector
			'upper': upper bound vector
			'x': solution vector
			'lambda': value of lambda
			'covariance': covariance matrix

	Returns
	-------
		False if infeasible
		True if feasible found
	"""
	n = alldata['n']
	lower = alldata['lower']
	upper = alldata['upper']
	x = alldata['x']
	print(lower, upper)

	x = np.copy(lower)

	sumx = np.sum(x) 
	
	if sumx > 1.0:
		print("Infeasible.")
		return False

	for j in range(0,n):
		print("lower", lower, "upper",upper)
		print(j, 'sum', sumx, sumx + upper[j] -lower[j])
		if sumx + (upper[j] - lower[j]) >= 1.0:
			x[j] = 1.0 - sumx + lower[j]
			print('done')
			break
		else:
			x[j] = upper[j]
			delta = upper[j] - lower[j]
			print(x[j], lower[j], upper[j], delta)
			sumx += upper[j] - lower[j]
			print(">>>>",j, x[j], sumx )

	sumx = np.sum(x) 
	
	if sumx < 1.0:
		print("Infeasible.")
		return False
		
	print(x)
	alldata['x'] = x
	return True

def breakexit(foo):
	stuff = input("("+foo+") break> ")
	if stuff == 'x' or stuff == 'q':
		sys.exit("bye")

def readdata(filename):
	"""
	Read data from the file
	Parameters
	---------
		filename: string of filename path
	Returns:
		alldata: a dictionary of the following form:
			'n': size of the QP problem
			'lower': lower bound vector
			'upper': upper bound vector
			'lambda': value of lambda
			'covariance': covariance matrix
			'x': initialized solution vector, same as 'lower'
	"""
	# read data
	try:
		f = open(filename, 'r')
	except IOError:
		print ("Cannot open file %s\n" % filename)
		sys.exit("bye")
	lines = f.readlines()
	f.close()
	line0 = lines[0].split()
	if len(line0) == 0:
		sys.exit("empty first line")
	n = int(line0[1])
	print("n = ", n)
	
	lower = np.zeros(n)
	upper = np.zeros(n)
	mu = np.zeros(n)
	x = np.zeros(n)
	covariance = np.zeros((n,n))
	
	numlines = len(lines)
	#crude python
	linenum = 5
	while linenum <= 5 + n-1:
		line = lines[linenum-1]
		thisline = line.split()
		print(thisline)
		index = int(thisline[0])
		lower[index] = float(thisline[1])
		upper[index] = float(thisline[2])
		mu[index] = float(thisline[3])        
		linenum += 1
	linenum = n + 6
	line = lines[linenum-1]
	thisline = line.split()
	print(thisline)
	lambdaval = float(thisline[1])
	print("lambda = ", lambdaval)
	linenum = n + 10
	while linenum <= n+10 + n-1:
		line = lines[linenum-1]
		thisline = line.split()
		print(thisline)
		i = linenum - n - 10
		print(i)
		for j in range(n):
			covariance[i,j] = float(thisline[j])
		linenum += 1
	print(covariance)

	alldata = {}
	alldata['n'] = n
	alldata['lower'] = lower
	alldata['upper'] = upper
	alldata['mu'] = mu
	alldata['covariance'] = covariance
	alldata['lambda'] = lambdaval
	alldata['x'] = x

	return alldata
