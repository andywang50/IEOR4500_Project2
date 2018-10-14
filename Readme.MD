IEOR 4500 Fall 2018 Project 2
Yashin Chen (yc3347), Chandana Dayapule (cd3037), Zhengrui Fu (zf2181), Pranathi Pothireddy(pp2665), Guoan Wang (gw2355)

python part:

	files needed: example.txt, func.py, power.py, Question1.py
	
	in console, python Question1.py example.txt
	
	First, use algo2 to solve the QP problem; the result is then compared to an outside library cvxopt.
	
	optimal turned out to be 11.915222033156242 (approximately).
	
	Solution:
	[0.05469998 0.10939043 0.43591295 0.39999664]
	
	Extra credit part:
	
	using only 1 spectra when doing the pca. (the number of spectra can be changed by changing the "num_spectra" variable in line 146 in Question1.py
	
	PCA is done using power method as in Project 1.
	
	Again, our result is compared to cvxopt library.
	
	optimal:
	9.364564918142761
	
	Solution:
	[0.0501096  0.13840027 0.41149188 0.39999826]
	
	
C part:
	files needed: example.txt, algo.h, algo.cpp, power.h, power.cpp, linalg.h, linalg.cpp, main.cpp
	
	Please use Visual Studio to run the code.
	
	PCA is done using power method as in Project 1 (we converted our Project 1 python code to C as well).
	
	Solutions in Question 1 and Question 2(Extra Credit) are approximately the same as in python, but the program runs way faster.
