'''
	Functions to ensure the calculation of the m-component visbilities for a given inclination and degree
'''

import numpy as np

def amplitude_ratio(l,beta):
	# beta is in degree
	angle=np.pi*beta/180.
	zz=function_rot(l,angle)
	zz=zz[:,l]**2
	norm=np.sum(zz)
	return zz/norm

def function_rot(l,beta):
	dim=2*l+1
	matr=np.zeros((dim,dim))
	for i in range(0,l+1):
	    for j in range(-i,i+1):
	    	matr[i+l,j+l]=dmm(l,i,j,beta)
	for i in range(-l,0+1):
		for j in range(i,-i+1):
			matr[i+l,j+l]=matr[-i+l,-j+l]*(-1e0)**(i-j)
	for j in range(0,l+1):
		for i in range(-j,j+1):
			matr[i+l,j+l]=dmm(l,j,i,-beta)
	for j in range(-l,0+1):
		if j+1 <= -j:
			step=-1
			maxi=j-1
		else:
			step=1
			maxi=j+1
		for i in range(-j, maxi, step):
			matr[i+l,j+l]=matr[-i+l,-j+l]*(-1e0)**(i-j)
	return matr

#************************
# compute dm1,m2(beta)  for  m1+m2>=0  and m1-m2>=0
#**************************
def dmm(l,m1,m2,beta):
	co=np.cos(beta/2e0)
	si=np.sin(beta/2e0)
	sum1=0e0
	if l-m1+1 < 0:
		step=-1
	else:
		step=1
	for s in range(0, l-m1+1, step):
		var=0e0
		var=combi(l+m2,l-m1-s)*combi(l-m2,s)*(-1e0)**(l-m1-s)
		var=var*co**(2e0*s+m1+m2)*si**(2e0*l-2e0*s-m1-m2)
		sum1=sum1+var
	sum1=sum1*np.sqrt(np.math.factorial(l+m1)*np.math.factorial(l-m1)*1e0)
	sum1=sum1/np.sqrt(np.math.factorial(l+m2)*np.math.factorial(l-m2)*1e0)
	return sum1

#************************************************
# Compute combinations n r  ,    n>r
#************************************************
def combi(n,r):
	c=np.math.factorial(n)/np.math.factorial(n-r)/np.math.factorial(r)
	return c