import numpy as np
'''
	A tiny function that allows you to give an estimate to the
	the slope and the ordinate at orgin for models fitting
	the power spectrum with a linear function of frequency
	This assumes you have already a guess estimate of the average
	aj coefficient and that you provide a guess on how much the value
	can depart from this average at most
	nu1, nu2: A minimum and maximum frequency range on which the estimate is made.
			  Should typically be your fitting range
	a_avg: The value of the average a-coefficient provided as a guess
	n: The multiplicative factor for a_avg that provides the maximum expected departure
	   For example, if you want the range that is associated to values of aj varying 
	   from 5*a_avg and a_avg/5, then n=5 
'''
def minmax(nu1, nu2, a_avg, n):
	k1=a_avg*(n-1./n)/(nu1-nu2)
	k2=-k1
	b_sols=[n*a_avg- k1*nu1, n*a_avg- k1*nu2, n*a_avg- k2*nu1, n*a_avg- k2*nu2]
	k_sols=[k1,k2]
	print(' solutions of linear type: a = k nu + b')
	print('slope range: [ {} , {} ]'.format(np.min(k_sols), np.max(k_sols)))
	print('ordinate @ origin range: [ {} , {} ]'.format(np.min(b_sols), np.max(b_sols)))


'''
	Same as minmax(), but with a direct range of a1,a2 provided
'''
def minmax_alt(nu1, nu2, a1, a2):
	k1=(a1-a2)/(nu1-nu2)
	k2=-k1
	b_sols=[a1- k1*nu1, a1- k1*nu2, a2- k2*nu1, a2- k2*nu2]
	k_sols=[k1,k2]
	print(' solutions of linear type: a = k nu + b')
	print('slope range: [ {} , {} ]'.format(np.min(k_sols), np.max(k_sols)))
	print('ordinate @ origin range: [ {} , {} ]'.format(np.min(b_sols), np.max(b_sols)))

	