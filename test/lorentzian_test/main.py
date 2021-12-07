import numpy as np
import matplotlib.pyplot as plt
from function_rot import *
from acoefs import nunlm_from_acoefs
from subprocess import Popen, PIPE
from iteration_utilities import deepflatten

def Qlm(l,m):
    Qlm=(l*(l+1) - 3*m**2)/((2*l - 1)*(2*l + 3))
    Qlm=Qlm*2./3
    return Qlm

def lorentzian_cpp(func_name, params, xmin, xmax, raw=False):
	err=False
	print('func_name  =', func_name)
	print('xmin       =', xmin)
	print('xmax       =', xmax)
	print('params     =', params)
	try:
		cmd=["./lorentzian_test", func_name, str(xmin), str(xmax)]
		#params_str=""
		for p in params:
			#params_str=params_str + "  " + str(p)
			cmd.append(str(p))
		print('cmd    =', cmd)
		process = Popen(cmd, stdout=PIPE, stderr=PIPE)
		(output, err) = process.communicate()
		exit_code = process.wait()
		output=output.decode("utf-8") 
		if raw == False:
			r=output.split('\n')
			freq=[]
			power=[]
			for line in r:
				line=line.strip()
				#print('B:', line,   ' len(line) =', len(line))
				try:
					if len(line) != 0:
						if line[0] != "!" and line[0] != "#":
							s=line.split()
							freq.append(float(s[0]))
							power.append(float(s[1]))
					else:
						print('Warning: line with no value returned (Ignored)')
				except:
					print('Raised an exception in lorentzian_cpp()')
					err=True
			return freq, power, err
		else:
			return output, err
	except: # Handling processes that does not exist, aka, when lorentzian_test file is not available
		error=True
		print("Error: Could not execute the lorentzian_test C++ program. The most likely explanation is that it is not in the current directory")
		return [], [], error

def get_func_infos():
	# Function that run a c++ program lorentzian_test in order to retrieve all function names that are hanlded by the program
	# It also returns a 2D array that provides the name of the inputs for each of those functions
	process = Popen(["./lorentzian_test", "get_all"], stdout=PIPE, stderr=PIPE)
	(output, err) = process.communicate()
	exit_code = process.wait()
	output=output.decode("utf-8") 
	r=output.split('\n')
	func_names=[]
	param_names_all=[]
	for line in r:
		param_names=[]
		line=line.strip() # remove white spaces at the end and at the begining
		try:
			s=line.split() # split according to spaces
			if s[0] != '#':
				func_names.append(s[0].strip()) # extract the name of the function
				for p_name in s[1:]: # extrract the names of all of the parameters and remove all potential spaces for the function in the line
					#print('p_name[-1]=', p_name[-1])
					if p_name[-1] != ',':
						param_names.append(p_name.strip())
					else:
						#print('len(p_name)=', len(p_name))
						param_names.append(p_name.strip()[0:len(p_name)-1])
				param_names_all.append(param_names) # Regroup in a list the list of names for the specific function (list of lists)
		except:
			err=True
	return func_names, param_names_all

def default_test_arguments(l):
	# function that list all of the arguments that can be fed into the 'lorentzian_test' C++ function
	# and set their default values. Names are set identical to their counterpart in the C++ function (see './lorentzian_test get_all' for details)
	# Usual function for the definition for the Lorentzian components in a global fit:
	# ---- Constants ----
	G=6.667e-8;
	Dnu_sun=135.1;
	R_sun=6.96342e5; #in km
	M_sun=1.98855e30; #in kg
	rho_sun=M_sun*1e3/(4*np.pi*(R_sun*1e5)**3 /3); #in g.cm-3
	# -------------------
	Dnu_star=135.1
	rho=(Dnu_star/Dnu_sun)**2. * rho_sun;
	# -------------------
	H_l=1    									# Maximum height
	fc_l=1000  									# Central frequency
	a1=10    									# a1 coeficient
	eta0=3./(4.*np.pi*G*rho)					# coefficient of distorsion due to the centrifugal force: eta0*Qlm(l,m) with Qlm = Dnl * [l(l+1) - 6m^2]/
	#eta=eta0*Dnl*Qlm(l,m)*(a1*1e-6)**2  		# coeficient of distorsion due to the centrifugal force (alternate form that includes a1) !!! WARNING: CHECK Dnl and Qlm in the THEORY... MAY BE THE SAME
	a2=a1/100									# a2 coeficient
	a3=a1/50			    					# a3 coeficient
	a4=a1/50									# a4 coeficient
	a5=a1/100			    					# a5 coeficient
	a6=0#a1/100                					# a6 coeficient
	f_s=a1                  					# rotational splitting (same as a1)
	asym=0 #50			        				# Asymetry of the Lorentzian normalised by width
	b=5. 				 						# Model with some activity in the form of a power law (Function that is not used anymore)
	alpha=1. 									# Model with some activity in the form of a power law (Function that is not used anymore)
	gamma_l=2               					# Width of the mode
	inclination=45          					# Stellar inclination
	V=amplitude_ratio(l, inclination)          	# Mode relative height. Varies with l and inclination
	handled_names=['H_l', 'fc_l', 'a1', 'eta0', 'a2', 'a3', 'a4', 'a5', 'a6', 'f_s', 'asym', 'gamma_l', 'inclination', 'V', 'l']
	return handled_names, [H_l, fc_l, a1, eta0, a2, a3,a4, a5, a6, f_s, asym, gamma_l, inclination, V, l]

def get_default_params(param_names, l):
	# Function that retrieve the default list of parameters and use it to build a suitable set of parameters 
	# for the requested model function (func_name) and as per specified the parameter names which can be retrieved using ./lorentzian_test get_all
	#
	handled_names, default_params=default_test_arguments(l)
	#
	params=[]
	for p in param_names:
		i=0
		found=False
		while found == False and i<len(handled_names):
			#print(p, '    ', handled_names[i],  '   bool: ', p == handled_names[i])
			if p == handled_names[i]:
				found=True
			i=i+1
		if found == True:
			params.append(default_params[i-1])
			#print("params: ", params)
			#print("      default_params[", i-1, "] =", default_params[i-1])
		else:
			print("Error in get_default_params while identifying the set of parameter to be used to test the function")
			print("      Please check that the 'default_test_arguments() contains the expected arguments for the function func_name that is attempted to compute")
			print("      Expected parameter names are: ", param_names)
			print("      Default list of arguments   : ", handled_names)
			print("      Execution stoped because it could not find the following parameter name: ", p)
			print("      The program will exit now")
			exit()
	return params


def get_params(param, param_names, param_class):
	# A function that ensure that we get all the variables that are necessary to compute
	# the visualisation are properly organized and provided. 
	# Ensure that all the required variable are identified and organised in a way that can be handled
	# Note that this function does not check whether there is a problem on the parameters
	# It only checks if param_class is adequately defined. Any parameter that was not found into 
	# the list of parameters will keep its default value, which can be used to detect errors
	error_class=True
	if param_class == 'Height':
		req_param_names=['l', 'H_l', 'H_lm', 'inclination']
		req_params=[-1, -1, [-1], -1]
		error_class=False
	if param_class == 'fc_l':
		req_param_names=['fc_l']
		req_params=[-1]
		error_class=False
	if param_class == 'splitting_a1etaa3': # Handling of all splitting models that involve a1, eta and a3
		req_param_names=['l', 'f_s', 'eta0', 'a3', 'b', 'alpha'] 
		req_params     =[-1,  -1,     -1   ,  -1 , -1, -1]
		error_class=False
	if param_class == 'splitting_a1a2a3':
		req_param_names=['l', 'f_s',  'a2',   'a3']
		req_params     =[-1 , -1   ,   -1 ,    -1]
		error_class=False
	if param_class == 'splitting_aj':
		req_param_names=['l', 'a1', 'a2',   'a3', 'a4', 'a5', 'a6',  'eta0']
		req_params     =[-1 , -1   ,   -1 ,  -1  ,  -1 ,  -1 ,  -1 ,  -1]
		error_class=False
	if param_class == 'splitting_Alm':
		req_param_names=['l', 'a1', 'a3' 'a5',  'theta0', 'delta', 'filter_type']
		req_params     =[-1 , -1   ,   -1 ,  -1  ,  -1 ,  -1, None ]
		error_class=False
	if param_class == 'splitting_a1l':
		req_param_names=['l', 'a1',  'a2',   'a3', 'eta0']
		req_params     =[-1 , -1   ,   -1 ,    -1,   -1]
		error_class=False
	if error_class == True:
		print("Error : Could not find the parameter in get_params: ")
		print("params     : ", param)
		print("param_names: ", param_names)
		print("param_class: ", param_class)
		exit()
	for i in range(len(req_param_names)):
		for j in range(len(param_names)):
			if param_names[j] == req_param_names[i]:
				req_params[i]=param[j]
	return req_params

def splittings_checks(split_params):
	error=False
	for s in split_params: # Verify that all parameters are provided
		if s == -1:
			error=True
			print('Warning: Error in split_params')
	return error

def Hlm_checks(l, H_l=-1, H_lm=[-1], inclination=-1):
	error=True
	if H_lm[0] == -1 and inclination !=-1: # Case where only inclination is provided
		error=False
		V=amplitude_ratio(l, inclination)
		Hlm=V*H_l
	if H_lm[0] != -1 and inclination == -1: # Case where the Visibilities are directly provided
		# Nothing to do except noting that there is no error
		error=False
	if error == True: # Any other scenario has a problem
		print('Error: Unexpected combination of parameter. ')
		print("       You must provide Hlm directly OR H_l + inclination")
		print("       The program will exit now")
		exit()
	return Hlm

def do_Hlm(params, param_names):
	# Performs the computation of Hlm depending on the circounstances
	l,H_l, H_lm, inc=get_params(params,param_names, 'Height')
	Hlm=Hlm_checks(l, H_l=H_l, H_lm=H_lm, inclination=inc)
	return Hlm

def do_fc_l(params, param_names):
	# Retrieve the central frequency
	fc_l=get_params(params,param_names, 'fc_l')
	if fc_l == -1:
		print('Error with the central frequency fc_l: It has to be set to a non-zero, positive value')
		print('      The program will exit now')
		exit()
	return fc_l

def nu_nlm(nu_nl, l, a1=0, a2=0, a3=0,a4=0,a5=0,a6=0, eta0=0, theta0=None, delta=None, filter_type=None):
	# Compute the nu(n,l,m) frequencies in all of the scenarii that are implemented in build_lorentzian routine.cpp
	# There is some fail checks in order to avoid to provide incompatible quantities
	# Allowed values:
	#		aj coefficients + eta0
	# 	OR
	#   	odds aj coefficients + theta0 AND delta AND filter_type 
	AR_term=False # The default is the Active Region effect is not implemented by lack of arguments
	#Dnl=2./3
	#
	if (a2 != 0 or a4 !=0 or a6 !=0) or eta0 !=0: # detect if any of the even coefficent is non-zero or if the Centrifugal force parameter is implemented
		forbid_Alm=True  # if it is the case, we flag as forbiden any model that includes theta0 and delta
	#
	error=False
	if forbid_Alm == False:
		if theta0 != None and delta != None and filter_type != None:
			AR_term=True
		else:
			error=True
			print('Error: You must provide theta0 AND delta AND filter_type inside nu_nlm')
			print('       The computation of the AR_term cannot be performed and the parameter will be ignored')
			#AR_term=False
	else:
		if theta0 != None or delta !=None or filter_type !=None:
			error=True
			print('Warning: Configuration of aj coefficients and/or eta, eta0 incompatible with Alm models detected')
			print('         The Alm coefficients will be ignored and set to 0')
			print("   nu_nl      :", nu_nl)
			print("   a1         : ", a1)
			print("   a2         : ", a2)
			print("   a3         : ", a3)
			print("   a4         : ", a4)
			print("   a5         : ", a5)
			print("   a6         : ", a6)
			print("   theta0     : ", theta0)
			print("   delta      : ", delta)
			print("   filter_type: ", filter_type)
			#AR_term=False
	nu_aj=nunlm_from_acoefs(nu_nl[0], l, a1=a1, a2=a2, a3=a3, a4=a4,a5=a5,a6=a6) # Use the acoefs.py core freq routine to get frequencies with the aj decomposition
	for m in range(-l, l+1):
		CF_term=eta0*Qlm(l,m)*(a1*1e-6)**2 
		if AR_term == True:
			AR_term=epsilon_nl*Alm(l, m, thetas[0], thetas[1], filter_type)	
		else:
			AR_term=0
		print("AR_term :", AR_term)
		print("CF_term :", CF_term,   '   eta0 =', eta0,  '    a1 =', a1)
		print("         Corresponding DR/R = ", eta0*(a1*1e-6)**2 * 2*np.pi**2 * 1e5)  # See my notes on Ipad for the derivation
		nu_aj[m+l]=nu_aj[m+l]+(CF_term + AR_term)*nu_nl[0]
	print("nu_aj :", nu_aj)
	print(" --------- ")
	#exit()
	return nu_aj, error

def do_splittings(params, param_names, func_name):
	# Handling of all the cases for splittings
	# Performs the computation of the frequencies for all modes depending on circounstances
	param_class=parameter_class(func_name)
	split_params=get_params(params,param_names, param_class)
	fc_l=do_fc_l(params, param_names)
	error=True
	if split_params[0] == -1:
		print('Error in do_splittings: The degree l was not provided')
		print('      split_params: ', split_params)
		print('      Please check the requirements in the program')
		print('      The program will stop now')
		exit()
	else:
		l=split_params[0]
		split_lm=np.zeros(2*l+1)
	#
	if param_class == 'splitting_a1etaa3': # Deal with cases a1,eta,a3 and also with the obselete function a1acta3
		split_lm, err=nu_nlm(fc_l, l, a1=split_params[1], a2=0, a3=split_params[3],a4=0,a5=0,a6=0, eta0=0, theta0=None, delta=None, filter_type=None)
		if (split_params[1] != -1) and (split_params[2] == -1) and (split_params[3] != -1) and (split_params[4] == -1) and (split_params[5] == -1): # models with a1, eta, a3 only
			error=False
		if (split_params[1] != -1) and (split_params[2] == -1) and (split_params[3] != -1) and (split_params[4] != -1) and (split_params[5] != -1): # models with a1, eta, a3 + b and alpha (obselete model)
			error=False
			b=split_params[4]
			alpha=split_params[5]
			for m in range(-l,l+1):
				split_lm[m+l]=split_lm[m+l] + Qlm(l,m)*b*pow(fc_l*1e-3,alpha)
		print("param_class == ", param_class, " IS NOT TESTED. YOU NEED TO WRITE/CHECK CODE FOR IT" )
		exit()
	#
	if param_class == 'splitting_a1etaAlma3':
		split_lm, err=nu_nlm(fc_l, l, a1=split_params[1], a2=0, a3=split_params[2],a4=0,a5=split_params[3],a6=0, eta0=0, theta0=split_params[4], delta=split_params[5], filter_type=split_params[6])
		print("param_class == ", param_class, " IS NOT TESTED. YOU NEED TO WRITE/CHECK CODE FOR IT" )
		exit()
	#
	if param_class == 'splitting_a1a2a3':
		error=splittings_checks(split_params)
		if error == False: # We can proceed only if all of the parameters were provided
			split_lm, err=nu_nlm(fc_l, l, a1=split_params[1], a2=split_params[2], a3=split_params[3],a4=0,a5=0,a6=0, eta0=0, theta0=None, delta=None, filter_type=None)
		print("param_class == ", param_class, " IS NOT TESTED. YOU NEED TO WRITE/CHECK CODE FOR IT" )
		exit()
	#
	if param_class == 'splitting_aj':
		error=splittings_checks(split_params)
		if error == False: # We can proceed only if all of the parameters were provided
			split_lm, err=nu_nlm(fc_l, l, a1=split_params[1], a2=split_params[2], a3=split_params[3],a4=split_params[4],a5=split_params[5],a6=split_params[6], eta0=split_params[7], theta0=None, delta=None, filter_type=None)
	#
	if param_class == 'splitting_a1l_etaa3': # Handle also the case _v2 which are about H_lm instead of H_l fitting
		if split_params[0] == 1:
			f_s=f_s1;
		if split_params[0] == 2:
			f_s=f_s2;
		if split_params[0] == 3:
			f_s=(f_s1 + f_s2)/2. # APPROXIMATION
		if split_params[0] > 3:
			print('Error in do_splittings(): l>3 is not allowed ')
			exit()
		split_lm,err=nu_nlm(fc_l, l, a1=split_params[1], a2=0, a3=split_params[3],a4=0,a5=0,a6=0, eta0=split_params[2], theta0=None, delta=None, filter_type=None)
		#split_lm=fc_l*(1. + eta*Qlm) + Pslm(1,l,m)*f_s + Pslm(3,l,m)*a3
		print("param_class == ", param_class, " IS NOT TESTED. YOU NEED TO WRITE/CHECK CODE FOR IT" )
		exit()
	#
	if param_class == 'splitting_a1l_a2a3':
		if l != 0:
			clm=Pslm(3,l,m)
			a2_terms=Pslm(2,l,m)*a2;
		if l == 1:
			f_s=f_s1;
		if l == 2:
			clm=(5*pow(m,3) - 17*m)/3.# a3 for l=2
			f_s=f_s2;
		if l == 3:
			clm=(pow(m,3)-7*m)/2 # a3 implemented on 30/04/2021
			f_s=(f_s1 + f_s2)/2. # APPROXIMATION      
		split_lm,err=nu_nlm(fc_l, l, a1=split_params[1], a2=0, a3=split_params[3],a4=0,a5=0,a6=0, eta0=split_params[2], theta0=None, delta=None, filter_type=None)
		#    split_lm=fc_l + m*f_s + a2_terms + clm*a3
		print("param_class == ", param_class, " IS NOT TESTED. YOU NEED TO WRITE/CHECK CODE FOR IT" )
		exit()
	if error == True:
		print('Error in do_splittings: Some required parameters for the splitting were not provided.')
		print('      split_params: ', split_params)
		print('      Please check the requirements in the program')
		exit()
	return split_lm, error # MUST BE A VECTOR OF SIZE 2l+1

def parameter_class(func_name):
	param_class='-1'
	if func_name == 'optimum_lorentzian_calc_a1etaa3' or func_name == 'optimum_lorentzian_calc_a1acta3':
		param_class='splitting_a1etaa3'
	if func_name =='optimum_lorentzian_calc_a1etaAlma3':
		param_class='splitting_a1etaAlma3'
	if func_name =='optimum_lorentzian_calc_a1a2a3':
		param_class='splitting_a1a2a3'
	if func_name =='optimum_lorentzian_calc_aj':
		param_class='splitting_aj'

	if param_class == "-1":
		print("Error in parameter_class(): Could not found the function name in the list of provided functions ")
		print("func_name:", func_name)
		print("param_class :", param_class)
	return param_class

def set_expectations(params, param_names, func_name):
	# This function allows to set what is the expectation value of a given
	# parameter directly in terms of frequency or power such that it can be
	# visualise/plot directly over the power spectrum
	error=True
	Hlm=do_Hlm(params, param_names)

	fc_l=do_fc_l(params, param_names)
#	print(param_class)
#	print(func_name)
	nu_split,error=do_splittings(params, param_names, func_name)
	#
	if error == True:
		print("Error in main.py: Could not find the function in the preset list of functions.")
		print("                  Cannot interpret inputs adequately. The program will exit now")
		exit()
	#
	xout=[fc_l, nu_split]
	xout=list(deepflatten(xout))
	yout=list(deepflatten(Hlm))
	xnames=['fc_l']
	for i in range(len(nu_split)):
		xnames.append('nu_split')	
	ynames=np.repeat('Hlm',len(Hlm))
	return xout, yout, xnames,ynames

def visualise(freq, power, params, xlines,ylines, xnames, ynames):
	# Function that take the power spectrum that is created with the lorentzian_cpp() function
	# and visualise it along with the expected frequency positions as per
	# defined by the params vector
	plt.plot(freq, power)

	#print("xlines:", xlines)
	#print("ylines:", ylines)
	#exit()
	for i in range(len(xlines)):
		x=xlines[i]
		if  xnames[i] == 'fc_l':
			col='green'
		else:
			col='blue'		
		plt.axvline(x=x, color=col, linestyle='--')
	for y in ylines:
		plt.axhline(y=y, color='red', linestyle='--')
	plt.show()


def main():
	# Main test program section:
	# 
	# -----  Initial setup -----
	# Frequency range on which we will draw the result
	xmin=950
	xmax=1050
	# Degrees that we want to test:
	l_all=[0,1,2,3]
	# --------------------------
	#
	# Retrieve all of the information on the models that need to be tested
	func_names, param_names_all=get_func_infos()
	#print('Available function names:')
	#for func in func_names:
	#	print('  ', func)
	#
	do_funcs=np.repeat(False, len(func_names)) # Variable that allows to control which models to test: Default is NO tests
	#do_funcs=np.repeat(True, len(func_names)) # Variable that allows to control which models to test: Default is ALL tests
	#
	do_funcs[func_names.index('optimum_lorentzian_calc_aj')]=1 # Single model testing
	#
	# We will loop on all of the models that are listed in func_names
	for i in range(len(func_names)):
		if do_funcs[i] == False:
			print(func_names[i], ' Skiping...')
		if do_funcs[i] == True: # Do only requested tests
			print(func_names[i], ' Testing...')
			for l in l_all:
				# Identify the values that we will need to use in function of the considered model
				# For this we use the information provided by get_func_infos about the required parameters
				params=get_default_params(param_names_all[i],l)
				# Compute the Lorentzian profile using the C++ program and according to the specified variables
				freq, power, err=lorentzian_cpp(func_names[i], params, xmin, xmax) 
				# Get information about the inputs and format them so that they can be plotted
				xlines, ylines, xnames, ynames=set_expectations(params, param_names_all[i], func_names[i])
				# Do plots that show the profile and the position of frequencies
				visualise(freq, power, params, xlines, ylines, xnames, ynames)
