import numpy as np
import matplotlib.pyplot as plt
from function_rot import *

def lorentzian_cpp(func_name, params, xmin, xmax, raw=False):
	try:
		params_str=""
		for p in params:
			params_str=params_str + " p"
		process = Popen(["./lorentzian_test", func_name, str(xmin), str(xmax), params_str], stdout=PIPE, stderr=PIPE)
		(output, err) = process.communicate()
		exit_code = process.wait()
		output=output.decode("utf-8") 
		if raw == False:
			r=output.split('\n')
			freq=[]
			power=[]
			for line in r:
				line=line.strip()
				try:
					if line[0] != "!" and line[0] != "#":
						s=line.split()
						freq.append(float(s[0]))
						power.append(float(s[1]))
				except:
					err=True
			return freq, power
		else:
			return output, err
	except: # Handling processes that does not exist, aka, when lorentzian_test file is not available
		error=True
		print("Error: Could not execute the lorentzian_test C++ program. The most likely explanation is that it is not in the current directory")
		return -1, error

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
			func_names.append(s[0]) # extract the name of the function
			for p_name in s[1:]: # extrract the names of all of the parameters and remove all potential spaces for the function in the line
				param_names.append(p_name.strip())
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
    rho_sun=M_sun*1e3/(4*pi*(R_sun*1e5)**3 /3); #in g.cm-3
    Dnl=0.75; 
   	# -------------------
	Dnu_star=135.1
	rho=(Dnu_star/Dnu_sun)**2. * rho_sun;
    # -------------------
	H_l=1    									# Maximum height
	fc_l=50  									# Central frequency
	a1=10    									# a1 coeficient
	eta0=(4./3.)*np.pi/(rho*G);  				# coefficient of distorsion due to the centrifugal force  !! WARNING WILL NEED CHECKS DUE TO Dnl missing here
	eta=eta0*Dnl*Qlm*(a1*1e-6)**2  				# coeficient of distorsion due to the centrifugal force (alternate form that includes a1) !!! WARNING: CHECK Dnl and Qlm in the THEORY... MAY BE THE SAME
	a2=a1/100									# a2 coeficient
	a3=a1/20			    					# a3 coeficient
	a4=a1/20									# a4 coeficient
	a5=a1/50			    					# a5 coeficient
	a6=a1/50                					# a6 coeficient
	f_s=a1                  					# rotational splitting (same as a1)
	asym=50			        					# Asymetry of the Lorentzian normalised by width
	gamma_l=2               					# Width of the mode
	inclination=45          					# Stellar inclination
	V=function_rot(l, inclination)          	# Mode relative height. Varies with l and inclination
	handled_names=['H_l', 'fc_l', 'a1', 'eta', 'eta0', 'a2', 'a3', 'a4', 'a5', 'a6', 'f_s', 'asym', 'gamma_l', 'inclination', 'V', 'l']
	return handled_names, [H_l, fc_l, a1, eta, eta0, a2, a3,a4, a5, a6, f_s, asym, gamma_l, inclination, V, l]

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
		while found == False or i<len(handled_names):
			if p == handled_names[i]:
				found=True
			i=i+1
		if found == True:
			params.append(default_params[i-1])
		else:
			print("Error in get_default_params while identifying the set of parameter to be used to test the function")
			print("      Please check that the 'default_test_arguments() contains the expected arguments for the function func_name that is attempted to compute")
			print("      Expected parameter names are: ", param_names)
			print("      Default list of arguments   : ", handled_names)
			print("      Execution stoped because it could not find the following parameter name: ", p)
			print("      The program will exit now")
			exit()
	return params

def set_expectations(params, param_names, func_name):
	# This function allows to set what is the expectation value of a given
	# parameter directly in terms of frequency or power such that it can be
	# visualise/plot directly over the power spectrum
	print("TBD")

def visualise(freq, power, params):
	# Function that take the power spectrum that is created with the lorentzian_cpp() function
	# and visualise it along with the expected frequency positions as per
	# defined by the params vector
	plt.plot(freq, power)



# Main test program section:
# 
# -----  Initial setup -----
# Frequency range on which we will draw the result
xmin=0
xmax=100
# Degrees that we want to test:
l_all=[0,1,2,3]
# --------------------------
#
# Retrieve all of the information on the models that need to be tested
func_names, param_names_all=get_func_infos()
#
# We will loop on all of the models that are listed in func_names
for i in range(len(func_names)):
	for l in l_all:
		# Identify the values that we will need to use in function of the considered model
		# For this we use the information provided by get_func_infos about the required parameters
		params=get_default_params(param_names_all[i],l)
		# Compute the Lorentzian profile using the C++ program and according to the specified variables
		freq, power=lorentzian_cpp(func_name, params, xmin, xmax) 
		# Do plots that show the profile and the position of frequencies