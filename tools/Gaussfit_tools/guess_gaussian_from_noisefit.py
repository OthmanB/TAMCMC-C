import numpy as np
import subprocess

def guess_Amaxnumax_from_noisefit(params_file, dir_getmodel='../', ):
	'''
		Use a fit of the noise background in order to define guesses 
		for any excess of power that could be due to modes
	'''

	f=open(params_file,'r')
	txt=f.read()
	f.close()

	s=txt.split("\n")

	header=s[0]
	plength=s[1]
	params=s[2]

	# impose gaussian Height parameter to 0
	params[7]=0

	# Use the in-built function of TAMCMC to get the median model with noise only
	cmd=[dir_getmodel +'./getmodel ', data_file , params_cfg_noise , modelname , file_noise_fit]
	subprocess.run(cmd)