# This file contains functions that can be used to generate 
# input files (.model and .data) for the analysis of the power spectrum
# by fitting a Gaussian mode envelope
import numpy as np
from scipy.io import readsav
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

# ------------------
# function that calculate initial guesses for noise and for the First fitting Step
def inital_param_S0(numax, Amp, freq, spec_reg, fmin, fmax):
	init_param=[]
	err=False

	Teff_sun=5777.
	Dnu_sun=135.
	numax_sun=3100.
	a0=0.77 # Stello et al. (2009) for the relation between Dnu and numax
	# ---- Values for the relation between Gamma and numax (Benomar 2013, Eq.1)
	#alpha=double(0.267) & beta=double(0.76) & aa=4.06d-42 & po=double(11)
	#K=alpha^(8./3.) * numax_sun^2 * Teff_sun / Dnu_sun^(8./3.)
	# ----------------- These are used as rough estimates ---------------------
	# ---- This is the average M/Teff/L of the observed stellar population ----
	M=1.25
	Teff=6000
	L=1.5
	b0=M**(0.5 -a0) * (Teff/Teff_sun)**(3-3.5*a0) / L**(0.75-a0)
	# -------------------------------------------------------------------------

	# ------- Setting the smooth coef for noise characterisation -------
	Dnu=Dnu_sun * b0 * ( numax/numax_sun )**a0 # Stello et al. (2009) Defines the Heavy smooth coef 
	coef=0.4 #1.2 # Defines the average smoothing (smooth=coef * Dnu/resol)
	tol=3 # used to remove modes of the PSF (removed window= tol *Dnu)
	# ------------------------------------------------------------------

	if (np.min(freq) > 10):
		print, 'Warning: This function was designed to receive full spectra (begining at 0)'
		print, 'Provide such a spectra so that we can proceed'
		print, 'The program will stop now'
		#exit()
		err=True
		return [], err

	resol=freq[2] - freq[1]
	smoothcoef=Dnu*coef/resol # in terms of sigma

	s_smooth=gaussian_filter1d(spec_reg, smoothcoef, mode='reflect', truncate=4) # a properly smooth spectra 
	s_smooth=s_smooth[np.where(np.bitwise_and(freq >= fmin, freq <= fmax))]
	f_smooth=freq[np.where(np.bitwise_and(freq >= fmin , freq <= fmax))]
	s_smooth=s_smooth[np.where(np.bitwise_or(freq <= (numax - Dnu*tol), freq >= (numax + Dnu*tol)))]
	f_smooth=f_smooth[np.where(np.bitwise_or(freq <= (numax - Dnu*tol), freq >= (numax + Dnu*tol)))]

	Fnyquist=np.max(freq)

	# **** estimate the harvey profile parameters ****

	if Fnyquist >= 1000.:
		stype='SC'
	else:
		stype='LC'

	# ---- B0 ----
	# SHORT CADENCE
	if stype == 'SC' :
		pos_0=np.min( np.where( freq >= (Fnyquist - 500.) )) # we take the last 500microHz to extract B0
		dnu_win=15. #15 microHz average
	# LONG CADENCE
	if stype == 'LC':
		pos_0=np.min( np.where( freq >= Fnyquist - 10.)) # we take the last 10microHz to extract B0
		dnu_win=5. #15 microHz average
	B0=np.mean(spec_reg[pos_0:])
	# ------------

	
	# ---- H1 ----
	# SHORT CADENCE
	nu1=np.min(f_smooth);
	pos=np.where(np.bitwise_and(f_smooth >= nu1 , f_smooth <= (nu1 + dnu_win) ))
	if stype == 'SC': 
		H1=np.mean(s_smooth[pos]) - B0 # VALID ONLY WHEN THE PSF IS NOT BF FILTERED !!!!
		increment=10.
	# LONG CADENCE
	if stype == 'LC':
		H1=np.max(s_smooth[pos]) - B0 # VALID ONLY WHEN THE PSF IS NOT BF FILTERED !!!!
		increment=3.
	# -------------

	# ---- tc1 ---
	nu2=np.min(f_smooth) + increment # second sample taken 10 microHz apart from 
	pos=np.where(np.bitwise_and(f_smooth >= nu2 , f_smooth <= nu2 + dnu_win))
	H1_2=np.mean(s_smooth[pos]) - B0
	cpt=0.
	while (H1_2/H1 >= 0.5) and cpt <= 60: # as soon as the intensity has not decreased by at least half...
		nu2=nu2 + increment
		pos=np.where(np.bitwise_and(f_smooth >= nu2 , f_smooth <= nu2 + dnu_win))
		H1_2=np.mean(s_smooth[pos]) - B0
		cpt=cpt+1
	p1=4.
	ac1=H1/H1_2
	tc1=((ac1 -1)/( ((nu2 + dnu_win/2.)*1e-3)**p1 - ac1*((nu1 + dnu_win/2.)*1e-3)**p1))**(1./p1)


	# ----- H2 -----
	if stype == 'SC':
		fmin=200 # remove lowest frequencies as we look for the second harvey profile
	if stype == 'LC':
		fmin=(3 - 1)**(1./p1) / tc1 * 1e3
	s_smooth=s_smooth[np.where(np.bitwise_and(f_smooth >= fmin , f_smooth <= fmax))]
	f_smooth=f_smooth[np.where(np.bitwise_and(f_smooth >= fmin , f_smooth <= fmax)) ]
	nu1=np.min(f_smooth);
	pos=np.where(np.bitwise_and(f_smooth >= nu1, f_smooth <= nu1 + dnu_win))
	if stype == 'SC':
		H2=np.mean(s_smooth[pos]) - B0 # VALID ONLY WHEN THE PSF IS NOT BF FILTERED !!!!
	# LONG CADENCE
	if stype == 'LC':
		H2=np.max(s_smooth) - B0 # VALID ONLY WHEN THE PSF IS NOT BF FILTERED !!!!

	# ---- tc2 ---
	# SHORT CADENCE
	if stype == 'SC':
		increment=25. # The incremental step to evaluate the noise background can be large for short cadence because the frequency range to explore is large
	# LONG CADENCE
	if stype == 'LC':
		increment=2 # The incremental step to evaluate the noise background must be small for long cadence
		nu2=np.min(f_smooth) + increment # second sample taken 50 microHz apart from 
		pos=np.where(np.bitwise_and(f_smooth >= nu2 , f_smooth <= nu2 + dnu_win))
		H2_2=np.mean(s_smooth[pos]) - B0
		while (H2_2/H2 <= 0.5) or H2_2 <=0:
			nu2=np.min(f_smooth) + increment # second sample taken 50 microHz apart from 
			pos=np.where(np.bitwise_and(f_smooth >= nu2 , f_smooth <= nu2 + dnu_win))
			H2_2=np.mean(s_smooth[pos]) - B0
			increment=increment*2
		increment=increment/2	
	cpt=0.
	force_exit=0 # Boolean switch that turns into 1 if Fnyquist is passed
	while (H2_2/H2 >= 0.5) and cpt < 40 and (force_exit == 0): # as soon as the intensity has not decreased by at least half...
		if nu2+increment <= Fnyquist:
			nu2=nu2 + increment
			pos=np.where(np.bitwise_and(f_smooth >= nu2 , f_smooth <= nu2 + dnu_win))
			H2_2=np.mean(s_smooth[pos]) - B0
		else:
			force_exit=1 # If Nyquist is reached without reaching half power, we consider That the curvature for H(nu) at H2 occurs at H(numax)=Hnumax. See a bit below
		cpt=cpt+1
	if force_exit == 0: # Case np.where we could find H2 without error
		p2=2.
		ac2=H2/H2_2
		tc2=((ac2 -1)/( ((nu2 + dnu_win/2.)*1e-3)**p2 - ac2*((nu1 + dnu_win/2.)*1e-3)**2))**(1./p2)
	else:  # Case np.where H2 was not found within the 0 - Fnyquist ==> Consider That the curvature for H(nu) at H2 occurs at H(numax)=Hnumax
		p2=2.5
		pos0=np.where(np.bitwise_and(freq >= numax-3*Dnu , freq <= numax+3*Dnu) ) # Use the expected Dnu to get the average power range
		Hnumax=np.mean(spec_reg[pos0])
		H2=Hnumax*2 # By definition of tc2
		tc2=(H1/(Hnumax - B0) - 1.)**(1./p1) * 1e3/numax
	if tc2 >= tc1:
		tc2=tc1/2
	pos0=np.where(np.bitwise_and(freq >= numax-3*Dnu , freq <= numax+3*Dnu )) 
	Ampnumax=np.mean(spec_reg[pos0])/3.
	# -------------
	H1=H1-H2 # Correct H1 by removing the H2 contribution

	#init_param=[H1,tc1,p1,H2,tc2,p2,B0,Amp/2.,numax, 3*Dnu] ; WAS WORKING FOR MS. BUT WANT TO TRY SOMETHING ELSE LEARNT FROM RGB ===> WILL NEED TESTING ON MS
	init_param=[H1,tc1,p1,H2,tc2,p2,B0,(Amp/2+Ampnumax)/2.,numax, 3*Dnu]

	fin=np.sum(np.where(np.isfinite(init_param) == 0))
	if fin != 0:
		print('We detected some invalid values in the power spectrum (NaN or Inf)!')
		if (np.isfinite(tc2) == 0) or (np.isfinite(H2) == 0) or (np.isfinite(p2) == 0): # IF THE INF COMES FROM THE SECOND HARVEY (MOST LIKELY) THEN use numax to deternp.mine parameters
			p2=2.5
			pos0=np.where(np.bitwise_and(freq >= numax-3*Dnu , freq <= numax+3*Dnu) ) # Use the expected Dnu to get the average power range
			Hnumax=mean(spec_reg[pos0])
			H2=Hnumax*2 # By definition of tc2
			tc2=(H1/(Hnumax - B0) - 1.)**(1./p1) * 1e3/numax
			H1=H1-H2 # Correct H1 by removing the H2 contribution
			init_param=[H1,tc1,p1,H2,tc2,p2,B0,(Amp/2+Ampnumax)/2.,numax, 3*Dnu]
		else:
			print('Cannot proceed. Emmergency exit')
			#exit()
			err=True
			return [], err

	if (H1-H2)/H2 < 0.3: 
		H2=H2/4 # Force to lower H2 if this happens to be too close than H1

	if np.min(init_param) < 0:
		print('init param contains negative terms. Need debug')
		#print('The program will exit now')
		#exit()
		err=True
	return init_param, err

def do_data_file(freq, spec_reg, fileout, rebin=1):
	'''
		Generate a .data file compatible with the MCMC code
	'''
	err=False
	header="# File auto-generated by init_fit.py\n"
	header=header+"# freq_range=[ {:16.8f} , {:16.8f} ]\n".format(np.min(freq), np.max(freq))
	header=header+"# rebin =" + str(rebin) + "\n"
	header=header+"!            frequency               power\n"
	header=header+"*            (microHz)     (ppm^2/microHz)\n"

	if rebin > 1: # Averaging using a moving window
		x=[]
		y=[]
		pos0=0
		do_exit=False
		while pos0 < len(freq)  and do_exit == False:
			#pos0=i*rebin
			pos1=pos0+ rebin
			if pos1 >= len(freq):
				pos1=len(freq)-1
				do_exit=True
			x.append(np.mean(freq[pos0:pos1]) )
			y.append(np.mean(spec_reg[pos0:pos1]))
			pos0=pos1
		x=np.array(x)
		y=np.array(y)
	else:
		x=freq
		y=spec_reg
	try:
		f=open(fileout, 'w')
	except:
		print("Could not open the file :", fileout)
		print("Check that the path exists")
		err=True
	if err == False:
		try:
			f.write(header)
			for i in range(len(x)):
				line="{:20.8}  {:20.8}".format(x[i], y[i])
				f.write(line + '\n')
		except:
			print("Unknown error while writing on file: ", fileout)
			err=True
	return err

def do_model_file(init_param, relax_param, name_param, prior_param, name_prior, fileout, header="# File auto-generated by init_fit.py\n"):
	'''
		Generate a .model file compatible with the MCMC code
		More specifically, designed to be read by config.cpp:Config::read_inputs_prior_Simple_Matrix()
		init_param: A vector of initial guesses. Typically [H1, tc1, p1, H2, tc2, p2, B0, Anp.max, numax, sigma]
		relax_param : A series of 0 and 1, with 1 implying that the parameter is a variable and 0 that it is a constant
		prior_name: Name of the prior
		prior_param: A matrix with dimension (len(init_param), 5) that contains the following for each parameters:
					[val1, val2, val3, val4]
					Note that if the prior has less than 4 parameters, extra fields should be set to -9999
					as per the standards in the MCMC code
	'''
	err=False
	try:
		f=open(fileout, 'w')
	except:
		print("Could not open the file :", fileout)
		print("Check that the path exists")
		err=True
	if err == False:
		try:
			f.write(header)
			f.write("!")
			for name in name_param:
				line=" {:16s}".format(name)
				f.write(line)
			f.write("\n")			
			for param in init_param:
				line="{:16.6f}".format(param)
				f.write(line)
			f.write("\n")
			for relax in relax_param:
				line="   {:14d}".format(relax)
				f.write(line)
			f.write("\n")
			f.write("!")
			for name in name_prior:
				line=" {:16s}".format(name)
				f.write(line)
			f.write("\n")			
			for j in range(len(prior_param[0,:])):
				line=""
				for i in range(len(init_param)):
					line=line + "{:16.6f}".format(prior_param[i,j])
				f.write(line + "\n")
				#print(line)
		except:
			print("Unknown error while writing on file: ", fileout)
			err=True
	return err


def mode_initial_setup(spec_file, KIC_number, numax_guess, numax_uncertainty_guess, Amax_guess, outdir, fmin=0, fmax=5000, rebin=1):
	'''
		This function has for main role to guess what should be the initial
		parameters for the fit of the noise background
		spec_file: a npz or sav file that contain the spectrum
		ID_number: KIC or TIC number
		numax_guess: a guess of numax
		Anp.max_guess: a guess of Anp.max
		numax_uncertainty_guess: a guess on the uncertainty for numax
		outdir: Directory where the .model and .data files are saved
		fmin / fmax (optional): The np.minimum frequency to consider
		rebin (optional): If set to a value >1, it will make a spectrum using averaged point with a window of size of rebin/
						  BEWARE: This changes the statistics and thus errors in fitting. e.g. rebin = 2 ==>> p=2
	'''
	
	d=readsav(spec_file)
	freq=d['freq']
	spec_reg=d['spec_reg']
	freq=np.reshape(freq, freq.size)
	spec_reg=np.reshape(spec_reg, spec_reg.size)

	init_param,err=inital_param_S0(numax_guess,Amax_guess,freq, spec_reg, fmin, fmax)

	if err == False:
		# --- Plots ---
		H1=init_param[0]/(1. + (init_param[1]*freq*1e-3)**init_param[2])
		H2=init_param[3]/(1. + (init_param[4]*freq*1e-3)**init_param[5])
		B0=init_param[6]

		resol=freq[2] - freq[1]
		scoef=numax_guess/100./resol
		smoothspec=gaussian_filter1d(spec_reg, scoef, mode='reflect', truncate=4) 
		plt.plot(freq, smoothspec, color='Grey')
		plt.plot(freq, H1, color='Red', linestyle='--')
		plt.plot(freq, H2, color='Blue', linestyle='--')
		plt.plot(freq, np.repeat(B0, len(freq)), color='Black', linestyle='--')
		plt.plot(freq, H1+H2+B0, color='Black')
		ymin=np.min(smoothspec)*0.8
		ymax=np.max(smoothspec)*1.2
		plt.yscale('log')
		plt.xscale('log')
		plt.xlim(0.1, np.max(freq))
		plt.ylim(ymin, ymax)
		plt.savefig(outdir + "/" + KIC_number +  "_GaussfitGuess.png", dpi=300)
		plt.close()

		print('------ Initial guesses -----')
		print('H1 = ', init_param[0])
		print('tc1 = ', init_param[1])
		print('p1 = ', init_param[2])
		print('H2 = ', init_param[3])
		print('tc2 = ', init_param[4])
		print('p2 = ', init_param[5])
		print('B0 = ', init_param[6])
		print('Anp.max =', init_param[7])
		print('numax =', init_param[8])
		print('sigma =', init_param[9])
		print('---------------------------')

		# -- Handling priors and intial guesses in a suitable way for the MCMC process --
		# Note that the priors were basically inpsired from the IDL code 'Preset-Analysis-v8.1'

		#The parameter sigma should be proportional to nu_np.max... # and is GUG(xnp.min, xnp.max, sig_np.min, sig_np.max)
		if init_param[8] <= 250: 
			psigma=[init_param[8]/7.,init_param[8]/2,init_param[8]/20,init_param[9]]
		if init_param[8] > 250 and init_param[8] 	<= 400: 
			psigma=[20.,100.,5.,100.]
		if init_param[8] > 400 and init_param[8] <= 700:
			psigma=[25.,120.,30.,100]	
		if init_param[8] > 700 and init_param[8] <= 1200:
			psigma=[40.,200.,30.,100.]
		if init_param[8] > 1200 and init_param[8] <= 2700:
			psigma=[90.,250.,30.,100.]
		if init_param[8] > 2700:
			psigma=[150.,300.,30.,110.]

		#std_at_numax=np.std(spec_reg[np.where(np.bitwise_and(freq > init_param[8] - init_param[9], freq <= init_param[8] + init_param[9]))])
		prior_param=np.zeros((len(init_param),4)) - 9999
		name_param=["H1",                        "tc1"       , "p1"         , "H2"             , "tc2"    , "p2"     ,     "B0"               , "Amax"         ,   "numax"                                 , "sigma"   ]
		prior_name=["Jeffreys",                 "Uniform"    , "Fix"        , "Jeffreys"       , "Uniform", "Uniform",  "Uniform"             , "Jeffreys"      ,  "Uniform"                                , "GUG"     ]
		relax_param=[  1      ,                     1        ,   0          ,    1             ,    1     ,     1    ,      1                 ,      1          ,       1                                   ,    1      ]
		prior_param[:,0]=[np.std(spec_reg)/5   ,  5        , init_param[2]  ,   B0             ,   0      ,    0.5   ,      0                 ,      B0         , numax_guess-numax_uncertainty_guess       , psigma[0] ]
		prior_param[:,1]=[np.max(spec_reg)*100 ,2e3/np.min(resol), -9999    , np.max(spec_reg)    ,   56  ,     5    , 10*np.mean(init_param[6]) , np.max(spec_reg) , numax_guess + numax_uncertainty_guess*0.5 , psigma[1] ]
		prior_param[:,2]=[-9999            , -9999         , -9999        , -9999            ,  -9999   ,  -9999   ,  -9999                 , -9999           , -9999                                     , psigma[2] ]
		prior_param[:,3]=[-9999            , -9999         , -9999        , -9999            ,  -9999   ,  -9999   ,  -9999                 , -9999           , -9999                                     , psigma[3] ]

		# **************
		err=do_data_file(freq, spec_reg, outdir + "/" + KIC_number+"_Gaussfit.data", rebin=rebin)
		err=do_model_file(init_param, relax_param, name_param, prior_param, prior_name, outdir + "/" + KIC_number + "_Gaussfit.model", 
   					      header="# File auto-generated by init_fit.py\n# Fit of Gaussian mode Envelope\n# ID:"+ str(KIC_number)+"\n")
		# **************

def fast_test():
	spec_file='/Volumes/home/import/usbdisk2tb-3.5inch/Level0-Inputs/ts_rgb/concatenated/TF_10528917.sav'
	numax_guess=76
	numax_uncertainty_guess=40
	Amax_guess=5000
	KIC_number='10528917'
	outdir='/Users/obenomar/tmp/'
	rebin=4
	mode_initial_setup(spec_file, KIC_number, numax_guess, numax_uncertainty_guess, Amax_guess, outdir, fmin=0, fmax=5000, rebin=rebin)
