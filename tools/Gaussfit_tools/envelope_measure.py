import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.io import readsav
from init_fit import do_data_file, do_model_file, mode_initial_setup
import os

def read_pow_PSD(filein):
    # read a file with .pow extension. This is files with 2 columns and the first character is "#" for comments/header
    freq = []
    PSD = []
    with open(filein, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                values = line.split()
                freq.append(float(values[0]))
                PSD.append(float(values[1]))
    return np.asarray(freq), np.asarray(PSD)

def read_Rafa_PSD(filein):
    d=fits.open(filein)
    return d[0].data[:,0]*1e6, d[0].data[:,1]

def boxcar_smooth(input_array, window_size):
    input_array = np.asarray(input_array)
    Nout = len(input_array)
    out = np.zeros(Nout)
    # Calculate the cumulative sum of the input array
    cumsum = np.cumsum(input_array)
    # Calculate the rolling sum using the cumulative sum
    out[window_size-1:Nout] = cumsum[window_size-1:] - cumsum[:-window_size+1]
    # Divide the rolling sum by the window size to get the mean
    out[window_size-1:Nout] /= window_size
    # Copy the input array values for the edges
    out[:window_size-1] = input_array[:window_size-1]
    out[Nout-window_size+1:] = input_array[Nout-window_size+1:]
    return out


def gauss_cdf(mu, sigma, x=None):
    import numpy as np
    if x is None:
        resol = 20 * sigma / 199
        x = np.arange(200) * resol - 10 * sigma + mu
    else:
        resol = np.mean(np.gradient(x))  # derivative of x
    logcdf = -np.log(np.sqrt(2. * np.pi) * sigma) - 0.5 * (x - mu) ** 2 / sigma ** 2
    maxi = np.sum(np.exp(logcdf))
    cdf = np.zeros(len(x))
    xcdf = cdf.copy()
    for i in range(len(x)):
        cdf[i] = np.sum(np.exp(logcdf[0:i])) / maxi
        xcdf[i] = x[i] + resol/2
    out = np.zeros((2, len(x)))
    out[0, :] = xcdf
    out[1, :] = cdf * 100
    return out

def read_data_file(datafile, extend_to_0=False, extend_type='max'):
# This function reads a data file, as per defined for the TAMCMC program
	f=open(datafile, 'r')
	data=f.read()
	f.close()

	txt=data.split('\n')
	count=0
	out=False
	while out == False:
		#print(txt[count])
		if txt[count][0] !='#':
			if txt[count][0] =='!':
				label=txt[count]
				count=count+1
			if txt[count][0] == '*':
				units=txt[count]
				count=count+1
			if txt[count][0] !='!' and txt[count][0] !='*':
				out=True
		else:
			count=count+1
	nmodels=len(txt[count].split()) - 2
	#
	x=np.zeros(len(txt))
	y=np.zeros(len(txt))
	if nmodels > 0:
		m=np.zeros((nmodels, len(txt)))
	else:
		print('No model found in the file... skipping it and returning -1')
		m=-1
	i=0
	for t in txt[count:]:
		line=t.split()
		if len(line)>0:
			x[i]=line[0]
			y[i]=line[1]
			for j in range(nmodels):
				m[j,i]=line[2+j]
			i=i+1
	x=x[0:i-1]
	y=y[0:i-1]
	if extend_to_0 == True:
		resol=x[2]-x[1]
		Nextra=int(np.floor(x.min()/resol)-1)
		if Nextra > 0:
			print('Warning: Extend to 0 is True: We will add datapoint using psf to the low range...')
			xextra=np.linspace(0, x.min()-resol, Nextra)
			xnew=np.zeros(len(x) + len(xextra))
			ynew=np.zeros(len(xnew))
			xnew[0:Nextra]=xextra
			xnew[Nextra:]=x
			true_xmin=x.min()
			if extend_type == 'max':
				print('extension using the max all the data points ')
				ynew[0:Nextra]=y[0:100].max()
			if extend_type == 'mean':
				print('extension using the mean of the lower freq 1microHz (or the first 10 bin) data points ')
				pos=np.where(x<=x.min()+0.5)
				if len(pos) <= 10:
					#print(' Warning: The lower 0.5 microHz range has less than 10 data points... using the 10 lowest frequency data points...')
					ynew[0:Nextra]=np.mean(y[0:10])
					#ynew[0:Nextra]=np.max(y[0:10])
				else:
					ynew[0:Nextra]=np.mean(y[pos])
			ynew[Nextra:]=y
			x=xnew
			y=ynew
		else:
			true_xmin=x.min()
	else:
		true_xmin=x.min()
	if nmodels >0:
		m=m[:,0:i-1]
	return x,y, m, extend_to_0, true_xmin


def envelope_measure(freq, spec_reg, fold_range, numax_step, fileout, sfactor=6, nu_rebin=None):
    if sfactor is None:
        sfactor = 1

    Teff_sun = 5777.
    Dnu_sun = 135.
    numax_sun = 3100.
    a0 = 0.77 
    alpha = np.float64(0.267)
    beta = np.float64(0.76)
    aa = 4.06e-42
    po = np.float64(11)
    K = alpha**(8./3.) * numax_sun**2 * Teff_sun / Dnu_sun**(8./3.)

    M = 1.25
    Teff = 6000
    L = 1.5
    b0 = M**(0.5 - a0) * (Teff/Teff_sun)**(3 - 3.5*a0) / L**(0.75 - a0)

    min_dpx = 4

    Ndata = len(freq)

    if max(freq) < 10:
        freq = freq * 1e6 

    fold_range[0] = max(fold_range[0], min(freq))
    fold_range[1] = min(fold_range[1], max(freq))

    resol = freq[1] - freq[0]

    x0 = freq
    s0 = spec_reg

    if nu_rebin is not None:
        if nu_rebin > resol:
            Nrebin = int((max(x0) - min(x0))/nu_rebin)
            x = np.interp(np.linspace(0, 1, Nrebin), np.linspace(0, 1, len(x0)), x0)
            s = np.interp(np.linspace(0, 1, Nrebin), np.linspace(0, 1, len(s0)), s0)
            print('Previous resolution:', resol, ' for a total of ' + str(Ndata) + ' points')
            print('New resolution:', nu_rebin, ' for a total of ' + str(Nrebin) + ' points')
            resol = nu_rebin
        else:
            x = x0
            s = s0
            raise Exception('rebin of the spectrum requested, but the provided a rebin resolution \
                    smaller than the current spectrum resolution! ')
    else:
        x = x0
        s = s0

    print('Selecting data within the requested folding range ' + 'fmin={}'.format(fold_range[0]) + \
            ' fmax={}'.format(fold_range[1]) + ' ...')
    zone = np.where(np.bitwise_and(x >= fold_range[0] , x <= fold_range[1]))[0]
    x = x[zone]
    s = s[zone]
    Nstep = int((fold_range[1] - fold_range[0])/numax_step)
    numax = numax_step * np.arange(Nstep+1) + fold_range[0]
    Dnu = Dnu_sun * b0 * (numax/numax_sun)**a0
    coef = 3.
    Gamma0 = np.float64(aa) * K**po * M**(2.*po/3.) * numax**(-2.*po*(3. - 4.*beta)/3.)
    coefsmoothG = 5

    print('Creating a slice-averaged spectrum centered on each possible numax....')
    print('   1. Setting the number of datapoints/slice using the numax_step parameter...')
    pslice = (x >= (numax[0] - Dnu[0]*coef)) & (x <= (numax[0] + Dnu[0]*coef))

    px0_0 = np.min(np.where(pslice)) 
    px1_0 = np.max(np.where(pslice))
    dpx = float(px1_0 - px0_0)
    print('   Initial dpx set to ' + str(dpx))
    if dpx < min_dpx:
        raise Exception('dpx is too small. Please either increase numax_step or reduce nu_rebin (if used).\n The minimum dpx is currently set to:' + min_dpx)
     
    print('2. Computing...')
    x_smooth=[]
    s_smooth=[]
    s_smooth_max=[]
    s_smooth_stddev=[]
    px0=px0_0 
    px1=px1_0
    for cpt in range(int(Nstep+1)-1):
        x_tmp = x[px0:px1]
        s_tmp = s[px0:px1]
        ss2 = boxcar_smooth(s0, int(coefsmoothG*Gamma0[cpt]/resol)) #if len(s0) > coefsmoothG*Gamma0[cpt]/resol else s0
        ss2 = ss2[zone]
        s2_tmp = ss2[px0:px1]
        x_smooth.append(np.mean(x_tmp))
        s_smooth.append(np.median(s_tmp))
        s_smooth_max.append(np.max(s2_tmp))
        s_smooth_stddev.append(np.std(s_tmp))
        pslice = np.where(np.bitwise_and(x >= (numax[cpt+1] - Dnu[cpt+1]*coef) , x <= (numax[cpt+1] + Dnu[cpt+1]*coef)))
        px0 = np.min(pslice)
        px1 = np.max(pslice)
    
    print('3. Extracting the guess for the maximum position and its significance....')
    x_r = x_smooth
    s_r= boxcar_smooth(np.asarray(s_smooth_max)/np.asarray(s_smooth), sfactor) #if len(s_smooth_max/s_smooth) > sfactor else s_smooth_max/s_smooth
    s_r0 = s_r
    s_r = (s_r - min(s_r))
    x_r=np.asarray(x_r)
    s_r=np.asarray(s_r)
    numax_guess0 = x_r[np.where(s_r == max(s_r))[0]]
    numax_guess0 = numax_guess0[0]
    uncertainty0 = np.interp(numax_guess0,numax, Dnu) * coef 
    rzone = np.where(np.bitwise_and(x_r >= (numax_guess0 - 2*uncertainty0) , (x_r <= numax_guess0 + 2*uncertainty0)))
    peak_coef = np.polyfit(x_r[rzone], s_r[rzone], 3)
    yfit = np.polyval(peak_coef, x_r[rzone])
    xselect = x_r[rzone]
    numax_guess = xselect[np.where(yfit == max(yfit))[0]][0]
    uncertainty = np.interp(numax_guess,numax, Dnu) * coef

    Gamma0_guess = float(aa) * K**po * M**(2.*po/3.) * numax_guess**(-2.*po*(3.-4.*beta)/3.)  # expected with a numax_guess
    s4Amax = boxcar_smooth(s0, int(coefsmoothG*Gamma0_guess/resol))  # define the optimum smooth a around numax
    Amax_guess = max(s4Amax[np.where(np.bitwise_and(x0 >= (numax_guess - uncertainty) , (x0 <= numax_guess + uncertainty)))])

    coefnoise = 5
    exits = 0
    while coefnoise >= 1 and exits == 0:
        noise_pos = np.where(np.bitwise_or(x_r <= (numax_guess - uncertainty*coefnoise) , (x_r >= numax_guess + uncertainty*coefnoise)))[0]
        if noise_pos == []: 
            coefnoise = coefnoise / 1.2
        else:
            nel = len(noise_pos)
            if nel > 10: 
                noise_mean = np.mean(s_r0[noise_pos])
                noise_stddev = np.std(s_r0[noise_pos])  
                exits = 1
            else:
                print('Warning during calculation of the peak significance: Less than 10 points in order to calculate the noise level')
                if nel < 2: 
                    coefnoise = coefnoise / 1.2
                    noise_pos = np.where((x_r >= numax_guess - uncertainty*coefnoise) & (x_r <= numax_guess + uncertainty*coefnoise))
                    print('Warning during calculation of the peak significance: Only few points to calculate noise level!')
                    print('    !!!!!!   PROCEED WITH EXTREME CAUTION WHEN INTERPRETING THE PEAK SIGNIFICANCE   !!!!!!')
                    noise_mean = np.mean(s_r0[noise_pos])
                    noise_stddev = np.std(s_r0[noise_pos])  
                else:
                    noise_mean = np.mean(s_r0[noise_pos])
                    noise_stddev = np.std(s_r0[noise_pos]) 
                    exits = 1

    ms_r = np.mean(s_r)
    noise_cdf = gauss_cdf(noise_mean, noise_stddev)  
    lvl_1sigma = np.interp(68.2, noise_cdf[0,:], noise_cdf[1,:])  # determine the 1 sigma detection level

    sig_gauss_sigma = [1, 2, 3, 4, 5, 6, 10, 20]  # significance in sigma
    sig_gauss_tab = [68.2, 95.4, 99.7, 99.99, 99.9999, 100, 100, 100]  # significance in percent

    significance = np.interp(Amax_guess/(lvl_1sigma + min(s_r0)), sig_gauss_sigma, sig_gauss_tab)
    sigs_r = lvl_1sigma

    print('------------------------------------------------------------------------------')
    print('       numax_guess              = ' + str(numax_guess) + ' microHz')
    print('       uncertainty (slice size) = '+ str(uncertainty) + ' microHz')
    print('       Amax_guess               = ' + str(Amax_guess) + ' ppm^2/microHz')
    print('       Significance probability = '+ str(significance) + ' %')
    print('       1-sigma detection level  = '+ str(lvl_1sigma) + ' ppm^2/microHz')
    print('------------------------------------------------------------------------------')

    s_r = (s_r - min(s_r))/max(s_r)
    yfit = yfit/max(yfit)

    # plot
    fig1, ax1 = plt.subplots(1, figsize=(12, 6), num=1, clear=True)  # create figure and axes
    ax1.plot(x_r, s_r, linewidth=2)  # plot data
    ax1.set_xlabel('Frequency (microHz)')
    ax1.set_ylabel('No unit')
    ax1.set_title('max(s)/average(s) calculated per slice of PSF. Peak significance: ' + str(significance) + '%')
    ax1.set_ylim([0,1])
    ax1.vlines(numax_guess, 0, 1, linestyles='dashed', color='blue', linewidth=3)  # add vertical line
    ax1.text(numax_guess + (fold_range[1]-fold_range[0])/100, 0.9, f'$\\nu_{{max}}$={numax_guess:.0f}+/-{uncertainty:.0f}', color='blue', fontsize=12)
    ax1.hlines(ms_r, min(x_r), max(x_r), linestyles='dashed', linewidth=3)  # add horizontal line
    ax1.hlines(sigs_r/max(s_r0), min(x_r), max(x_r), linestyles='dashed', color='red', linewidth=3)  # add horizontal line
    ax1.plot(x_r[rzone], yfit, color='green', linewidth=3)  # overlay plot

    # Second plot: original spectra in search range
    fig2, ax2 = plt.subplots(1, figsize=(12, 6), num=2, clear=True)  # create figure and axes
    ax2.plot(x0, s4Amax, linewidth=2)  # plot data
    ax2.set_xlabel('Frequency (microHz)')
    ax2.set_ylabel('ppm^2/microHz')
    ax2.set_title('PSF at optimal smooth for modes enhancement: Local view')
    ax2.set_xlim([fold_range[0], fold_range[1]])
    ax2.set_ylim([0, Amax_guess*2])
    ax2.vlines(numax_guess, 0, Amax_guess*2, linestyles='dashed', color='blue', linewidth=3)  # add vertical line
    ax2.text(numax_guess + (fold_range[1]-fold_range[0])/100, Amax_guess*2 - 10.*(Amax_guess*2)/100, f'$\\nu_{{max}}$={numax_guess:.0f}+/-{uncertainty:.0f}', color='blue', fontsize=12)
    ax2.hlines(Amax_guess, numax_guess-uncertainty, numax_guess+uncertainty, linestyles='dashed', color='blue', linewidth=3)  # add horizontal line
    ax2.text(numax_guess + (fold_range[1]-fold_range[0])/100, Amax_guess + (Amax_guess*2)/100, f'$A_{{max}}$={Amax_guess:.2f}', color='blue', fontsize=12)

    # Third plot: The whole original spectra
    fig3, ax3 = plt.subplots(1, figsize=(12, 6), num=3, clear=True)  # create figure and axes
    ax3.plot(x0, s4Amax, linewidth=2)  # plot data
    ax3.set_xlabel('Frequency (microHz)')
    ax3.set_ylabel('ppm^2/microHz')
    ax3.set_title('PSF at optimal smooth for modes enhancement: Global view')
    ax3.set_xlim([0, max(x0)])
    ax3.set_ylim([0, Amax_guess*10])
    ax3.vlines(numax_guess, 0, Amax_guess*10, linestyles='dashed', color='blue', linewidth=3)  # add vertical line
    ax3.hlines(Amax_guess, numax_guess-uncertainty, numax_guess+uncertainty, linestyles='dashed', color='blue', linewidth=3)  # add horizontal line

    fig1.savefig(fileout + "_1.jpg")
    fig2.savefig(fileout + "_2.jpg")
    fig3.savefig(fileout + "_3.jpg")
    return numax_guess, uncertainty, Amax_guess, significance

def make_guess_Kallinger2014(freq, spec_reg, outdir, ID, Amax_guess, numax_guess, numax_uncertainty_guess, rebin=1, do_Kallinger_model=False, do_data=False):
    Fnyquist=np.max(freq)
    
    if Fnyquist >= 1000.:
        stype='SC'
    else:
        stype='LC'

    # ---- B0 ----
    # SHORT CADENCE
    if stype == 'SC' :
        pos_0=np.min( np.where( freq >= (Fnyquist - 500.) )) # we take the last 500microHz to extract B0
    # LONG CADENCE
    if stype == 'LC':
        pos_0=np.min( np.where( freq >= Fnyquist - 10.)) # we take the last 10microHz to extract B0
    B0=np.mean(spec_reg[pos_0:])
    # ------------

    #The parameter sigma should be proportional to nu_np.max... # and is GUG(xnp.min, xnp.max, sig_np.min, sig_np.max)
    if numax_guess <= 250: 
        psigma=[numax_guess/7.,numax_guess/2,numax_guess/20,numax_uncertainty_guess]
    if numax_guess > 250 and numax_guess 	<= 400: 
        psigma=[20.,100.,5.,100.]
    if numax_guess > 400 and numax_guess <= 700:
        psigma=[25.,120.,30.,100]	
    if numax_guess > 700 and numax_guess <= 1200:
        psigma=[40.,200.,30.,100.]
    if numax_guess > 1200 and numax_guess <= 2700:
        psigma=[90.,250.,30.,100.]
    if numax_guess > 2700:
        psigma=[150.,300.,30.,110.]
    # ------ Values as per defined by Kallinger+2014, Table 2 ------
    # ---> a0
    k_Agran0=[3335, 9]
    s_Agran0=[-0.564, 0.002]
    k_taugran0=[836, 4]
    s_taugran0=[-0.886, 0.002]
    k_Agran=[k_Agran0[0], "k_Agran", "Gaussian", 1, k_Agran0[0], 10*k_Agran0[1], -9999, -9999]
    s_Agran=[s_Agran0[0], "s_Agran", "Gaussian", 1, s_Agran0[0], 10*s_Agran0[1], -9999, -9999]
    k_taugran=[k_taugran0[0],"k_taugran", "Gaussian", 1,  k_taugran0[0], 10*k_taugran0[1], -9999, -9999]
    s_taugran=[s_taugran0[0],"s_taugran", "Gaussian", 1,  s_taugran0[0], 10*s_taugran0[1], -9999, -9999]
    c0=[2, "c0", "GUG", 1, 2, 4, 0.1, 0.1]
    ''' We will not use this as we do not need a precise description of the very low frequencies for the mode fitting. The global parameters in Gaussian are used instead.
    a0_val=k_Agran*  numax_guess**s_Agran
    deriv_a0_k, deriv_a0_s, deriv_a0_numax=derivatives_power_law(k_Agran, s_Agran, numax_guess)
    err_a0=np.sqrt((deriv_a0_k*k_Agran[1])**2 + (deriv_a0_s*s_Agran[1])**2 + (deriv_a0_numax*numax_uncertainty_guess)**2)
    a0=[a0_val[0], "a0", "Gaussian", 1, a0_val[0], 5*err_a0, -9999, -9999]
    '''
    #
    # ---> a1
    k_a1=[3382,9]
    s_a1=[-0.609, 0.002]
    a1_val=k_a1[0]*  numax_guess**s_a1[0]
    deriv_a1_k, deriv_a1_s, deriv_a1_numax=derivatives_power_law(k_a1[0], s_a1[0], numax_guess)
    err_a1=np.sqrt((deriv_a1_k*k_a1[1])**2 + (deriv_a1_s*s_a1[1])**2 + (deriv_a1_numax*numax_uncertainty_guess)**2)
    a1=[a1_val, "a1", "Gaussian", 1, a1_val, 5*err_a1, -9999, -9999]
    #a1=[a1_val[0], "a1", "Uniform", 1, a1_val[0]/2, a1_val[0]*2, -9999, -9999]
    # ---> a2
    a2=[a1_val, "a2", "Gaussian", 1, a1_val, 5*err_a1, -9999, -9999]
    # ---> b1 :
    k1=[0.317, "k1", "Gaussian", 1, 0.317, 0.002, -9999, -9999]
    s1=[0.970, "s1", "Gaussian", 1, 0.970, 0.002, -9999, -9999]
    c1=[2    , "c1", "GUG"     , 1, 2   ,  4    , 0.1  , 0.1  ]
    # ---> b2 :
    k2=[0.948, "k2", "Gaussian", 1, 0.948, 0.003, -9999, -9999]
    s2=[0.992, "s2", "Gaussian", 1, 0.992, 0.002, -9999, -9999]
    c2=[2    , "c2", "GUG"     , 1, 2   ,  4    , 0.1  , 0.1  ]
    # ---> N0 :
    N0=[B0    , "N0", "Uniform" ,1, 0   , 10*B0 , -9999, -9999]
    # ---> Amax :
    Amax=[Amax_guess/10, "Amax", "Jeffreys", 1, B0/10, np.max(spec_reg), -9999, -9999]
    # ---> numax :
    numax_min=numax_guess-2*numax_uncertainty_guess
    numax_max=np.max(freq)
    numax=[numax_guess, "numax", "Uniform", 1, numax_min, numax_max, -9999, -9999]
    # ---> Gauss_sigma :
    Gauss_sigma=[(psigma[0] + psigma[1])/2, "Gauss_sigma", "GUG", 1, psigma[0], psigma[1], psigma[2], psigma[3]]
    # ---> Mass : # This parameter has been show to make no sense in the context of individual fitting
    #Mass=[1.2, "Mass", "Uniform", 1, 0.6, 4.0, -9999, -9999]
    # ---> mu_numax :
    mu_numax=[0.01, "mu_numax", "Uniform_abs", 1, 0.00, 0.05*numax_guess, -9999, -9999]
    # ---> omega_numax :
    omega_numax=[0.0025*numax_guess/5, "omega_numax", "Uniform", 1, 0.0, 0.05*numax_guess/5, -9999, -9999]
    # -------- Regrouping --------
    init_param= [k_Agran[0]    , s_Agran[0]  , k_taugran[0], s_taugran[0], c0[0], a1[0], a2[0], k1[0], s1[0], c1[0], k2[0], s2[0], c2[0], N0[0], Amax[0], numax[0], Gauss_sigma[0], mu_numax[0], omega_numax[0]]
    name_param= [k_Agran[1]    , s_Agran[1]  , k_taugran[1], s_taugran[1], c0[1], a1[1], a2[1], k1[1], s1[1], c1[1], k2[1], s2[1], c2[1], N0[1], Amax[1], numax[1], Gauss_sigma[1], mu_numax[1], omega_numax[1]]
    prior_name= [k_Agran[2]    , s_Agran[2]  , k_taugran[2], s_taugran[2], c0[2], a1[2], a2[2], k1[2], s1[2], c1[2], k2[2], s2[2], c2[2], N0[2], Amax[2], numax[2], Gauss_sigma[2], mu_numax[2], omega_numax[2]]
    relax_param=[k_Agran[3]    , s_Agran[3]  , k_taugran[3], s_taugran[3], c0[3], a1[3], a2[3], k1[3], s1[3], c1[3], k2[3], s2[3], c2[3], N0[3], Amax[3], numax[3], Gauss_sigma[3], mu_numax[3], omega_numax[3]]
    Nparams=len(name_param)
    prior_param=np.zeros((Nparams,4)) - 9999
    for j in range(4):
         prior_param[:,j]=[k_Agran[4+j]    , s_Agran[4+j]  , k_taugran[4+j], s_taugran[4+j], c0[4+j],  a1[4+j], a2[4+j], k1[4+j], s1[4+j], c1[4+j], k2[4+j], s2[4+j], c2[4+j], N0[4+j], Amax[4+j], numax[4+j], Gauss_sigma[4+j], mu_numax[4+j], omega_numax[4+j]]
    if do_data == True:
         err=do_data_file(freq, spec_reg, outdir + "/" + ID +"_KGaussfit.data", rebin=rebin)
    if do_Kallinger_model == True:
         hdr=      "# File auto-generated by envelope_measure.py\n"
         hdr=hdr + "# Fit of Gaussian mode Envelope with Kallinger2014 noise function\n"
         hdr=hdr + "# ID:"+ str(ID)+"\n# rebin="+str(rebin)+"\n"
         hdr=hdr + "# model_fullname= model_Kallinger2014_Gaussian\n" 
         err=do_model_file(init_param, relax_param, name_param, prior_param, prior_name, outdir + "/" + ID + "_KGaussfit.model", 
   			np.min(freq), np.max(freq), header=hdr)

def derivatives_power_law(k, s, numax):
     deriv_k=numax**s
     deriv_s=k*np.log(numax)* numax**s
     deriv_numax=k*s* numax**(s-1)
     return deriv_k, deriv_s, deriv_numax

def do_envelope(spec_file, outdir, do_Kallinger_model=False, do_Harvey_model=False, do_data=False, search_range=[500,4400], numax_step=10, rebin_resol=1):
    sfactor=6
    nu_rebin=None
    # get the extension of spec_file
    directory = os.path.dirname(spec_file)
    filename = os.path.basename(spec_file)
    ID = os.path.splitext(filename)[0]
    extension = spec_file.split('.')[-1]
    passed=False
    if extension == 'sav': # my legacy format
        d=readsav(spec_file)
        freq=d['freq']
        spec_reg=d['spec_reg']
        passed=True
    if extension == 'data':
        freq, spec_reg, models, extended, true_xmin=read_data_file(spec_file, extend_to_0=True)
        allowed_fmin=5
        if np.min(freq) > allowed_fmin:
             raise("Error: You must provide a full spectrum, not a section of it. The min(freq) should exceed 5")
        passed=True
    if extension == "fits": # Rafa concatenated lightcurves with in-painting
        try:
            freq, spec_reg=read_Rafa_PSD(spec_file)
            passed=True
        except:
             raise("Format not recognized. When using .fits, only Rafa PSD are accepted.")
    if extension == "pow": # KASOC concatenated lightcurves
         freq, spec_reg=read_pow_PSD(spec_file)
         passed=True
    if passed == False:
         raise("Format not recognized. Only .sav, .data .pow and .fits are accepted.")

    # Compute the level of rebining in bins, using the requested frequency resolution
    resol=freq[10] - freq[9]
    rebin=np.max([1,int(rebin_resol/resol)])
    fileout_jpg=outdir + "/" + ID + "_guess"
    numax_guess, uncertainty, Amax_guess, significance=envelope_measure(freq, spec_reg, search_range, 
                                                    numax_step, fileout_jpg, sfactor=sfactor, nu_rebin=nu_rebin)
    make_guess_Kallinger2014(freq, spec_reg, outdir, ID, Amax_guess, numax_guess, uncertainty, 
                             rebin=rebin, do_Kallinger_model=do_Kallinger_model, do_data=do_data)
    if do_Harvey_model == True:
         spec_file_4_harvey=outdir + "kplr004671239_kasoc-psd_llc_v1_KGaussfit.data"
         mode_initial_setup(spec_file, ID, numax_guess, numax_guess, Amax_guess, outdir, 
                            fmin=0, fmax=5000, rebin=1, do_S1=None, datatype='data',
                            do_data=False, do_model=True, extend_to_0=False)
    return ID, rebin


def do_all_envelopes(indir, outdir, extension='.data', do_Kallinger_model=True, do_data=True, 
                     rebin_resol=1, search_range=[500,4400], numax_step=10):
    # Check if outdir is a subdirectory of indir
    common_path = os.path.commonpath([indir, outdir])
    if common_path == indir:
        raise ValueError("outdir cannot be a subdirectory of indir")

    # Rest of the code...
    print(" ...Scanning for files with extension " + extension + "...")
    spec_files=[]
    for root, dirs, files in os.walk(indir):
        for f in files:
            if f.endswith(extension):
                spec_files.append(os.path.join(root, f))
    print(" ...Number of found files: ", len(spec_files))
    # Now do the envelope for each file
    print(" ...Processing found files...")
    index=1
    list_IDs=[]
    for spec_file in spec_files:
        print("[{}/{}] {}".format(index, len(spec_files), spec_file))
        ID, rebin=do_envelope(spec_file, outdir, do_Kallinger_model=do_Kallinger_model, do_data=do_data, 
                              rebin_resol=rebin_resol,search_range=search_range, numax_step=numax_step)
        index=index+1
        list_IDs.append("{:<80} {};".format(ID+"_KGaussfit",rebin))
    
    print("List to copy/past at the end of the config_presets.cfg")
    print("   table_ids= {} , 2; Number of lines and number of columns".format(len(list_IDs)))
    for l in list_IDs:
        print(l)

#outdir="/Users/obenomar/Work/dev/test_Kallinger2014/Data/inputs/Kallinger2014_Gaussian/LC_CORR_FILT_INP_ASCII_REBINNED/"
#indir="/Users/obenomar/Work/dev/test_Kallinger2014/Data/inputs/Kallinger2014_Gaussian/LC_CORR_FILT_INP_ASCII/"
print("  ----- Program to compute the noise + envelope parameters of a MCMC fit with TAMCMC ----")
indir = input("Enter the input directory (indir): ")
outdir = input("Enter the output directory (outdir): ")
do_all_envelopes(indir, outdir, extension='.data', do_Kallinger_model=True, do_data=False, rebin_resol=1, search_range=[250, 4400])