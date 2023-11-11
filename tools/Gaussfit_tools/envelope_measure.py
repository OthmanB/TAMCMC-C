import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.io import readsav
from init_fit import do_data_file, do_model_file
import os

def read_Rafa_PSD(filein):
    d=fits.open(filein)
    return d[0].data[:,0]*1e6, d[0].data[:,1]


'''
def boxcar_smooth(input_array, window_size):
    input_array=np.asarray(input_array)
    Nout=len(input_array)
    out=np.zeros(Nout)
    for i in range(Nout):
        if (i >= int(window_size-1)/2) and (i<= Nout - int(window_size + 1)/2):
            for j in range(0, window_size-1):
                #print(i, j)
                out[i]=out[i] + input_array[i+j-int(window_size/2)]/window_size
        else:
            out[i]=input_array[i]
    return out
'''
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
		print('Warning: Extend to 0 is True: We will add datapoint using psf to the low range...')
		resol=x[2]-x[1]
		Nextra=int(np.floor(x.min()/resol)-1)
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
			#print('resol = ', resol)
			#print("Ndata_low_range=",len(pos[0]))
			if len(pos) <= 10:
				#print(' Warning: The lower 0.5 microHz range has less than 10 data points... using the 10 lowest frequency data points...')
				ynew[0:Nextra]=np.mean(y[0:10])
				#ynew[0:Nextra]=np.max(y[0:10])
			else:
				ynew[0:Nextra]=np.mean(y[pos])
		ynew[Nextra:]=y
		x=xnew
		y=ynew
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
        #if px1-px0 < 10:
        #    scoefmax = 2
        #else:
        #    scoefmax = 1.*(px1-px0)/50
        #if px1-px0 < 2:
        #    scoefmax = 1
        s_smooth_stddev.append(np.std(s_tmp))
        pslice = np.where(np.bitwise_and(x >= (numax[cpt+1] - Dnu[cpt+1]*coef) , x <= (numax[cpt+1] + Dnu[cpt+1]*coef)))
        px0 = np.min(pslice)
        px1 = np.max(pslice)
    
    # write on a debug file the quantity x_smooth and s_smooth and x_smooth_max in a 3 column format

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
    fig1, ax1 = plt.subplots()  # create figure and axes
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
    fig2, ax2 = plt.subplots()  # create figure and axes
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
    fig3, ax3 = plt.subplots()  # create figure and axes
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

    init_param=      [  3710.     ,   -0.613       ,   -0.26      ,     0.317      ,    0.970       ,        2         ,       0.948      ,       0.992          ,         2.        ,        B0        ,  Amax_guess/10   ,    numax_guess-numax_uncertainty_guess  ,  (psigma[0] + psigma[1])/2     ,     1.2          ,      0.01           ,         1]
    name_param=      [ "ka"        ,    "sa"       ,   "t"        ,   "k1"         ,    "s1"        ,     "c1"         ,      "k2"        ,       "s2"           ,      "c2"         ,      "N0"        ,    "Amax"        ,    "numax"                              ,   "Gauss_sigma"                ,    "Mass"        ,     "mu_numax"      ,    "omega_numax"]
    prior_name=      ["Gaussian"   ,   "Gaussian"  ,  "Gaussian"  ,   "Gaussian"   ,   "Gaussian"   ,     "Uniform"    ,     "Gaussian"   ,     "Gaussian"       ,     "Uniform"     ,     "Uniform"    ,    "Jeffreys"    ,    "Uniform"                            ,     "GUG"                      ,    "Uniform"     ,     "Uniform_abs"    ,    "Uniform" ]
    relax_param=     [  1         ,      1          ,       1      ,     1          ,    1           ,        1         ,        1         ,        1             ,         1         ,         1        ,       1          ,      1                                  ,       1                        ,     1            ,         1            ,        1 ]
    prior_param=np.zeros((len(name_param),4)) - 9999
    prior_param[:,0]=[  3710.     ,   -0.613       ,   -0.26      ,     0.317      ,    0.970       ,        2         ,       0.948      ,       0.992          ,         2.        ,        0.        ,     B0/10        ,    numax_guess-numax_uncertainty_guess  ,     psigma[0]                  ,     0.6          ,      0.01           ,         5]
    prior_param[:,1]=[   21       ,     0.002      ,   0.03       ,    0.002       ,    0.002       ,        5         ,       0.003      ,       0.002          ,         5         ,      10*B0       , np.max(spec_reg) ,    np.max(freq)                         ,     psigma[1]                  ,     4.0          ,        100          ,      20  ]
    prior_param[:,2]=[-9999.000000, -9999.000000   , -9999.000000 , -9999.000000   , -9999.000000   , -9999.000000     ,  -9999.000000    ,   -9999.000000       ,    -9999.000000   ,  -9999.000000    ,  -9999.000000    ,         -9999.0000                      ,     psigma[2]                  ,  -9999.000000    ,  -9999.00000        , -9999.000]
    prior_param[:,3]=[-9999.000000, -9999.000000   , -9999.000000 , -9999.000000   , -9999.000000   , -9999.000000     ,  -9999.000000    ,   -9999.000000       ,    -9999.000000   ,  -9999.000000    ,  -9999.000000    ,         -9999.0000                      ,     psigma[3]                  ,  -9999.000000    ,  -9999.00000        , -9999.000]
    if do_data == True:
         err=do_data_file(freq, spec_reg, outdir + "/" + ID +"_KGaussfit.data", rebin=rebin)
    if do_Kallinger_model == True:
         err=do_model_file(init_param, relax_param, name_param, prior_param, prior_name, outdir + "/" + ID + "_KGaussfit.model", 
   			np.min(freq), np.max(freq), header="# File auto-generated by envelope_measure.py\n# Fit of Gaussian mode Envelope with Kallinger2014 noise function\n# ID:"+ str(ID)+"\n# rebin="+str(rebin)+"\n")


def do_envelope(spec_file, outdir, do_Kallinger_model=False, do_data=False, search_range=[500,4400], numax_step=10, rebin=1):
    sfactor=6
    nu_rebin=None
    # get the extension of spec_file
    directory = os.path.dirname(spec_file)
    filename = os.path.basename(spec_file)
    ID = os.path.splitext(filename)[0]
    extension = spec_file.split('.')[-1]
    passed=False
    if extension == 'sav':
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
    if extension == "fits":
        try:
            freq, spec_reg=read_Rafa_PSD(spec_file)
            passed=True
        except:
             raise("Format not recognized. When using .fits, only Rafa PSD are accepted.")
    if passed == False:
         raise("Format not recognized. Only .sav, .data and .fits are accepted.")
    
    fileout_jpg=ID + "_guess"
    numax_guess, uncertainty, Amax_guess, significance=envelope_measure(freq, spec_reg, search_range, 
                                                    numax_step, fileout_jpg, sfactor=sfactor, nu_rebin=nu_rebin)
    make_guess_Kallinger2014(freq, spec_reg, outdir, ID, Amax_guess, numax_guess, uncertainty, 
                             rebin=rebin, do_Kallinger_model=do_Kallinger_model, do_data=do_data)

def test_envelope():
    outdir="/Users/obenomar/Work/dev/TAMCMC-C-v1.86.4/test/inputs/Kallinger2014_Gaussian/LC_CORR_FILT_INP_ASCII/rebinned/"
    spec_file="/Users/obenomar/Work/dev/TAMCMC-C-v1.86.4/test/inputs/Kallinger2014_Gaussian/LC_CORR_FILT_INP_ASCII/kplr008379927_91_COR_PSD_filt_inp.data"
    do_envelope(spec_file, outdir, do_Kallinger_model=True, do_data=True, rebin=10)
    print("Done")