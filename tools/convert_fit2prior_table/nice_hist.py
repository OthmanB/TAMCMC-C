import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def round_any(number, increment, offset):
    return np.ceil((number - offset) / increment ) * increment + offset;

def nice_hist(samples, stats, ax=None, intervals=[True,True], binning=30, color=['black', 'gray', 'darkgray'], alpha=None, rotate=False, max_norm=None, label=[], yscale=1.2, linewidth=1):
	'''
		histograms with confidence intervals shown withtray areas
		samples: The data samples
		stats: Vector with [-2s, -1s, med, 1s, 2s]
		intervals: Boolean Vector [1s, 2s]:
			If 1s is True: The 1sigma interval is shown
			If 2s is True: The 2sigma interval is shown
		binning: Control the binning of the pdf
		color: Choose the set of color to show the (1) line of the pdf, (2) the 2sigma confidence interval and (3) the 1sigma confidence interval
		rotate: If set to true, invert the x and y axis 
		max_norm: Must be 'None' or any Real value. If it is a real value, the pdf will be normalised such that 'max_norm' is the maximum value of the pdf 
		yscale: Default is 1.2. This is a multiplicative factor of the plot area. It defines manually the maximum y-value of the plot area. 
	'''
	if rotate == True:
		orientation='horizontal'
	else:
		orientation='vertical'
	if ax == None:
		fig_1d, ax = plt.subplots(1, figsize=(12, 6))
	if max_norm == None:
		yvals, xvals, patches=ax.hist(samples,linestyle='-', bins=binning,histtype='step', color=color[0], density=True, orientation=orientation, linewidth=linewidth)#, label=label)
	else:
		fig_1d, ax_dummy = plt.subplots(1, figsize=(12, 6))
		yvals, xvals = np.histogram(samples,bins=np.linspace(np.min(samples), np.max(samples), binning))
		yvals=yvals*max_norm/np.max(yvals)
		ax.plot(xvals[0:-1], yvals,linestyle='-', ds='steps-mid', color=color[0], linewidth=linewidth)#, label=label)	
    #
	ymax=max(yvals)*yscale
	#
	xvals=xvals[0:len(yvals)]
	pos_1sig_interval=np.where(np.bitwise_and(xvals >= stats[1], xvals<=stats[3]))
	pos_2sig_interval=np.where(np.bitwise_and(xvals >= stats[0], xvals<=stats[4]))
	if rotate == False:
		if intervals[1] == True:
				ax.fill_between(xvals[pos_2sig_interval],yvals[pos_2sig_interval], color=color[1], alpha=alpha, step='post', interpolate=True)
		if intervals[0] == True:
			ax.fill_between(xvals[pos_1sig_interval],yvals[pos_1sig_interval], color=color[2], alpha=alpha,  step='post', interpolate=True)
		f = interpolate.interp1d(xvals, yvals, kind='cubic')
		ax.plot([stats[2], stats[2]], [0, f(stats[2])], color=color[0], linestyle='--', linewidth=linewidth)
		ax.set_ylim(0, ymax)
	else:
		if intervals[1] == True:
				ax.fill_betweenx(xvals[pos_2sig_interval],yvals[pos_2sig_interval], color=color[1], alpha=alpha, step='post', interpolate=True)
		if intervals[0] == True:
			ax.fill_betweenx(xvals[pos_1sig_interval],yvals[pos_1sig_interval], color=color[2], alpha=alpha,  step='post', interpolate=True)
		f = interpolate.interp1d(xvals, yvals, kind='cubic')
		ax.plot([0, f(stats[2])],[stats[2], stats[2]], color=color[0], linestyle='--', linewidth=linewidth)
		ax.set_xlim(ymax, 0)
	return xvals, yvals

def nice_hist_matrix(samples, stats, labels, binning=30, posplots=None, 
		    file_out='plots_pdfs.jpg', extra_data=[], interpolation_2d=None):
	'''
		Function to plot a nice correlation map between parameters
		samples: A 2D array with each line being sample of a given parameter: of the size (dim, Nsamples)
		stats: A 2D array with each line being the statitical summary of a given parameter (confidence intervals for -2s, 1s, med, 1s, 2s): of the size (dim, 5)
		binning: Level of binning used for the plot
		posplot: Allows to control the order of appearance of the variables in the graphical matrix
	'''
	Nticks=5
	intervals=[True,True]
	colors=['black', 'gray', 'darkgray']
	shift=1  # used to correct from the fact that we do not need redondant plots...
	dim=len(samples[:,0])
	#Nsamples=len(samples[0,:])
	if dim <2:
		print('Error: You need at least two sets of samples to plot their correlation matrix')
		exit()
	if posplots == None:
		posplots=np.zeros((2,dim), dtype=int)
		posplots[0,:]=np.linspace(0, dim-1, dim)
		posplots[1,:]=np.linspace(0, dim-1, dim)

	fig_2d, axes = plt.subplots(dim, dim, figsize=(12, 6))#, subplot_kw={'xticks': [], 'yticks': []})

	# Preparing the samples in order to plots the extra overplotted pdfs
	# Expected data format is a tri-list, contain itself within a list (of pdfs that need to plotted): 
	# 		- index list pointing to the concerned variable: []
	#		- Extra samples for each extra data serie: np.array((dim_extra_samples, Nsamples))
	#		- Extra statistics for each extra data serie : np.array((dim_extra_samples, 5))
	if extra_data != []:
		#extra_colors=['red', 'tomato','lightsalmon']
		extra_colors=['red', 'red','lightsalmon']
		Nextra=1
		'''
		Nextra=len(extra_data) # Number of extra pdfs that we are going to add
		do_extra=[]
		dim_extra=[]
		Nsamples_extra=[]
		extra_samples=np.zeros((Nextra, dim, Nsamples_extra))
		extra_stats=np.zeros((Nextra, dim, len(extra_data[2][0,:])))
		for k in range(Nextra):
			do_extra.append(np.repeat(False, dim))
			Nsamples_extra.append(len(extra_data[k][1][0,:])) # Get the dimension of the extra series (the number of series)
			dim_extra.append(len(extra_data[k][0])) # Get the number of samples of the extra series
			for i in range(0,dim):
				for j in range(0, dim_extra):
					if extra_data[k][0][j] == i:
						extra_samples[k][i,:]=extra_data[k][1][j,:]
						extra_stats[k][i,:]=extra_data[k][2][j,:]
						do_extra[k][i]=True
		'''
		Nsamples_extra=len(extra_data[1][0,:]) # Get the dimension of the extra series (the number of series). All series must have the same number of samples
		extra_samples=np.zeros((dim, Nsamples_extra))
		extra_stats=np.zeros((dim, len(extra_data[2][0,:])))
		do_extra=np.repeat(False, dim)
		dim_extra=len(extra_data[0]) # Get the number of samples of the extra series
		for i in range(0,dim):
			for j in range(0, dim_extra):
				if extra_data[0][j] == i:
					extra_samples[i,:]=extra_data[1][j,:]
					extra_stats[i,:]=extra_data[2][j,:]
					do_extra[i]=True
	# --------------------------------------------
	# -------- Filling upper line of PDF ---------
	# -------------------------------------------- 
	for i in range(0,dim-1):
		xvals, yvals=nice_hist(samples[posplots[0,i],:], stats[posplots[0,i],:], ax=axes[0,i+1], intervals=intervals, binning=binning, color=colors, alpha=None)
		# Process extra data, if any:
		if extra_data != []:
			for k in range(Nextra):
				if do_extra[i] == True:
					#xvals_extra, yvals_extra=nice_hist(extra_samples[k][i,:], extra_stats[k][i,:], ax=axes[0,i+1], intervals=intervals, binning=binning, color=extra_colors, alpha=None, max_norm=max(yvals))
					xvals_extra, yvals_extra=nice_hist(extra_samples[i,:], extra_stats[i,:], ax=axes[0,i+1], intervals=intervals, binning=binning, color=extra_colors, alpha=None, max_norm=max(yvals))
		#
		xmin=np.floor(np.min(xvals))
		xmax=np.ceil(np.max(xvals))
		step=(xmax-xmin)/Nticks
		axes[0,i+1].set_xlabel(labels[posplots[0,i]], fontsize=14)#, fontsize=20)
		axes[0,i+1].xaxis.set_major_locator(MultipleLocator(round_any(step, 10, 0)))
		axes[0,i+1].xaxis.set_minor_locator(AutoMinorLocator(2))
		axes[0,i+1].tick_params(which='minor', length=2)
		axes[0,i+1].xaxis.set_ticks_position('top')
		axes[0,i+1].xaxis.set_label_position('top')
		axes[0,i+1].xaxis.set_tick_params(labeltop=True, labelbottom=False)
		axes[0,i+1].set_yticks([])
 	# ---------------------------------------------
	# -------- Filling left column of PDF ---------
	# ---------------------------------------------
	for j in range(dim-1, 0, -1):
		xvals, yvals=nice_hist(samples[posplots[1,j],:], stats[posplots[1,j],:], ax=axes[j,0], intervals=intervals, binning=binning, color=colors, alpha=None, rotate=True)
		# Process extra data, if any:
		if extra_data != []:
			for k in range(Nextra):
				if do_extra[j] == True:
					#xvals_extra, yvals_extra=nice_hist(extra_samples[k][j,:], extra_stats[k][j,:], ax=axes[j,0], intervals=intervals, binning=binning, color=extra_colors, alpha=None, rotate=True, max_norm=max(yvals))
					xvals_extra, yvals_extra=nice_hist(extra_samples[j,:], extra_stats[j,:], ax=axes[j,0], intervals=intervals, binning=binning, color=extra_colors, alpha=None, rotate=True, max_norm=max(yvals))
		#
		xmin=np.floor(np.min(xvals))
		xmax=np.ceil(np.max(xvals))
		step=(xmax-xmin)/Nticks
		axes[j,0].set_ylabel(labels[posplots[1,j]], fontsize=14)#, fontsize=20)
		axes[j,0].yaxis.set_major_locator(MultipleLocator(round_any(step, 10, 0)))
		axes[j,0].yaxis.set_minor_locator(AutoMinorLocator(2))
		axes[j,0].tick_params(which='minor', length=2)
		axes[j,0].yaxis.set_ticks_position('left')
		axes[j,0].yaxis.set_label_position('left')
		axes[j,0].yaxis.set_tick_params(labelleft=True, labelright=False)
		axes[j,0].set_xticks([])
	# -----------------------------------------------
	# ------ Filling lines/columns with 2D PDF ------
	# -----------------------------------------------	
	for i in range(0, dim-1):
		for j in range(0, dim):
			if i < j:
				axes[j,i+1].hist2d(samples[posplots[0,i],:], samples[posplots[1,j],:], bins=binning, cmap=plt.cm.jet)
				axes[j,i+1].set_xticks([])
				axes[j,i+1].set_yticks([])
				#axes[j,i+1].set_xlabel('('+labels[posplots[0,i]] + ',' + labels[posplots[1,j]]+')', fontsize=8)#, fontsize=20)
			else:
				if (j !=0 and i!=0):
					axes[j,i+1].axis('off')
	axes[0,0].axis('off')
	fig_2d.savefig(file_out, dpi=300)

