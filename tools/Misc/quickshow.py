import os
import numpy as np
from read_outputs_tamcmc import bin2txt, getmodel_bin, read_datafile
import matplotlib.pyplot as plt


def quickshow(x,y,m, xm=None, c=['red', 'orange', 'blue', 'cyan', 'purple'], do_loglog=False, fileout=None):
	try:
		nmodels=len(m[:,0])
	except:
		nmodels=1
		if xm is None:
			mnew=np.zeros((nmodels, len(x)))
		else:
			mnew=np.zeros((nmodels, len(xm)))
		mnew[0,:]=m
		m=mnew
	fig, ax= plt.subplots(2, 1, sharex=True, num=1, clear=True)
	ax[0].plot(x, y, color='gray', label="Data")
	for j in range(nmodels):
		m_int=np.interp(x, xm, m[j, :])
		residuals = y / m_int
		ax[1].plot(x, residuals, color=c[j], label="Residuals " + str(j))
		if xm is None:
			ax[0].plot(x,m[j,:], color=c[j], label="Model " + str(j))
		else:
			ax[0].plot(xm,m[j,:], color=c[j], label="Model " + str(j))
	if do_loglog == True:
		ax[0].set_xscale('log')
		ax[0].set_yscale('log')
	ax[0].set_xlabel(r'Frequnecy $(\mu$Hz)')
	ax[0].set_ylabel(r'Power $(ppm/\mu$Hz)')
	ax[0].set_title('Quick Show Data vs Model')
	ax[0].legend()
	ax[1].plot(x, np.ones(len(x)), color='black', linestyle='--')
	ax[1].set_xlabel(r'Frequency $(\mu$Hz)')
	ax[1].set_ylabel('Residuals')
	plt.tight_layout()
	if fileout == None:
		plt.show()
	else:
		fig.savefig(fileout, dpi=300)

def do_show(dir_mcmc, process_name, model_name, fileout, cpp_path="../../bin/",
				phase="A", chain=0, first_index=10000, last_index=-1, period=1, 
				single_param_index=-1,erase_tmp=True,
				do_loglog=False):
	'''
		dir_mcmc: High-level Directory where the outputs of the TAMCMC are
		process_name: Name of the process that was analysed with TAMCMC. The program expects that the data for that specific process to be in
				dir_mcmc + process_name
		model_name: Name of the model that was used during the TAMCMC analysis
		cpp_path: directory where to find bin2txt and getmodel. By default, it is the place where binary files are created upon compilation
		phase: Either "B" for Burn-in files, "L" for Learning files or "A" for Acquire files
		chain: chain index. The coldest chain (target chain) is 0
		first_index: First sample index to be read
		last_index: Last index to be read. If set to -1, will read until the last index of the MCMC process
		period: specify the interval between to read index. One means that all samples will be read. Higher number allow faster read (and data reduction)
		single_param_index: If set to -1, all parameters will be read. Otherwise, only the parameter with the given index will be read
		erase_tmp: If set to True, will erase the temporary files created by the program
	'''
	current_dir=os.getcwd()
	cpp_version="1.85.0"
	outdir=current_dir + "/tmp/"
	file_data=dir_mcmc + process_name + "/inputs_backup/data_backup.zip"
	print(" ... Reading data file...")
	data, hdr=read_datafile(file_data)
	try:
		os.mkdir(outdir)
	except FileExistsError:
		pass
	print('... Gather posterior samples and extract plength...')
	samples, labels, isfixed, plength=bin2txt(dir_mcmc, process_name, phase=phase, chain=chain, 
					first_index=first_index, last_index=last_index, period=period, single_param_index=single_param_index,
	    			erase_tmp=erase_tmp, cpp_path=cpp_path, cpp_version=cpp_version, outdir=outdir, get_plength=True)
	print('... Computing the median for each parameter...')
	median=np.median(samples, axis=0)
	print("median: ", median)
	print("... Make the model...")
	resol=data[10,0]-data[9,0]
	ok=np.isnan(data[:,0])
	data=data[~ok,:]
	xr=[data[0,0], data[-1,0], resol]
	xm, m=getmodel_bin(model_name, median, plength, xr, cpp_path=cpp_path, outdir='tmp/', read_output_params=False, data_type='range')
	print("...Showing the plot...")
	quickshow(data[:,0],data[:,1],m, xm=xm, c=['red', 'orange', 'blue', 'cyan', 'purple'], do_loglog=do_loglog, fileout=fileout)

def do_all_show(dir_mcmc, phase, model_name, do_loglog=False):
	#dir_mcmc="/Users/obenomar/Work/dev/test_Kallinger2014/TAMCMC-C-v1.86.4/test/outputs/Kallinger2014_Gaussian/LC_CORR_FILT_INP_ASCII/"
	#phase="L"
	#model_name="model_Kallinger2014_Gaussian"
	#do_loglog=True
	# Get a List of all subdirectory within dir_mcmc
	list_subdir = [subdir for subdir in os.listdir(dir_mcmc) if not subdir.startswith('.')]
	# sort the list of subdir
	list_subdir.sort()
	# Loop over the list of subdir
	for process_name in list_subdir:
		print("    -----  Processing : ", process_name)
		fileout=dir_mcmc + "/" + process_name + "/" + process_name+"_bestfit.jpg"
		do_show(dir_mcmc, process_name, model_name, fileout, phase=phase, do_loglog=do_loglog)
		print("    -> File saved at : ", fileout)

print("  ----- Program to visualise outputs ----")
dir_mcmc = input("Enter the output directory for the tamcmc data: ")
phase = input("Enter the phase of the analysis (B/L/A): ")
model_name = input("Enter the model name for the analysis: ")

#dir_mcmc="/Users/obenomar/Work/dev/test_Kallinger2014/Data/outputs/Kallinger2014_Gaussian/LC_CORR_FILT_INP_ASCII_REBINNED/"
#phase="L"
#model_name="model_Kallinger2014_Gaussian"
do_loglog=True
do_all_show(dir_mcmc, phase, model_name, do_loglog=do_loglog)