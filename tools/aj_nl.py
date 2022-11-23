import os
import numpy as np
import matplotlib.pyplot as plt
from  scipy.io import readsav
#from plotly.subplots import make_subplots
from matplotlib import gridspec

def gen_data_filename(d, index):
	err=False
	if index < 10:
		s=d+'00' + str(index) + '.sav'
	if index >= 10 and index < 100:
		s=d+'0' + str(index) + '.sav'
	if index >=100 and index < 1000:
		s=d + str(index) + '.sav'
	if index >=1000:
		s='Cannot handle index numbers greater than 999'
		err=True
	return s, err

def read_sav(dir, ind):
	s, err=gen_data_filename(dir, ind)
	if err == True:
		print(s)
		exit()
	else:
		r=readsav(s)
		print(r['variable_name'])
	return r['param'], len(r['param'])

def read_parameters_length(dir, param_file):
	plength=[]
	with open(dir + param_file, 'r') as f:
		# Skip initial comments that starts with #
		for line in f:
				plength.append(line) # Each line contains a single value

	plength=np.array(plength, dtype=int)
	return plength

def read_freqs(d, i0, Nf_el):
	lmax=3
	r, Nsize=read_sav(d, i0)
	nu_l0=np.zeros((Nf_el[0], Nsize))
	nu_l1=np.zeros((Nf_el[1], Nsize))
	nu_l2=np.zeros((Nf_el[2], Nsize))
	nu_l3=np.zeros((Nf_el[3], Nsize))

	for en in range(Nf_el[0]):
		r, Nsize=read_sav(d, i0 + en)
		nu_l0[en, :]=r
	for en in range(Nf_el[1]):
		r, Nsize=read_sav(d, i0 + Nf_el[0]+ en)
		nu_l1[en, :]=r
	for en in range(Nf_el[2]):
		r, Nsize=read_sav(d, i0 + Nf_el[0] + Nf_el[1] + en)
		nu_l2[en, :]=r
	for en in range(Nf_el[3]):
		r, Nsize=read_sav(d, i0 + Nf_el[0] + Nf_el[1] + Nf_el[2] + en)
		nu_l3[en, :]=r
	return nu_l0, nu_l1, nu_l2, nu_l3


def read_rot(d, i0_a1, i0_a1cosi, i0_a1sini, i0_a2, i0_a3):
	if i0_a1 != -1:
		a1_param, Nsize=read_sav(d, i0_a1)
	else:
		a1cosi, Nsize=read_sav(d, i0_a1cosi)
		a1sini, Nsize=read_sav(d, i0_a1sini)
		a1_param=a1cosi**2 + a1sini**2 # Despite the variable name, what we have here is sqrt(a1.cosi) and sqrt(a1.sini)
	Na2=3
	a2_param=np.zeros((Na2, Nsize))
	for i in range(Na2):
		a2, Nsize=read_sav(d, i0_a2 + i)
		a2_param[i, :]=a2
	a3_param, Nsize=read_sav(d, i0_a3)
	return a1_param, a2_param, a3_param

def read_rot_aj(d, i0_a1, Nsize, Np=2, Nj=6):
	#Np=2 # two parameters per aj
	#Nj=6 # a1,a2,a3,a4,a5,a6
	#Nsize : Number of samples
	aj_param=np.zeros((Nj, Np, Nsize))
	j=0
	p=0
	for i in range(Np*Nj):
		aj_tmp, Nsize=read_sav(d, i0_a1 + i)
		aj_param[j,p, :]=aj_tmp
		p=p+1
		if p==Np:
			p=0
			j=j+1
	return aj_param

def read_inc(d, i0_inc, i0_a1cosi, i0_a1sini):
	if i0_inc !=-1:
		inc, Nsize=read_sav(d, i0_inc)
	else:
		a1cosi, Nsize=read_sav(d, i0_a1cosi)
		a1sini, Nsize=read_sav(d, i0_a1sini)
		inc=np.arctan(a1sini**2/a1cosi**2)*180./np.pi # Despite the variable name, what we have here is sqrt(a1.cosi) and sqrt(a1.sini)
	return inc

def write_aj_stats(fileout, l, nu, aj, append=False):
	Nj=len(aj[:,0,0])
	Nk=6
	header="# Table of aj coefficients in function of nu(n,l)\n# Col(0):l, Col(1):nu, Col(2+5(j-1))-Col(7 + 5*(j-1)):aj (for P(aj)=[2.25,16,50,84,97.75])\n" #a2 ~ 
	if append == False:
		f=open(fileout, 'w')
		f.write(header)
	else:
		f=open(fileout, "a")	
	for i in range(len(l)):
		string="{}   {}   ".format(l[i], nu[i]) 
		for j in range(Nj):
			for k in range(Nk):
				string=string + "{}   ".format(a[j,i,k])
		string=string+'\n'
		f.write(string)
	f.close()

def write_a2_stats(fileout, l, nu, a2, append=False):
	header="# Table of a2 coefficient in function of nu(n,l)\n# Col(0):l, Col(1):nu, Col(2)-Col(7):a2 (for P(a2)=[2.25,16,50,84,97.75])\n" #a2 ~ 
	if append == False:
		f=open(fileout, 'w')
		f.write(header)
	else:
		f=open(fileout, "a")	
	string="{}  {}  "
	for i in range(len(a2[0,:])):
		string=string + "{}  "
	string=string+"\n"
	for i in range(len(l)):
		f.write(string.format(l[i], nu[i], a2[i,1], a2[i,2], a2[i,3], a2[i,4], a2[i,5]))
	f.close()

def write_aj_raw_stats(fileout, nu_nl1, nu_nl2, nu_nl3, aj_mean_samples, Nf_el, rootdir):
	lmax=4
	f=open(fileout,"w")
	f.write("# Summary Table of RAW fit outputs for aj coefficients using a linear fit of the aj(nu,l) coefficients"+"\n")
	f.write("# Created using a2_nl.py : dostar_aj()"+" Version 27 Dec 2021\n")
	f.write("# Original content used to generate this file in : " + rootdir + "\n")
	f.write("! l = ")
	for el in range(1,lmax):
		for en in range(Nf_el[el]):
			f.write(" {0:1d}".format(el))
	f.write("\n")
	f.write("! nu_nl_obs = ")
	for nu in nu_nl1:
		f.write(" {0:10.6f}".format(nu))
	for nu in nu_nl2:
		f.write(" {0:10.6f}".format(nu))
	for nu in nu_nl3:
		f.write(" {0:10.6f}".format(nu))
	f.write("\n")
	f.write('# Mean and standard deviation of fitted parameters\n')
	f.write("# Col(0):a1, Col(1):a2, Col(2):a3, Col(3):a4, Col(4):a5, Col(5):a6, Col(6):err_a1, Col(7):err_a2, Col(8):err_a3, Col(8):err_a4, Col(8):err_a5, Col(8):err_a6\n")		
	for j in range(0,6):
		f.write(" {0:10.6f}".format(np.mean(aj_mean_samples[j,:])))
	for j in range(0,6):
		f.write(" {0:10.6f}".format(np.std(aj_mean_samples[j,:])))
	f.write("\n")
	f.close()
	
	
def write_aj_samples(fileout,l, nu, aj_samples, j=2): # If j is not specified, assumed that we deal with a2:
	header="# a{} Samples\n#l={}\n#nu={}\n".format(j,l,nu)
	f=open(fileout, 'w')
	f.write(header)
	for i in range(len(aj_samples)):
		f.write("{}\n".format(aj_samples[i]))
	f.close()

def get_files_list(rootdir):
	'''
		A tiny function that scan a directory in order to find all of 
		its files
	'''
	files = [f.path for f in os.scandir(rootdir) if f.is_file()]
	return files

def get_dirs_list(rootdir):
	'''
		A tiny function that scans a directory in order to find all of
		its subdirectories
	'''
	dirs= [d.path for d in os.scandir(rootdir) if d.is_dir()]
	return dirs
	
def get_files(rootdir, extension):

	f=get_files_list(rootdir)

	files=[]
	for ff in f:
		ext=ff.split('.')[-1]
		if ext == extension:
			files.append(ff)
	return files

def compute_a2nl_model_MS_Global_a1a2a3_HarveyLike(a2_terms, nu_l1, nu_l2, nu_l3, Nf_el):

	Nsize=len(a2_terms[0,:])
	if Nf_el[1] !=0:
		a2_l1=np.zeros((Nf_el[1], Nsize))
	else:
		a2_l1=np.zeros((1,1))
	if Nf_el[2] !=0:
		a2_l2=np.zeros((Nf_el[2], Nsize))
	else:
		a2_l2=np.zeros((1,1))
	if Nf_el[3] !=0:
		a2_l3=np.zeros((Nf_el[3], Nsize))
	else:
		a2_l3=np.zeros((1,1))
	
	# THIS IS THE IMPLEMENTATION OF model_MS_Global_a1a2a3_HarveyLike MODELS ONLY...
	for en in range(Nf_el[1]):
		a2_l1[en,:]=a2_terms[0,:] + a2_terms[1,:]*(1e-3*nu_l1[en]) + a2_terms[2,:]*(1e-6*nu_l1[en]**2); 
	for en in range(Nf_el[2]):
		a2_l2[en,:]=a2_terms[0,:] + a2_terms[1,:]*(1e-3*nu_l2[en]) + a2_terms[2,:]*(1e-6*nu_l2[en]**2); 
	for en in range(Nf_el[3]):
		a2_l3[en,:]=a2_terms[0,:] + a2_terms[1,:]*(1e-3*nu_l3[en]) + a2_terms[2,:]*(1e-6*nu_l3[en]**2); 

	return a2_l1, a2_l2, a2_l3

def compute_a2nl_model_MS_Global_aj(aj_terms, nu_l1, nu_l2, nu_l3, Nf_el):
	# aj_terms must be a Tensor of type (0:jmax-1, 0:Na-1, 0:Nsamples-1)
	# with jmax the max aj coefficient (jmax=6 usually) and Na, the max number of parameters (usually Na=2 for a linear fit)
	Nsize=len(aj_terms[0, 0,:])
	Np=len(aj_terms[0,:,0])
	lmax=3
	Nj=6
	if Nf_el[1] !=0:
		aj_l1=np.zeros((6,Nf_el[1], Nsize))
	else:
		aj_l1=np.zeros((6, 1,1))
	if Nf_el[2] !=0:
		aj_l2=np.zeros((6,Nf_el[2], Nsize))
	else:
		aj_l2=np.zeros((6,1,1))
	if Nf_el[3] !=0:
		aj_l3=np.zeros((6,Nf_el[3], Nsize))
	else:
		aj_l3=np.zeros((6,1,1))
	#
	# THIS IS THE IMPLEMENTATION OF model_MS_Global_aj MODELS ONLY...
	aj_mean=np.zeros((Nj,Nsize))
	Ntot=0
	for j in range(0,Nj):
		for el in range(1,lmax+1):
			for en in range(Nf_el[el]):
				tmp=np.zeros(Nsize)	
				#print('a{}0 = aj_terms[{},0,:] = {}'.format(j+1, j, aj_terms[j,0,:]))
				#print('a{}1 = aj_terms[{},1,:] = {}'.format(j+1, j, aj_terms[j,1,:]))
				if el == 1:
					aj_l1[j, en,:]=aj_terms[j,0,:] + aj_terms[j,1,:]*(1e-3*nu_l1[en]) 
					tmp=aj_l1[j,en,:]
					#print("l = {}      :       mean(a{}) = mean(aj_terms[{},0,:]) = {} ".format(el, j+1, j, np.mean(aj_terms[j,0,:])))
					Ntot=Ntot+1
				if el == 2:
					aj_l2[j, en,:]=aj_terms[j,0,:] + aj_terms[j,1,:]*(1e-3*nu_l2[en]) 
					tmp=aj_l2[j,en,:]
					#print("l = {}      :       mean(a{}) = mean(aj_terms[{},0,:]) = {} ".format(el, j+1, j, np.mean(aj_terms[j,0,:])))
					Ntot=Ntot+1
				if el == 3:
					aj_l3[j, en,:]=aj_terms[j,0,:] + aj_terms[j,1,:]*(1e-3*nu_l3[en])
					tmp=aj_l3[j,en,:]
					#print("l = {}      :       mean(a{}) = mean(aj_terms[{},0,:]) = {} ".format(el, j+1, j, np.mean(aj_terms[j,0,:])))
					Ntot=Ntot+1
				aj_mean[j,:]=aj_mean[j,:]+tmp
	#print("Ntot=", Ntot/Nj)
	#exit()
	return aj_l1, aj_l2, aj_l3, aj_mean/(Ntot/Nj)

def compute_confidence_intervals(l1_samples, l2_samples, l3_samples, Nf_el):
	conf_intervals=[2.25,16,50,84,97.75]
	if Nf_el[1] != 0:
		l1_stats=np.zeros((Nf_el[1], len(conf_intervals)))
	else:
		l1_stats=np.zeros((1,1))
	if Nf_el[2] != 0:
		l2_stats=np.zeros((Nf_el[2], len(conf_intervals)))
	else:
		l2_stats=np.zeros((1,1))
	if Nf_el[3] != 0:
		l3_stats=np.zeros((Nf_el[3], len(conf_intervals)))
	else:
		l3_stats=np.zeros((1,1))

	for en in range(Nf_el[1]):
		r=make_stats(l1_samples[en,:], confidence=conf_intervals) # Get the confidence intervals by making a cdf
		l1_stats[en,:]=r
	for en in range(Nf_el[2]):
		r=make_stats(l2_samples[en,:], confidence=conf_intervals) # Get the confidence intervals by making a cdf
		l2_stats[en,:]=r
	for en in range(Nf_el[3]):
		r=make_stats(l3_samples[en,:], confidence=conf_intervals) # Get the confidence intervals by making a cdf
		l3_stats[en,:]=r

	return l1_stats, l2_stats, l3_stats

def make_stats(samples, confidence=[2.25,16,50,84,97.75]):
	N=len(samples)
	s=np.sort(samples)
	cdf = 100.*np.array(range(N))/float(N) # in %
	r=np.interp(confidence, cdf, s)
	#plt.plot(s, cdf)
	#plt.plot([r[0], r[0]], [0, 100])
	#plt.plot([r[1], r[1]], [0, 100])
	#plt.plot([r[2], r[2]], [0, 100])
	#plt.plot([r[3], r[3]], [0, 100])
	#plt.show()
	return r

def make_error_from_stats(stats):
	try:
		err=np.zeros((2, len(stats[:,0]))) 
		err[0,:]=stats[:, 2] - stats[:, 1]
		err[1,:]=stats[:, 3] - stats[:, 2]
	except:
		try:
			err=np.zeros(2)
			err[0]= stats[2] - stats[1]
			err[1]= stats[3] - stats[2]
		except:
			print("Error : Invalid format of stats in make_error_from_stats()")
			print("        stats = ", stats)
			exit()
	return err

def dostar_aj(rootdir, fileout='Results'):	
	Np=2 # Number of parameter per aj
	Nj=6 # Number of aj terms
	# 
	# Getting all of the indexes required for the process
	print("1.  Preparing data..")
	param_file='plength.txt'
	plength=read_parameters_length(rootdir + 'Files/', param_file) # Read the parameters_length file and retrieves plength
	Nf_el=plength[2:6] # Elements 2,3,4 and 5
	i0_freq=sum(plength[0:2]) # Sum of elements 0 and 1 which are Nmax and lmax
	#print("plength: ", plength)
	a1_test, Nsize=read_sav(rootdir + 'Files/', sum(plength[0:6]))
	#print("len(a1_test) =", len(a1_test))
	#print("a1_test =", a1_test)	
	i0_a1=sum(plength[0:6]) # Position after Nf_el list
	i0_inc=sum(plength[0:-2]) # The last parameter is before the extra parameters that are at the end ([-1] position)
	
	# Get the Frequencies samples
	print("2. Gathering frequencies...")
	nu_l0_samples, nu_l1_samples, nu_l2_samples, nu_l3_samples=read_freqs(rootdir + 'Files/', i0_freq, Nf_el)
	# Get the rotation parameters in form of samples
	print("3. Gathering rotation parameters...")
	aj_param_samples=read_rot_aj(rootdir + 'Files/', i0_a1, Nsize, Np=Np, Nj=Nj)

	# Get the inclination in form of samples
	print("4. Gathering inclination parameters...")
	inc_samples=read_inc(rootdir + 'Files/', i0_inc, -1, -1)

	# Compute the a2_nl on the form of samples
	print("5. Compute aj_nl...")
	aj_l1_samples, aj_l2_samples, aj_l3_samples, aj_mean_samples=compute_a2nl_model_MS_Global_aj(aj_param_samples, nu_l1_samples, nu_l2_samples, nu_l3_samples, Nf_el)
	aj_l1_samples=aj_l1_samples*1e3 
	aj_l2_samples=aj_l2_samples*1e3 
	aj_l3_samples=aj_l3_samples*1e3 
	aj_mean_samples=aj_mean_samples*1e3
	# Extract basic statistics from the a2 samples and frequency samples
	print("6. Get stats using the samples of all of the parameters...")
	nu_l1_stats, nu_l2_stats, nu_l3_stats=compute_confidence_intervals(nu_l1_samples, nu_l2_samples, nu_l3_samples, Nf_el)
	aj_l1_stats=[]
	aj_l2_stats=[]
	aj_l3_stats=[]
	aj_mean_stats=[]
	for j in range(0,Nj):
		l1_stats, l2_stats, l3_stats=compute_confidence_intervals(aj_l1_samples[j,:,:], aj_l2_samples[j,:,:], aj_l3_samples[j,:,:], Nf_el)
		mean_stats=make_stats(aj_mean_samples[j,:])
		aj_l1_stats.append(l1_stats)
		aj_l2_stats.append(l2_stats)
		aj_l3_stats.append(l3_stats)
		aj_mean_stats.append(mean_stats)
	inc_stats=make_stats(inc_samples)
	#	
	#exit()
	# Generate pdfs for all the calculated quantities
	print("7. Plots...")
	print("     - aj In function of frequency...")
	fig, axs = plt.subplots(Nj)
	el=1
	for j in range(0,Nj):
		#print("j = ", j)
		#print(' aj_l{}_stats[{}] = {}'.format(el,j+1, aj_l1_stats[j]))
		axs[j].set_xlabel('Frequency (microHz)')
		axs[j].set_ylabel('a{}_nl (nHz)'.format(j+1))
		axs[j].plot([np.min(nu_l1_stats[:,2]),np.max(nu_l1_stats[:,2])], [0,0], linestyle='dashed')
		yerr=make_error_from_stats(aj_l1_stats[j][:,:])
		xerr=make_error_from_stats(nu_l1_stats)
		axs[j].errorbar(nu_l1_stats[:,2], aj_l1_stats[j][:, 2], xerr=xerr, yerr=yerr)
	plt.savefig(rootdir+'/'+fileout+'.png')
	print("     - average aj...")
	fig, axs = plt.subplots(Nj)
	for j in range(0,Nj):
		print('a{}_mean_stats = {}'.format(j+1, aj_mean_stats[j]))
		axs[j].hist(aj_mean_samples[j,:], bins=50)
		axs[j].set_xlabel('<a{}_nl> (nHz)'.format(j+1))
		axs[j].set_ylabel("PDF")
		#axs[j].plot([0,0], [,0], linestyle='dashed')
		xerr=make_error_from_stats(aj_mean_stats[j])
	plt.savefig(rootdir+'/'+fileout+'_pdf.png')
	#plt.show()
	#exit()
	print("7. Writing data on aj...")
	print("     - Table of aj raw parameters with gaussian representation...")
	write_aj_raw_stats(rootdir+"/aj_raw.txt", nu_l1_stats[:,2], nu_l2_stats[:,2], nu_l3_stats[:,2], aj_mean_samples, Nf_el, rootdir)
	#print("    - Table with all frequencies and individual samples")
	#for j in range(0,Nj):
	#	for l in range(1,3):
	#		if l == 1:
	#			aj=aj_l1_stats[j,:]  #aj_l1_samples[j,:,:]
	#			aj_s=aj_l1_samples[j,:,2]
	#			nu=nu_l1_stats[:,2] # Median values only
	#			append=False
	#		if l == 2:
	#			aj=aj_l2_stats[j,:]
	#			aj_s=aj_l2_samples[j,:,2]
	#			nu=nu_l2_stats[:,2] # Median values only
	#			append=True
	#		if l == 3:
	#			a2=a2_l3_stats
	#			a2_s=a2_l3_samples
	#			nu=nu_l3_stats[:,2] # Median values only
	#			append=True
	#		for i in range(Nf_el[l]):
	#			write_aj_stats(rootdir+"a2stats_n"+str(i)+".txt", np.repeat(l, len(nu)), nu, a2, append=append)
	#			write_aj_samples(rootdir+"Files/a2samples_n"+str(l)+"_"+str(i)+".txt", l, nu[i], a2_s[i,:])


def convert_ajraw2data(aj_file, outfile=None):
	'''
		A function that generates a .data file with the 
		indexes j as x-axis, aj coefficients as y-axis and 
		uncertainties as z-axis
		It also creates a .model file with the rest of the inputs
		that can be directly read using the matrix function from TAMCMC
	'''



def dostar_a1a2a3(rootdir, fileout='Results'):
# Going to be obselete and replaced by aj models

	#rootdir='/Users/obenomar/tmp/TRASH/a2-fits/products/1111/kplr003427720_kasoc-psd_slc_v1_1111/'
	#rootdir='/Users/obenomar/tmp/TRASH/a2-fits/products/1111/kplr008379927_kasoc-psd_slc_v2_1111/'
	
	# Getting all of the indexes required for the process
	print("1.  Preparing data..")
	param_file='plength.txt'
	plength=read_parameters_length(rootdir + 'Files/', param_file) # Read the parameters_length file and retrieves plength
	Nf_el=plength[2:6] # Elements 2,3,4 and 5
	i0_freq=sum(plength[0:2]) # Sum of elements 0 and 1 which are Nmax and lmax
	#print("plength: ", plength)
	a1_test, Nsize=read_sav(rootdir + 'Files/', sum(plength[0:6]))
	#print("len(a1_test) =", len(a1_test))
	#print("a1_test =", a1_test)	
	if  len(a1_test) != 1 : # If we deal with a model that fits directly a1 and inc
		i0_a1=sum(plength[0:6]) # Position after Nf_el list
		i0_inc=sum(plength[0:-2]) # The last parameter is before the extra parameters that are at the end ([-1] position)
		i0_a1cosi=-1
		i0_a1sini=-1
	else: # If we deal with a model that fits a1.sin(inc) and a1.cos(inc)
		i0_a1=-1
		i0_inc=-1
		i0_a1cosi=sum(plength[0:6]) +3
		i0_a1sini=sum(plength[0:6]) +4
	i0_a2=sum(plength[0:6]) + 6 # a2 parameters starts after the asymetry parameter (and is made of 3 parameters)
	i0_a3=sum(plength[0:6]) + 2 # a3 is after the asphericity parameter and two position further from a1  

	# Get the Frequencies samples
	print("2. Gathering frequencies...")
	nu_l0_samples, nu_l1_samples, nu_l2_samples, nu_l3_samples=read_freqs(rootdir + 'Files/', i0_freq, Nf_el)
	# Get the rotation parameters in form of samples
	print("3. Gathering rotation parameters...")
	a1_samples, a2_param_samples, a3_samples=read_rot(rootdir + 'Files/', i0_a1, i0_a1cosi, i0_a1sini, i0_a2, i0_a3)
	a3_samples=a3_samples*1e3 # Conversion in nHz

	# Get the inclination in form of samples
	print("4. Gathering inclination parameters...")
	inc_samples=read_inc(rootdir + 'Files/', i0_inc, i0_a1cosi, i0_a1sini)

	# Compute the a2_nl on the form of samples
	print("5. Compute a2_nl...")
	a2_l1_samples, a2_l2_samples, a2_l3_samples=compute_a2nl_model_MS_Global_a1a2a3_HarveyLike(a2_param_samples, nu_l1_samples, nu_l2_samples, nu_l3_samples, Nf_el)
	a2_l1_samples=a2_l1_samples*1e3 
	a2_l2_samples=a2_l2_samples*1e3 
	a2_l3_samples=a2_l3_samples*1e3 

	# Extract basic statistics from the a2 samples and frequency samples
	print("6. Get stats using the samples of all of the parameters...")
	a2_l1_stats, a2_l2_stats, a2_l3_stats=compute_confidence_intervals(a2_l1_samples, a2_l2_samples, a2_l3_samples, Nf_el)
	nu_l1_stats, nu_l2_stats, nu_l3_stats=compute_confidence_intervals(nu_l1_samples, nu_l2_samples, nu_l3_samples, Nf_el)
	a1_stats=make_stats(a1_samples)
	a3_stats=make_stats(a3_samples)
	inc_stats=make_stats(inc_samples)
	#print(a2_l1_stats*1e3)

	#exit()
	# Generate pdfs for all the calculated quantities
	print("6. Plots...")
	fig = plt.figure(constrained_layout=True)
	gs = fig.add_gridspec(2, 3)
	f_ax1 = fig.add_subplot(gs[0, :]) # This plot is taking the whole line (3 blocks, upper line)
	f_ax1.set_xlabel('Frequency (microHz)')
	f_ax1.set_ylabel('a2_nl (nHz)')
	f_ax1.plot([np.min(nu_l1_stats[:,2]),np.max(nu_l1_stats[:,2])], [0,0], linestyle='dashed')
	yerr=make_error_from_stats(a2_l1_stats)
	xerr=make_error_from_stats(nu_l1_stats)
	f_ax1.errorbar(nu_l1_stats[:,2], a2_l1_stats[:, 2], xerr=xerr, yerr=yerr)

	f_ax2 = fig.add_subplot(gs[1, 0]) # This plot is on the bottom left corner: [1, 0]
	f_ax2.set_xlabel('a1 (microHz)')
	f_ax2.set_ylabel('PDF')
	f_ax2.hist(a1_samples, bins=50)

	f_ax3 = fig.add_subplot(gs[1, 1]) # This plot is on the bottom one step right from the left: [1, 1]
	f_ax3.set_xlabel('inc (deg)')
	f_ax3.set_ylabel('PDF')
	f_ax3.hist(inc_samples, bins=50)

	f_ax4 = fig.add_subplot(gs[1, 2]) # This plot is on the bottom two step right from the left: [1, 2]
	f_ax4.set_xlabel('a3 (nHz)')
	f_ax4.set_ylabel('PDF')
	f_ax4.hist(a3_samples, bins=50)
	plt.savefig(fileout+'.png')
	#plt.show()

	print("7. Writing data on a2...")
	for l in range(1,3):
		if l == 1:
			a2=a2_l1_stats
			a2_s=a2_l1_samples
			nu=nu_l1_stats[:,2] # Median values only
			append=False
		if l == 2:
			a2=a2_l2_stats
			a2_s=a2_l2_samples			
			nu=nu_l2_stats[:,2] # Median values only
			append=True
		if l == 3:
			a2=a2_l3_stats
			a2_s=a2_l3_samples
			nu=nu_l3_stats[:,2] # Median values only
			append=True
		for i in range(Nf_el[l]):
			write_a2_stats(rootdir+"a2stats_n"+str(i)+".txt", np.repeat(l, len(nu)), nu, a2, append=append)
			write_aj_samples(rootdir+"Files/a2samples_n"+str(l)+"_"+str(i)+".txt", l, nu[i], a2_s[i,:])

def doallstars(rootdir):

	failed=[]
	dirs=get_dirs_list(rootdir)
	for d in dirs:
		try:
			dostar(d+'/', fileout=d+'/a1a2a3inc_sum.png')	
		except:
			failed.append(d)
	if len(failed) !=0:
		print("Directories that could not be processed due to a failure:")
		print(failed)
		
if __name__ == "__main__":
	print("Please use doallstars([dir]) to get a1,a2,a3 and inc summary for a group of stars")
	print("Or use dostar([dir]) for a single star analysis")
	
