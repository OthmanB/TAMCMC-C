import numpy as np
import matplotlib.pyplot as plt
from  scipy.io import readsav
#from plotly.subplots import make_subplots
from matplotlib import gridspec
import a2_nl as rotfit # Most reading routine are in there
import os 

def get_dir_list(rootdir):
	'''
		A tiny function that scan a directory in order to find all of 
		its subdirectories
	'''
	dirs = [f.path for f in os.scandir(rootdir) if f.is_dir()]
	return dirs

# Converts nusini and nucosi fit parameters into a splitting a1 and an
# inclination i
def nusini_nucosi_disociator(data_dir, do_plots=True):

	#data_dir=rootdir + 'Files/'

	# Getting all of the indexes required for the process
	print("1.  Preparing data..")
	param_file='plength.txt'
	plength=rotfit.read_parameters_length(data_dir, param_file) # Read the parameters_length file and retrieves plength
	#Nf_el=plength[2:6] # Elements 2,3,4 and 5
	#i0_freq=sum(plength[0:2]) # Sum of elements 0 and 1 which are Nmax and lmax
	a1_test, Nsize=rotfit.read_sav(data_dir, sum(plength[0:6]))
	if  len(a1_test) != 1 : # If we deal with a model that fits directly a1 and inc
		print("No conversion required: The model that was used to fit those data fits directly a1 and inclination")
		print("We invite you to look at the file number: ", i0_a1, "  and ", i0_inc, " into the directory ", rootdir + 'Files/', ' to see the pdfs of a1 and of the inclination')
		print("The program will exit now")
		exit()
	else: # If we deal with a model that fits a1.sin(inc) and a1.cos(inc)
		i0_a1cosi=sum(plength[0:6]) +3
		i0_a1sini=sum(plength[0:6]) +4


	print("2.  Reading data..")
	a1cosi_samples, Nsize=rotfit.read_sav(data_dir, i0_a1cosi)
	a1sini_samples, Nsize=rotfit.read_sav(data_dir, i0_a1sini)
	print("3.  Converting data..")
	a1_samples=a1cosi_samples**2 + a1sini_samples**2 # Despite the variable name, what we have here is sqrt(a1.cosi) and sqrt(a1.sini)
	inc_samples=np.arctan(a1sini_samples**2/a1cosi_samples**2)*180./np.pi # Despite the variable name, what we have here is sqrt(a1.cosi) and sqrt(a1.sini)

	if do_plots == True:
		print("4. Plots...")
		plt.hist(a1_samples, bins=50)
		plt.xlabel('a1 (microHz)')
		plt.ylabel('PDF')
		plt.savefig(data_dir + 'a1.png')
		plt.hist(inc_samples, bins=50)
		plt.xlabel('inc (deg)')
		plt.ylabel('PDF')
		plt.savefig(data_dir + 'inc.png')
		plt.close('all')
	return a1_samples, inc_samples

def main(rootdir):
	data_dir=rootdir + '/Files/'
	a1_samples, inc_samples=nusini_nucosi_disociator(data_dir, do_plots=True)

if __name__ == "__main__":
	rootdir='/Users/obenomar/tmp2/fromdalma/RGB_depressed/data/products/'
	dirs=get_dir_list(rootdir)
	for d in dirs:
		print('Processing ', d, '...')
		main(d)
