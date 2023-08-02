import shutil
import os 
import sys
import matplotlib.pyplot as plt
import numpy as np
from read_outputs_tamcmc import bin2txt
from nice_hist import nice_hist,nice_hist_matrix

def make_prior_table(index1, binning, output_dir, dir_tamcmc_outputs, process_name, phase='A', chain=0, 
		     first_index=0, last_index=-1, period=1,
		     index2=None, cpp_prg="../bin/", convert_a1cosisini=False):
	'''
		Program that allows you to create a 1D or 2D table in a format suitable for the TAMCMC 
		tabulated priors.
		index1: Index of the a parameter that will be extracted and for which a prior table will be made. 
			Note that the index starts at 1.
		output_dir: Directory in which the table and plot (used for checks) of the table will be created
		index2: If provided, will be used to create a 2D table, by using jointly index1 parameter values
			Note that the index starts at 1.
		output_dir : Output directory where to save the data
		binning: The number of bins in the histogram. For a 1D prior, must be an integer. For a 2D prior it
			can be either an integer (in that case Nx_bin=Ny_bin) or an array with two values [Nx_bin, Ny_bin]
		cpp_prg: The directory in which the 'bin2txt' program is
		convert_a1cosisini: If set to True, when the index1 and index2 correspond to parameters of names :
					- sqrt(splitting_a1).cosi    and   sqrt(splitting_a1).sini
	'''
	# Define and create a temporary directory for the extracted binary data
	outdir_tmp=output_dir + '/tmp/'
	try:
		if not os.path.exists(outdir_tmp):
			os.makedirs(outdir_tmp)
		smcmc1, label1, isfixed1=bin2txt(dir_tamcmc_outputs, process_name, phase=phase, chain=chain, 
				first_index=first_index, last_index=last_index, period=period, single_param_index=index1,
				erase_tmp=True, cpp_path=cpp_prg, cpp_version="1.85.0", outdir=outdir_tmp, get_plength=False)	
		stats1=[np.min(smcmc1), 
	  			np.median(smcmc1) - np.std(smcmc1), 
				np.median(smcmc1), 
				np.median(smcmc1) + np.std(smcmc1), np.max(smcmc1),
				np.max(smcmc1)]
		if isfixed1 == True:
			print("Error: The requested index1 is a fixed parameter. This is not suitable for building a prior", file=sys.stderr)
			print("      index1        : ", index1, file=sys.stderr)
			print("      parameter name: ", label1[0], file=sys.stderr)
			raise ValueError("Exiting the program...")
		
		if index2 != None:
			smcmc2, label2, isfixed2=bin2txt(dir_tamcmc_outputs, process_name, phase=phase, chain=chain, 
				first_index=first_index, last_index=last_index, period=period, single_param_index=index2,
				erase_tmp=True, cpp_path=cpp_prg, cpp_version="1.85.0", outdir=outdir_tmp, get_plength=False)	
			stats2=[np.min(smcmc2), 
					np.median(smcmc2) - np.std(smcmc2), 
					np.median(smcmc2), 
					np.median(smcmc2) + np.std(smcmc2), np.max(smcmc2),
					np.max(smcmc2)]		
			if isfixed2 == True:
				print("Error: The requested index2 is a fixed parameter. This is not suitable for building a prior", file=sys.stderr)
				print("      index2        : ", index2, file=sys.stderr)
				print("      parameter name: ", label2[0], file=sys.stderr)
				raise ValueError("Exiting the program...")
		
		# Case of a 1D table
		if index2 == None:
			# Definitions
			file_out=output_dir + "/" + process_name + "_0.priors"
			file_out_plot=file_out + ".jpg"
			# Conversion of sqrt(splitting_a1).cosi, qrt(splitting_a1).sini is impossible without index1 and index2 
			if convert_a1cosisini == True and (label1 == "sqrt(splitting_a1).cosi" or label1 == "sqrt(splitting_a1).sini"):
				print("Warning: convert_a1cosisini = True but only index1 was provided")
				print("         Converting sqrt(splitting_a1).cosi, qrt(splitting_a1).sini would require index1 and index2 to be set")
				print("         Ignoring this parameter and making the 1D table as it is")
			# Make the plot of the table
			fig, ax = plt.subplots(1, figsize=(12, 6))
			xvals, yvals=nice_hist(smcmc1, stats1, ax=ax, intervals=[True,True], 
				binning=binning, color=['black', 'gray', 'darkgray'], alpha=None, rotate=False, 
				max_norm=None, label=[], yscale=1.2, linewidth=1)
			fig.savefig(file_out_plot, dpi=300)
			# Make the table
			err_status=write_table(xvals, yvals, label1[0], None, file_out, process_name)
		else: # Case of a 2D table
			# Definitions
			file_out=output_dir + "/" + process_name + "_1.priors"
			file_out_plot=file_out + ".jpg"
			# Conversion of sqrt(splitting_a1).cosi, qrt(splitting_a1).sini if requested and necessary
			if convert_a1cosisini == True and \
					(label1[0] == "sqrt(splitting_a1).cosi" and label2[0] == "sqrt(splitting_a1).sini") or \
					(label1[0] == "sqrt(splitting_a1).sini" and label2[0] == "sqrt(splitting_a1).cosi"):
				a1=smcmc1**2 + smcmc2**2
				if (label1[0] == "sqrt(splitting_a1).cosi" and label2[0] == "sqrt(splitting_a1).sini"):
					inc = np.arctan(smcmc2/smcmc1) * 180./np.pi
				if (label1[0] == "sqrt(splitting_a1).sini" and label2[0] == "sqrt(splitting_a1).cosi"):
					inc = np.arctan(smcmc1/smcmc2) * 180./np.pi
				label1[0]="Splitting_a1"
				label2[0]="Inclination"
				smcmc1=a1
				smcmc2=inc
			# Make the table based on np.histogram2d() values
			try:
				Nx_bins = binning[0]
				Ny_bins = binning[1]
			except:
				Nx_bins = binning
				Ny_bins = binning

			x_min = np.min(smcmc1)
			x_max = np.max(smcmc1)
			y_min = np.min(smcmc2)
			y_max = np.max(smcmc2)
			x_bins = np.linspace(x_min, x_max, Nx_bins+1)
			y_bins = np.linspace(y_min, y_max, Ny_bins+1)

			hist2d, x_edges, y_edges = np.histogram2d(smcmc1[:,0], smcmc2[:,0], bins=(x_bins, y_bins), density=True)
			xvals = 0.5 * (x_edges[1:] + x_edges[:-1])	
			yvals = 0.5 * (y_edges[1:] + y_edges[:-1])
			zvals = hist2d.T
			err_status=write_table(xvals, yvals, label1[0], None, file_out, process_name, zvals=zvals, label2=label2[0], unit2=None)
			# Make the 2D plot: Will use ax.hist2d which is supposed to do the exact same as np.histogram2d
			samples=np.zeros((2,len(smcmc1[:,0])))
			samples[0,:]=smcmc1[:,0]
			samples[1,:]=smcmc2[:,0]
			stats=np.zeros((2,6))
			stats[0,:]=stats1
			stats[1,:]=stats2
			labels=[label1[0], label2[0]]
			#fig, ax = plt.subplots(1, figsize=(12, 6))
			fig, ax = plt.subplots()
			im = ax.imshow(zvals, aspect='auto', origin='lower', extent=[x_min, x_max, y_min, y_max])
			plt.colorbar(im)  # Add colorbar if needed

			# Set labels and title
			ax.set_xlabel(label1[0])
			ax.set_ylabel(label2[0])
			ax.set_title('Prior Map')
			fig.savefig(file_out_plot, dpi=300)
	except ValueError as e:
		print("Error while making tabulated priors ")
		print("Error : ", e)
		if os.path.exists(outdir_tmp):
			shutil.rmtree(outdir_tmp)
		exit()
	# Upon normal exit, we erase the temporary directory and its content
	if os.path.exists(outdir_tmp):
		shutil.rmtree(outdir_tmp)
	return 0


def write_table(xvals, yvals, label1, unit1, output_file, process_name, zvals=None, label2=None, unit2=None):
	'''
		Function that handle the formating of the outputs in a 1D or 2D table
		depending on the provided inputs
	'''
	err=False
	if zvals is None: 
		dim=1
	else:
		dim=2
	header="# Tabulated priors for " + process_name + "\n"
	header=header + "# Table created by TAMCMC/tools/convert_fit2prior_table/make_priors_tables.py\n"
	header=header + "# dimensions: " + str(dim)
	if dim == 1:
		labels="!\t" + label1 +"\t PDF"  
		if unit1 is None:
			units="*\t (Unspecified) \t (no_unit)"
		else:
			units="*\t" + unit1 + "\t (no_unit)"
		data=""
		for i in range(len(xvals)):
			data=data + "{:.12f}\t{:.12f}\n".format(xvals[i], yvals[i])
	else:
		labels = "!\t{}\t{}\t PDF".format(label1, label2)
		if unit1 is None and unit2 is None:
			units = "*\t (Unspecified) \t (Unspecified) \t (no_unit)"
		elif unit1 is None and unit2 is not None:
			units = "*\t (Unspecified) \t {}\t (no_unit)".format(unit2)
		elif unit1 is not None and unit2 is None:
			units = "*\t{}\t (Unspecified) \t (no_unit)".format(unit1)
		else:
			units = "*\t{}\t{}\t (no_unit)".format(unit1, unit2)
		data="{:12}\t".format("NA")
		# x-axis
		for i in range(len(xvals)):
			data=data + "{:.12f}\t".format(xvals[i])
		data=data+"\n"
		# y-axis followed by the map in matricial form z(x,y)
		for i in range(len(yvals)):
			data=data + "{:.12f}\t".format(yvals[i])
			for j in range(len(xvals)):
				data = data + "{:.12f}\t".format(zvals[i][j])
			data=data + "\n"
	try:
		with open(output_file, "w") as f:
			f.write(header + "\n")
			f.write(labels + "\n")
			f.write(units + "\n")
			f.write(data)
	except IOError as e:
		print("Error: Unable to open or write to output file", file=sys.stderr)
		print("       Could not write the prior table", file=sys.stderr)
		print("       Debug required", file=sys.stderr)
		err=True
	return err
