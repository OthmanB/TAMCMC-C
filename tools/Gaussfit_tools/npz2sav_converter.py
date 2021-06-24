'''
---------------
This file regroups a collection of functions that can be used to perform conversion into various format of the 
spectrum
---------------
'''

import subprocess
import os, fnmatch
import numpy as np

def version():
	return 'v0.1'

# TAKEN FROM prepare_lc_kepler (https://github.com/OthmanB/LCconcat_kepler/prepare_lc_kepler.py)
def get_files_list(dir_in, extension='.npz', prefix='', ID='*'):
	'''
	This routine returns all files that have some extension and prefix and optionally a specific kic encoded in a 9-numerical byte format
	dir_in: The directory to look in for files. This must be a full path
	extension: The extension to look for. This is designed to provide files for prepare_lc_kepler::concatenate_kepler() so that it is fits by default
	prefix: the prefix that defines the names of the files. Kepler data from MAST and KASOC are preceded of 'kplr', which is the default prefix
	ID: The ID number of a specific star. Should not be set if the data structure is 1-star, 1 directory. Instead use the default
	'''
	listOfFiles = os.listdir(dir_in)
	if ID != '*': 
		if extension == '.fits' and prefix == 'kplr': # Ensure that the search incorporate the date sufix that is always in the quarter filenames 
			ID=ID+'-*'
		else:
			ID=ID +'*' # Case where the function is used in another context than Kepler quarters
	pattern = prefix + ID + extension
	matching=[]
	for entry in listOfFiles:
		if fnmatch.fnmatch(entry, pattern):
			matching.append(entry)
			#print(entry)
	return matching


def npzTF2idlsav(dir_npz, npzfile, dir_sav=''):
	'''
		A small program that calls idl in order to save the npz file with the spectrum  
		into a sav file.
		Because the idlpy bridge is not always available (no pip install)
		here we pass by ascii temporary file
		Typically, the conversion is as follow:
			1. Read the npz file
			2. Write the npz content temporary files: Support only the spectrum content
			3. Call an IDL subroutine inside python that will save the file in sav format 
	'''
	if dir_sav == '':
		dir_sav=dir_npz

	current_dir=os.getcwd()
	try:
		idldir=os.environ['IDL_DIR']
	except:
		print('Error: IDL_DIR environment variable does not exist')
		print('       Please set it in the shell script before using this program')
		print('The program will exit now')
		exit()
	ascii_file=npzTF2ascii(dir_npz, npzfile, dir_ascii=dir_sav)
	subprocess.call(["mv", ascii_file, "tmp_in.ascii"])
	subprocess.call([idldir + "/bin/idl", "npz2sav_script.pro"])  # calling idl
	file_list=get_files_list(current_dir, extension='_out.sav', prefix='')
	#print('file_list=', file_list)
	if file_list[0] == '':
		print('Error: Could not find the temporary output sav file')
		print('The program will exit now')
		exit()
	else:
		# If everything went well, we can move the temporary sav file into 
		# the directory dedicated for this sav file
		savfile=npzfile.split('.')[0] + '.sav'
		subprocess.call(["mv", file_list[0], dir_sav+savfile])
	subprocess.call(["rm", "tmp_in.ascii"])

def npzTF2ascii(dir_npz, npzfile, dir_ascii=''):
	
	if dir_ascii == '':
		dir_ascii=dir_npz

	core_file=npzfile.split('.')[0] # remove the extension

	d=np.load(dir_npz+npzfile)
	freq=d['freq']
	power=d['power']
	id_number=str(d['id_number'])
	infos=str(d['infos'])
	#freq=freq, power=power, id_number=id_number, config=config, infos=infos
	ascii_file=dir_ascii + core_file + '.ascii'
	#print(ascii_file)
	try:
		f=open(ascii_file, 'w+')
		f.write('#' + str(infos) + '\n')
		f.write('id: ' + str(id_number) + '\n')
		for i in range(len(freq)):
			string='{:f}   {:f}'.format(freq[i], power[i])
			f.write(string + '\n')
		f.close()
		print('Content of ', npzfile, ' written into', ascii_file)
	except:
		print('Error: prepare_lc_kepler.py::npzTF2ascii(): Failed to write npz spectrum content into ascii file')
		print('The program will exit now')
		#exit()
	return ascii_file

def npzLC2ascii(dir_npz, npzfile, dir_ascii=''):
	
	if dir_ascii == '':
		dir_ascii=dir_npz

	core_file=npzfile.split('.')[0] # remove the extension

	d=np.load(dir_npz+npzfile)
	time=d['lightcurve'][0,:]
	flux=d['lightcurve'][1,:]
	id_number=str(d['id_number'])
	infos=str(d['infos'])
	#freq=freq, power=power, id_number=id_number, config=config, infos=infos
	ascii_file=dir_ascii + core_file + '.ascii'
	print(ascii_file)
	try:
		f=open(ascii_file, 'w+')
		f.write('#' + str(infos) + '\n')
		f.write('id: ' + str(id_number) + '\n')
		for i in range(len(time)):
			string='{:f}   {:f}'.format(time[i], flux[i])
			f.write(string + '\n')
		f.close()
		print('Content of ', npzfile, ' written into', ascii_file)
	except:
		print('Error: prepare_lc_kepler.py::npzLC2ascii(): Failed to write npz lightcurve content into ascii file')
		print('The program will exit now')
		#exit()
	return ascii_file

def npzLC2idlsav(dir_npz, npzfile, dir_sav=''):
	'''
		A small program that calls idl in order to save the npz file with the spectrum  
		into a sav file.
		Because the idlpy bridge is not always available (no pip install)
		here we pass by ascii temporary file
		Typically, the conversion is as follow:
			1. Read the npz file
			2. Write the npz content temporary files: Support only the spectrum content
			3. Call an IDL subroutine inside python that will save the file in sav format 
	'''
	if dir_sav == '':
		dir_sav=dir_npz

	current_dir=os.getcwd()
	try:
		idldir=os.environ['IDL_DIR']
	except:
		print('Error: IDL_DIR environment variable does not exist')
		print('       Please set it in the shell script before using this program')
		print('The program will exit now')
		exit()
	ascii_file=npzLC2ascii(dir_npz, npzfile, dir_ascii=dir_sav)
	subprocess.call(["mv", ascii_file, "tmp_in.ascii"])
	subprocess.call([idldir + "/bin/idl", "npzLC2sav_script.pro"])  # calling idl
	file_list=get_files_list(current_dir, extension='_out.sav', prefix='')
	#print('file_list=', file_list)
	if file_list[0] == '':
		print('Error: Could not find the temporary output sav file')
		print('The program will exit now')
		exit()
	else:
		# If everything went well, we can move the temporary sav file into 
		# the directory dedicated for this sav file
		savfile=npzfile.split('.')[0] + '.sav'
		subprocess.call(["mv", file_list[0], dir_sav+savfile])
	subprocess.call(["rm", "tmp_in.ascii"])

print('npz2sav_converter.py, version ', version(), ' loaded')
