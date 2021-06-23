'''
	A procedure that handle the generation of the power spectrum
	It creates a copy of the file into a sav format for compatibility
	with my IDL codes
'''
from npz2sav_converter import *

dir_in='/Volumes//homes/dataonly/Kepler/Siddarth_products_newonly/'
dir_out='/Volumes//homes/dataonly/Kepler/Siddarth_products_newonly/sav_files/'

all_files=get_files_list(dir_in, extension='_TF.npz', prefix="")
for k in all_files:
	print(' Converting ', k , ' ...')
	npzTF2idlsav(dir_in, k, dir_sav=dir_out)
