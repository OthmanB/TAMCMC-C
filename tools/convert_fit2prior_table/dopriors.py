import numpy as np
import make_priors_tables as m
import os

def scan_directory(directory, prefix):
    subdirectories_found = []
    for path, subdirs, files in os.walk(directory):
        for subdir in subdirs:
            if subdir.startswith(prefix):
                #subdirectories_found.append(os.path.join(path, subdir))
                subdirectories_found.append(subdir)
    if len(subdirectories_found) == 0:
        raise Exception("No subdirectories found with the specified prefix.")
    else:
        return subdirectories_found

def extract_variable_names(param_hdr_file):
    '''
        A small function that extract the variable_names for a parameter hdr file
        created by the TAMCMC process.
    '''
    # Split the file content by newline to get each line as a separate element
    lines = param_hdr_file.split('\n')   
    variable_names = None
    # Loop through each line
    for line in lines:
        # Check if the line starts with the key 'variable_names'
        if line.startswith('! variable_names='):
            # Remove the key and whitespace to get only the variable names
            variable_names = line[len('! variable_names='):].strip().split()
            break
     # Check if variable_names is still None after the loop
    if variable_names is None:
        raise ValueError("'! variable_names=' not found in the param_hdr_file")
    return variable_names

def find_occurrence_index(param_list, string_to_find):
    index = -1
    for i, param_name in enumerate(param_list):
        if param_name == string_to_find:
            if index != -1:
                raise ValueError(f"Multiple occurrences of {string_to_find} found in param_names")
            index = i
    if index == -1:
        raise ValueError(f"{string_to_find} not found in param_names")
    return index


def dopriors(datadir, param_name1, outdir, param_name2=None, prefix='kplr', 
             binning=50, phase='A', chain=0, first_index=0, last_index=-1, 
             period=1, cpp_prg="../../bin/"):
    '''
       Program that look into a directory for all of the outputs of MCMC runs and 
       scans the [process_name]_[analysis_step]_params.hdr in order to resolve the indexes
       of the requested parameters (param_name1 and param_name2) required to build a prior table. 
       This is then used to iteratively write the priors tables of each process into an output directory 
       datadir: The directory containing all of the posterior results for various processes
       param_name1 and param_name2: Name of the parameter that has to be resolved. 
       if param_name2 is not specified, a 1D table will be created. Otherwise a 2D table will be made
       outdir: Output directory where the tables are written
       prefix: Specify the prefix expected for the subdirectories within datadir and that contain the posteriors
    '''
    # Found directories that match the requested syntax
    # It is assumed here that each directory correspond to a process_name,
    # folowing specifications of the TAMCMC program
    process_names=scan_directory(datadir, prefix)
    for i in range(process_names):
        # Important variables definitions
        dir_tamcmc_outputs=datadir + process_names[i]
        param_hdr_file=dir_tamcmc_outputs + '/outputs/' + \
            process_names[i] + "_" + phase + "_params_chain-" + str(chain) + ".hdr"
        
        # Extract from the header the parameter names
        param_names=extract_variable_names(param_hdr_file)
        # Get the occurence location for param_name1 and param_name2
        index1=find_occurrence_index(param_names, param_name1)
        if param_name2 != None:
            index2=find_occurrence_index(param_names, param_name2)
        else:
            index2=None
        # Compute the prior according to parameters above
        m.make_prior_table(index1, binning, outdir, dir_tamcmc_outputs, 
             process_names[i], phase=phase, chain=chain, 
		     first_index=first_index, last_index=last_index, period=period,
		     index2=index2, cpp_prg=cpp_prg)