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

def extract_variable_names(param_hdr_file, merge=True):
    '''
        A small function that extract the variable_names, constant_names and relax for a parameter hdr file
        created by the TAMCMC process.
        If merge == True, it also reconstruct the full vector of values
    '''
    # Read all of the file at once without formating
    try:
        with open(param_hdr_file, "r") as f:
            data = f.read()
    except IOError:
        print("Error: The file {} could not be opened!".format(param_hdr_file))
        exit()
    # Split the file content by newline to get each line as a separate element
    lines = data.split('\n')   
    variable_names = None
    constant_names = None
    relax = None
    # Loop through each line
    for line in lines:
        # Check if the line starts with the key 'variable_names'
        if line.startswith('! variable_names='):
            # Remove the key and whitespace to get only the variable names
            variable_names = line[len('! variable_names='):].strip().split()
        # Check if the line starts with the key 'constant_names'
        if line.startswith('! constant_names='):
            # Remove the key and whitespace to get only the constant names
            constant_names = line[len('! constant_names='):].strip().split()
        # Check if the line starts with the key 'relax'
        if line.startswith('! relax='):
            # Remove the key and whitespace to get only the relax array
            relax = line[len('! relax='):].strip().split()
     # Check if variable_names is still None after the loop
    if variable_names is None:
        raise ValueError("'! variable_names=' not found in the param_hdr_file")
    if merge == False:
        return variable_names, constant_names, relax
    else:
        full_vector = []
        k=0
        l=0
        for i in range(len(relax)):
            if relax[i] == '1':
                full_vector.append(variable_names[k])
                k=k+1
            else:
                full_vector.append(constant_names[l])
                l=l+1
        return full_vector, relax
    
def find_occurrence_index(param_list, string_to_find, silent=False):
    index = -1
    for i, param_name in enumerate(param_list):
        #print(i , param_name, string_to_find)
        if param_name == string_to_find:
            if index != -1:
                raise ValueError("Multiple occurrences of {} found in param_names".format(string_to_find))
            index = i
    if index == -1:
        #raise ValueError("{} not found in param_names".format(string_to_find))
        if silent == False:
            print("{} not found in param_names".format(string_to_find))
    return index


def dopriors(datadir, param_name1, outdir, param_name2=None, prefix='kplr', 
             binning=50, phase='A', chain=0, first_index=0, last_index=-1, 
             period=1, cpp_prg="../../bin/", search_alternate_names=True):
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
       search_alternate_names: If True and that the param_name1 or param_name2 was not found, consider alternative names
            known to be used in other version of the TAMCMC code. Supported alternate names:
                - Splitting_a1, a1_0
                - sqrt(splitting_a1).cosi, sqrt(splitting_a1).sini ==> Conversion is performed to a1,Inc in that case
    '''
    # Definitions of alternative names. 
    # It is an array of lists. Each array entry provides a combination of all equivalent names
    alternative_names=[
        ["splitting_a1", "Splitting_a1", "a1_0", "sqrt(splitting_a1).cosi"],
        ["Inclination", "inclination",  "sqrt(splitting_a1).sini"]
        ]
    # Found directories that match the requested syntax
    # It is assumed here that each directory correspond to a process_name,
    # folowing specifications of the TAMCMC program
    process_names=scan_directory(datadir, prefix)
    for i in range(len(process_names)):
        print(" [{}] Processing {}...".format(i, process_names[i]))
        # Important variables definitions
        dir_tamcmc_outputs=datadir + process_names[i]
        param_hdr_file=dir_tamcmc_outputs + '/outputs/' + \
            process_names[i] + "_" + phase + "_params.hdr"
        # Extract from the header the parameter names
        param_names, relax=extract_variable_names(param_hdr_file)
        # Get the occurence location for param_name1 and param_name2
        index1=find_occurrence_index(param_names, param_name1, silent=True)
        if relax[index1] == '0': # fixed parameters cannot be considered as valid
            index1=-1
        # Look if alternate names do lead to index1 != -1 (ie names that can be found)
        # If it is, then check if an available alternative name is not suitable (but only if the user requests it)
        if index1 == -1 and search_alternate_names == False:
            raise ValueError("'{}' not found in param_names".format(param_name1))
        else:     
            # Search within the variable alternative_names if nothing is equivalent to param_name1 
            # 1. First check if alternative_names do at least contain the currently searched param_name1
            index_ref=-1
            for k in range(len(alternative_names)):
                index_ref=find_occurrence_index(alternative_names[k], param_name1)
                if index_ref !=-1:
                    break
            if index_ref == -1:
                print("Available alternative names: ")
                for a in alternative_names:
                    print("  - {}", a)
                raise ValueError("No alternative names found for '{}' in the above list!".format(param_name1))
            # 2 Second, if we found it, we use the identified table (indexed with k) that list the alternates 
            for alt_name in alternative_names[k]:
                index = find_occurrence_index(param_names, alt_name, silent=True)
                if index != -1:
                    if relax[index] == '1': # Only update the index1 if the found parameter is a variable
                        index1 = index
                        break
            if index1 == -1:
                raise ValueError("'{}' not found in param_names, despite having looked for alternative names! Check for a typo".format(param_name1))
        if param_name2 != None:
            index2=find_occurrence_index(param_names, param_name2, silent=True)
            if relax[index2] == '0': # fixed parameters cannot be considered as valid
                index2=-1
            if index2 == -1 and search_alternate_names == False:
                raise ValueError("{} not found in param_names".format(param_name2))
            else:     
                # Search within the variable alternative_names if nothing is equivalent to param_name2
                # 1. First check if alternative_names do at least contain the currently searched param_name2
                index_ref=-1
                for k in range(len(alternative_names)):
                    index_ref=find_occurrence_index(alternative_names[k], param_name2, silent=True)
                    if index_ref !=-1:
                        break
                if index_ref == -1:
                    print("Available alternative names: ")
                    for a in alternative_names:
                        print("  - {}", a)
                    raise ValueError("No alternative names found for '{}' in the above list!".format(param_name2))
                # 2. Second, if we found it, we use the identified table (indexed with k) that list the alternates 
                for alt_name in alternative_names[k]:
                    index = find_occurrence_index(param_names, alt_name, silent=True)
                    if index != -1:
                        if relax[index] == '1': # Only update the index2 if the found parameter is a variable
                            index2 = index
                            break
                if index2 == -1:
                    raise ValueError("{} not found in param_names, despite having looked for alternative names! Check for a typo".format(param_name2))
        else:
            index2=None
        # 3. NOTE: any specific rule that is required. For example converting sqrt(a1).cosi, sqrt(a2).sini has to be 
        # inside make_prior_table
        #
        # Compute the prior according to parameters above
        m.make_prior_table(index1, binning, outdir, datadir, 
             process_names[i], phase=phase, chain=chain, 
		     first_index=first_index, last_index=last_index, period=period,
		     index2=index2, cpp_prg=cpp_prg, convert_a1cosisini=True)
        
def dopriors_quicktest():
    datadir='/Volumes/homes/dataonly/Work-2023/ajAlm_project/data_tests/outputs_1001/'
    #datadir='/Volumes/homes/dataonly/Work-2023/16CygAB_analysis_sqrta1cosi/'
    param_name1='a1_0'
    param_name2='Inclination'
    outdir='/Volumes/homes/dataonly/Work-2023/tmp/'
    period=7
    dopriors(datadir, param_name1, outdir, param_name2=param_name2, prefix='kplr', 
             binning=50, phase='A', chain=0, first_index=0, last_index=-1, 
             period=period, cpp_prg="../bin/", search_alternate_names=True)