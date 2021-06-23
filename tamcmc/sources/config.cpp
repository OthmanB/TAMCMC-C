/*
 * Config.cpp
 *
 * Contains all kind of functions
 * used to process and/or encapsulate data
 * 
 *  Created on: 20 Mar 2016
 *      Author: obenomar
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "io_ms_global.h"
#include "string_handler.h"

Config::Config(std::string current_path, std::string cfg_file_in, std::string cfg_file_errors, 
			   std::string cfg_models_ctrl_file_in, std::string cfg_priors_ctrl_file_in, std::string cfg_likelihoods_ctrl_file_in,
			   std::string cfg_primepriors_ctrl_file_in){ // The constructor

	Data_Basic listoutputs;
	
	working_dir=current_path;
	cfg_file=cfg_file_in;
	errordefault_file=cfg_file_errors;
	cfg_models_ctrl_file=cfg_models_ctrl_file_in;
	cfg_priors_ctrl_file=cfg_priors_ctrl_file_in;
	cfg_likelihoods_ctrl_file=cfg_likelihoods_ctrl_file_in;
	cfg_primepriors_ctrl_file=cfg_primepriors_ctrl_file_in;
	
	// ---- Read the configuration files and perform the setup ----
	const bool verbose=0;
	std::cout << "       - Reading the file with the default configuration: " << std::endl;
	std::cout << "                 " << cfg_file_in << "..." << std::endl;
	read_cfg_file(verbose); // read the default configuration file... because it is not necessarily the final configuration, do not verbose
	std::cout << "       - Reading the file with the default error values on parameters: " << std::endl;
	std::cout << "                 " << cfg_file_errors << "..." << std::endl;
	read_defautlerrors(verbose); // read the values that are used to initialized the covariance matrix
	std::cout << "       - Reading the file listing model functions: " << std::endl;
	std::cout << "                 " << cfg_models_ctrl_file << "..." << std::endl;
	listoutputs=read_listfiles(cfg_models_ctrl_file, verbose);
	modeling.models_case_list_ctrl=listoutputs.vecXi;
	modeling.models_list_ctrl=listoutputs.strarr;
		
	std::cout << "       - Reading the file listing (meta)prior functions: " << std::endl;
	std::cout << "                 " << cfg_priors_ctrl_file << "..." << std::endl;
	listoutputs=read_listfiles(cfg_priors_ctrl_file, verbose);
	modeling.priors_case_list_ctrl=listoutputs.vecXi;
	modeling.priors_list_ctrl=listoutputs.strarr;
	
	std::cout << "       - Reading the file listing likelihood functions: " << std::endl;
	std::cout << "                 " << cfg_likelihoods_ctrl_file << "..." << std::endl;
	listoutputs=read_listfiles(cfg_likelihoods_ctrl_file, verbose);
	modeling.likelihoods_case_list_ctrl=listoutputs.vecXi;
	modeling.likelihoods_list_ctrl=listoutputs.strarr;

	std::cout << "       - Reading the file listing primitive prior functions : " << std::endl;
	std::cout << "                 " << cfg_primepriors_ctrl_file << "..." << std::endl;
	listoutputs=read_listfiles(cfg_primepriors_ctrl_file, verbose);
	modeling.primepriors_case_list_ctrl=listoutputs.vecXi;
	modeling.primepriors_list_ctrl=listoutputs.strarr;

	modeling.slice_ind=0; // Default there is only one slice of data that is analysed
	
//	for(int i=0; i< modeling.primepriors_case_list_ctrl.size(); i++){
//		std::cout << 	 modeling.primepriors_case_list_ctrl[i] << "   " << modeling.primepriors_list_ctrl[i] << std::endl;
//	}
	data.data.xrange.resize(2);
	data.data.xrange.setConstant(-9999); // Initialize the data range to a dummy value easily recognizable
	// ------------------------------------------------------------
}

Config::Config(){ // The empty constructor

}

void Config::setup(const int slice_ind){ 
	
	long imin, imax;


    // ---- Read the data ----
    std::string delimiter=" ";
    std::cout << "Data file: " << data.data_file << std::endl;
    Data_Nd data_in=read_data_ascii_Ncols(data.data_file, delimiter, data.verbose_data);
    data.data_all=data_in; // save the whole data file into the configuration class
    
    // ---- Reading the model-specific configuration files ----
	modeling.slice_ind=slice_ind;
    std::cout << " ---------- " << std::endl;
    read_inputs_files(); // Here we read the configuration files (e.g. the .MCMC file)
    std::cout << "       - Converting prior names into integers..." << std::endl;
    modeling.inputs.priors_names_switch=convert_priors_names_to_switch(modeling.inputs.priors_names);
    std::cout << "       - Converting model function names into integers..." << std::endl;
    modeling.model_fct_name_switch=convert_model_fct_name_to_switch(modeling.model_fct_name);
    std::cout << "       - Converting  likelihood function names into integers..." << std::endl;
    modeling.likelihood_fct_name_switch=convert_likelihood_fct_name_to_switch(modeling.likelihood_fct_name);
    std::cout << "       - Converting prior function names into integers..." << std::endl;
    modeling.prior_fct_name_switch=convert_prior_fct_name_to_switch(modeling.prior_fct_name);
	
	std::cout << "       - Handling data inputs..." << std::endl;
    // ----- Define which columns are containing the x values, the y values and ysig_ind ----
	if(data.data.xrange[0] == -9999 && data.data.xrange[1] == -9999){ // Case where no range was given in the cfg file ==> Take all
		imin=0;
		imax=data_in.data.rows();
	} else{ // Case where a range was set ==> select the appropriate data range
		imin=0;
		while(imin<data_in.data.rows() && data_in.data(imin, data.x_col) <data.data.xrange[0]){
			imin=imin+1.;
			//std::cout << "imin=" << imin << "   data_in.data(imin, data.x_col)= " << data_in.data(imin, data.x_col) << std::endl;
        }
		if(imin >= data_in.data.rows()){
			std::cout << "Warning: Found that xmin > max(data.x)" << std::endl;
			std::cout << "         The requested xrange[0] is inconsistent with the given data" << std::endl;
			std::cout << "         Cannot proceed further. Check the configuration file and the data file " << std::endl;				std::cout << "         The program will stop now" << std::endl;
			exit(EXIT_FAILURE);
		}
		//std::cout << "imin=" << imin << std::endl;
		imax=imin;
		while(imax<data_in.data.rows() && data_in.data(imax, data.x_col)<data.data.xrange[1]){
			imax=imax+1.;
		}
		if(imax > data_in.data.rows()){
			imax=data_in.data.rows();
			std::cout << "Warning: Found that xmax > max(data.x) OR xmax < xmin" << std::endl;
			std::cout << "         ==> xmax fixed to max(data.x)" << std::endl;
			std::cout << " Stop for debug purpose" <<std::endl;
			exit(EXIT_SUCCESS);
			std::cout << "         Proceeding... " << std::endl;
		}
	}

	if(data.x_col >=0){
		data.data.x=data_in.data.col(data.x_col).segment(imin, imax-imin);
		data.data.xlabel=data_in.labels[data.x_col];
		data.data.xunit=data_in.units[data.x_col];		
	}
	//std::cout << "data.data.x.minCoeff()=" << data.data.x.minCoeff() << std::endl;
	//std::cout << "data.data.x.maxCoeff()=" << data.data.x.maxCoeff() << std::endl;
	
	if(data.y_col >=0){
		std::cout << "Importing data...";
		if(data_in.data.cols() > data.y_col){
			data.data.y=data_in.data.col(data.y_col).segment(imin, imax-imin);
			std::cout << "labels...";
			data.data.ylabel=data_in.labels[data.y_col];
			std::cout << "units";
			data.data.yunit=data_in.units[data.y_col];
			std::cout << "...Done" << std::endl;
		}else{
			std::cout << "Fatal Error: The data file has only " << data_in.data.cols() << " columns but you specified that y-data are in column " << data.y_col+1 << " (y_col=" << data.y_col << ")" << std::endl;
			std::cout << "             Cannot proceed. The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
	} else{
		std::cout << "Fatal error: You need to specify a column for y-data" << std::endl;
		std::cout << "             Cannot proceed. The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	if(data.ysig_col >= 0){
		data.data.sigma_y=data_in.data.col(data.ysig_col).segment(imin, imax-imin);
	} else {
		std::cout << "Warning: Config::Config.data.data.sigma_y was not specified in the data file" << std::endl;
		std::cout << "         In case of a chi_square fit, this variable is required" << std::endl;
		std::cout << "         To avoid code failure, this value will be fixed to 1" << std::endl;
		std::cout << "         If you wish to use sigma_y != 1, then add a column in you data file: " << data.data_file << std::endl;
		VectorXd tmp(data.data.x.size());
		data.data.sigma_y=tmp.setConstant(1);
	}
	data.data.header=data_in.header;
	data.data.Nx=data.data.x.size();


    // --- Finishing the Setup ---
	if(outputs.do_restore_last_index == 1 && outputs.do_restore_proposal == 0){
		std::cout << "Warning: Problem in the configuration file" << std::endl;
		std::cout << "         If do_restore_last_index is 1, do_restore_proposal must also be 1" << std::endl;
		std::cout << "         Here, this is not the case: do_restore_last_index=1 but do_restore_proposal=0" << std::endl;
		std::cout << "         The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	if(outputs.do_restore_variables == 1 || outputs.do_restore_proposal == 1){
		read_restore_files(); // This fils the structure restored_vals ONLY if do_restore (proposal or variables) are set to true
		//std::cout << "   [4] Performing consistency checks..." << std::endl;
		std::cout << std::endl << "       Note: All values listed above superseed the configuration as set in the configuration file" <<std::endl;
	} else {
		restored_vals.iteration=0;
	}
	
	std::cout << "Initial configuration done" << std::endl;
}

void Config::reset(){

	working_dir=""; 
	cfg_file=""; // The configuration file for everything
	MALA.target_acceptance=0.;
	MALA.c0=0;
	MALA.epsilon1=0;
	MALA.epsi2=0;
	MALA.A1=0;
	MALA.Nt_learn.resize(0);
	MALA.periods_learn.resize(0);
	// ---------------- Controls for the Langevin scheme  -----------------
	MALA.use_drift=0; // Not operational
	MALA.delta=0; // Not operational
	MALA.delta_x=0; // Not operational
	// ----------- Controls for the parallel tempering ---------
	MALA.Nchains=0;
	MALA.dN_mixing=0;
	MALA.lambda_temp=0;
	// --------- Experimental extra parameters ----------
	MALA.proposal_type="";


	// --------------- Model parameters --------------
	modeling.slice_ind=0;
	modeling.model_fct_name="";
	modeling.likelihood_fct_name="";
	modeling.prior_fct_name="";
	modeling.likelihood_params=0; // Parameters that could define the likelihood. e.g., in the case of a chi(2,2p) statistics (chi22p function), we need p
	modeling.cfg_model_file=""; // Contains the initial guesses and the priors in an ASCII format. To each model_fct_name, a given format is expected
	modeling.inputs.inputs_names.resize(0);
	modeling.inputs.priors_names.resize(0);
	modeling.inputs.inputs.resize(0);
	modeling.inputs.relax.resize(0);
	modeling.inputs.priors.resize(0,0);
	modeling.inputs.plength.resize(0);

	// all the configuration required for the Data
	data.verbose_data=0; // Should we show them on screen before proceeding?
	data.data_file=""; // The file containing all the data
	data.type_data="";
	data.x_col=0;
	data.y_col=0;
	data.ysig_col=0;
	data.data_all.data.resize(0,0); // this is the raw outputs of the input data file
	data.data_all.header.resize(0); // this is the raw outputs of the input data file
	data.data_all.labels.resize(0); // this is the raw outputs of the input data file
	data.data_all.units.resize(0); // this is the raw outputs of the input data file
	data.data.x.resize(0); // this is a formated outputs suitable for the MCMC analysis
	data.data.y.resize(0);
	data.data.sigma_y.resize(0);
	data.data.Nx=0;
	data.data.xlabel="";
	data.data.ylabel="";
	data.data.xunit="";
	data.data.yunit="";
	data.data.header.resize(0);

	outputs.Nsamples=0;
	outputs.Nbuffer=0;

	outputs.erase_old_files=0;

	// --------- Variables that decide which outputs we do write ------------
	outputs.get_params=0;  //if set to 1, save the parameters of the fitted model. This option might always be set to 1
	outputs.get_statcriteria=0; //if set to 1, save the statistical parameters (logLikelihood, logPrior, logPosterior, ...). Recommended to always set it to 1	
	outputs.get_proposal_params=0;  //if set to 1, save the proposal parameters (sigmas, mus, covarmats) of the parallel chains
	outputs.get_parallel_tempering=0; //if set to 1, save the parallel tempering information (swaped chain, probability, ...)
	outputs.get_models=0;  //if set to 1, save the models for all chains. WARNING: Might greatly slow-down the process and will take a HUGE amount of space

	//  -----------  Variables defining the output files. Do not write the filename full path and extension -----------
	outputs.output_root_name="";  // The root of the names of the outputs file(e.g. the star identifier such as the KIC number)
	outputs.dir_out="";  // Output directory for all outputs
	outputs.file_out_format=""; // PLAIN ASCII or BINARY or HDF5 (NOT IMPLEMENTED YET)

	outputs.params_txt_fileout="";  // extension name for the txt file in output for the parameters
	outputs.proposal_txt_fileout=""; // extension name for the txt file in output for the proposal parameters
	outputs.parallel_tempering_txt_fileout=""; // extension name for the txt file in output for the parallel tempering parameters
	outputs.model_txt_fileout=""; // extension name for the txt file in output for the models
	outputs.stat_txt_fileout=""; // extension name for the txt file in output for the statistical info such as the logLikelihood, logPrior and logPosterior
	outputs.acceptance_txt_fileout=""; // extension name for the txt file with the acceptance rates

	//outputs.hdf5_fileout=""; // extension name for the hdf5 file. All data are put in one single HDF5 if this option is chosen
		
	// ------------ Information required for the completion of an aborted run OR for pursuing a complete run  -------------
	outputs.do_restore_proposal=0;
	outputs.do_restore_variables=0;
	outputs.do_restore_last_index=0;
	outputs.restore_dir="";
	outputs.restore_file_in="";
	outputs.restore_file_out="";

	// --------- Variables that decide which outputs we do plot ------------
	diags.chains_diags=0;
	diags.evidence_diags=0;
	diags.pdfs_diags=0;

	diags.output_dir="";
	diags.output_root_name="";
	diags.file_chains_diags="";
	diags.file_evidence_diags="";
	diags.file_pdfs_diags="";

	diags.model_initial_diags=0;
	diags.model_buffer_diags=0;
	diags.model_final_diags=0;

	diags.file_model_init_diags="";
	diags.file_model_buffer_diags="";
	diags.file_model_final_diags="";

	diags.data_scoef1=0; 
	diags.data_scoef2=0;
	diags.show_original_data=0;
	diags.evidence_interpolation_factor=1;
	diags.Nclasses=0;

	// Contains the values that were restored but that do not belong to this class (variables, sigma, mu and the covarmat)
	restored_vals.variable_names.resize(0);
	restored_vals.iteration=0; // by default it is 0. If we restore from a previous run, this is changed to the restore file value
	restored_vals.Nchains=0;
	restored_vals.Nvars=0;
	restored_vals.vars.resize(0,0);
	restored_vals.sigma.resize(0);
	restored_vals.mu.resize(0,0);
	restored_vals.covarmat=initialize_3dMatrix(0, 0, 0); // depth, Nrows, Ncols
	std::cout << restored_vals.covarmat[0]->size() << std::endl;

	std::cout << "Initial configuration parameters, re-initialised to 0/0-sized values/vectors" << std::endl;

}

void Config::read_inputs_priors_MS_Global(){

	bool verbose=1;
	MCMC_files iMS_global;
	Input_Data in_vals;
	//std::cout << "Before read_MCMC" << std::endl;
	std::cout << "  - Reading the MCMC file: " << modeling.cfg_model_file << "..." << std::endl;
	iMS_global=read_MCMC_file_MS_Global(modeling.cfg_model_file, 0); // Read the MCMC file, with verbose=0 here
	data.data.xrange=iMS_global.freq_range; // Load the wished frequency range into the data structure (contains the spectra)
	std::cout << "   - Preparing input and priors parameters..." << std::endl;
    in_vals=build_init_MS_Global(iMS_global, verbose, data.data_all.data(2, data.x_col)-data.data_all.data(1, data.x_col)); // Interpret the MCMC file and format it as an input structure
	in_vals.priors_names_switch=convert_priors_names_to_switch(in_vals.priors_names); // Determine the switch cases from the prior names
	modeling.inputs=in_vals;
	modeling.model_fct_name=in_vals.model_fullname;
	std::cout << "Setup according to the MCMC configuration file finished" << std::endl;
}

void Config::read_inputs_priors_asymptotic(){

	bool verbose=1;
	MCMC_files i_asymptotic;
	Input_Data in_vals;
	//std::cout << "Before read_MCMC" << std::endl;
	std::cout << "  - Reading the MCMC file: " << modeling.cfg_model_file << "..." << std::endl;
	i_asymptotic=read_MCMC_file_asymptotic(modeling.cfg_model_file, 0); // Read the MCMC file, with verbose=0 here.. 
	data.data.xrange=i_asymptotic.freq_range; // Load the wished frequency range into the data structure (contains the spectra)
	std::cout << "   - Preparing input and priors parameters..." << std::endl;
    in_vals=build_init_asymptotic(i_asymptotic, verbose, data.data_all.data(2, data.x_col)-data.data_all.data(1, data.x_col)); // Interpret the MCMC file and format it as an input structure
	in_vals.priors_names_switch=convert_priors_names_to_switch(in_vals.priors_names); // Determine the switch cases from the prior names
	modeling.inputs=in_vals;
	modeling.model_fct_name=in_vals.model_fullname;
	std::cout << "Setup according to the MCMC configuration file finished" << std::endl;
}


void Config::read_inputs_priors_local(){

	const bool verbose=1;
	MCMC_files i_local;
	Input_Data in_vals;
		
	//std::cout << "Before read_MCMC" << std::endl;
	std::cout << "  - Reading the MCMC file: " << modeling.cfg_model_file << "..." << std::endl;
	i_local=read_MCMC_file_local(modeling.cfg_model_file, modeling.slice_ind, 0); // Read the MCMC file, with verbose=0 here
	
	data.data.xrange=i_local.freq_range; // Load the wished frequency range into the data structure (contains the spectra)
	
//	std::cout << "Stop in Config::read_inputs_priors_local: Need to be implemented from here..." << std::endl;
//	exit(EXIT_SUCCESS);
	
	std::cout << "   - Preparing input and priors parameters..." << std::endl;
    in_vals=build_init_local(i_local, verbose, data.data_all.data(2, data.x_col)-data.data_all.data(1, data.x_col)); // Interpret the MCMC file and format it as an input structure
	in_vals.priors_names_switch=convert_priors_names_to_switch(in_vals.priors_names); // Determine the switch cases from the prior names
	modeling.inputs=in_vals;
	modeling.model_fct_name=in_vals.model_fullname;
	std::cout << "Setup according to the MCMC configuration file finished" << std::endl;
}


VectorXi Config::convert_priors_names_to_switch(const std::vector<std::string> p_names){
/* 
 * This function convert the p_names string into a set of pre-defined integer values
 * This is specifically intended to convert priors_names into a limited set of cases
 * that can be handled by the 'switch(x)' and case statement. The switch statement
 * is faster than the if/else and is prefered when used in loops.
 * See stats_dictionary.cpp to find the list of available functions for the priros
*/

	bool passed;
	int i, j;
	const int Npriors=modeling.primepriors_list_ctrl.size();
	VectorXi switch_names(p_names.size()); // The output vector

	for (i = 0; i< p_names.size(); i++){	
		passed=0;
		for(j = 0; j<Npriors;j++){
			if (p_names[i] ==  modeling.primepriors_list_ctrl[j]){
				switch_names[i]= modeling.primepriors_case_list_ctrl[j];
				passed=1;
			}
		}
		if (passed == 0){
			msg_handler("", "prior_fctname", "Config::convert_priors_names_to_switch()", p_names[i], 1);
		}
	}
	return switch_names;


}


int Config::convert_model_fct_name_to_switch(const std::string model_name){
/* 
 * This function convert the model_name string into a pre-defined integer value
 * This is specifically intended to convert model_name into a limited set of cases
 * that can be handled by the 'switch(x)' and case statement. The switch statement
 * is faster than the if/else and is prefered when used in loops.
 * See model_def.cpp for the actual use of the switch statement.
*/

	bool passed;
	int i;
	const int Nmodels=modeling.models_list_ctrl.size();
	int switch_name=-9999; // The output case index

	passed=0;

	for (i = 0; i< Nmodels; i++){	
		if (model_name ==  modeling.models_list_ctrl[i]){
			switch_name= modeling.models_case_list_ctrl[i];
			passed=1;
		}
	}
    if (passed == 0){
		msg_handler("", "model_name", "Config::convert_model_fct_name_to_switch()", model_name, 1);
	}	
	//std::cout << model_name << "    ==> " << switch_name << std::endl;
	//exit(EXIT_SUCCESS);
	return switch_name;
}

int Config::convert_model_fct_name_to_switch(const std::string model_name, const Data_Basic models_ctrl){
/* 
 * This function convert the model_name string into a pre-defined integer value
 * This is specifically intended to convert model_name into a limited set of cases
 * that can be handled by the 'switch(x)' and case statement. The switch statement
 * is faster than the if/else and is prefered when used in loops.
 * See model_def.cpp for the actual use of the switch statement.
 * HERE : models_list_ctrl is a list of models usually listed in a text fle and read with Config::read_listfiles()
*/

	bool passed;
	int i;
	const int Nmodels=models_ctrl.vecXi.size();
	int switch_name=-9999; // The output case index

	passed=0;

	for (i = 0; i< Nmodels; i++){	
		if (model_name ==  models_ctrl.strarr[i]){
			switch_name= models_ctrl.vecXi[i];
			passed=1;
		}
	}
    if (passed == 0){
		msg_handler("", "model_name", "Config::convert_model_fct_name_to_switch()", model_name, 1);
	}	
	//std::cout << model_name << "    ==> " << switch_name << std::endl;
	//exit(EXIT_SUCCESS);
	return switch_name;
}

int Config::convert_prior_fct_name_to_switch(const std::string prior_name){
/* 
 * This function convert the prior_name string into a pre-defined integer value
 * This is specifically intended to convert prior_fct_name into a limited set of cases
 * that can be handled by the 'switch(x)' and case statement. The switch statement
 * is faster than the if/else and is prefered when used in loops.
 * See model_def.cpp for the actual use of the switch statement.
*/
	bool passed;
	int i;
	const int Npriors=modeling.priors_list_ctrl.size();
	int switch_name=-9999; // The output case index

	passed=0;

	for (i = 0; i< Npriors; i++){
		if (prior_name ==  modeling.priors_list_ctrl[i]){
			switch_name= modeling.priors_case_list_ctrl[i];
			passed=1;
		}
	}
	if (passed == 0){
		msg_handler("", "prior_name", "Config::convert_prior_fct_name_to_switch()", prior_name, 1);
	}
	//std::cout << prior_name << "    ==> " << switch_name << std::endl;
	//exit(EXIT_SUCCESS);

	return switch_name;
}

int Config::convert_likelihood_fct_name_to_switch(const std::string likelihood_name){
/* 
 * This function convert the model_name string into a pre-defined integer value
 * This is specifically intended to convert model_name into a limited set of cases
 * that can be handled by the 'switch(x)' and case statement. The switch statement
 * is faster than the if/else and is prefered when used in loops.
 * See model_def.cpp for the actual use of the switch statement.
*/
	
	bool passed;
	int i;
	const int Npriors=modeling.likelihoods_list_ctrl.size();
	int switch_name=-9999; // The output case index

	passed=0;

	for (i = 0; i< Npriors; i++){	
		if (likelihood_name ==  modeling.likelihoods_list_ctrl[i]){
			switch_name= modeling.likelihoods_case_list_ctrl[i];
			passed=1;
		}
	}
	if (passed == 0){
		msg_handler("", "likelihood_name", "Config::convert_likelihood_fct_name_to_switch()", likelihood_name, 1);
	}
	//std::cout << likelihood_name << "    ==> " << switch_name << std::endl;

	return switch_name;

}

Data_Nd Config::read_data_ascii_Ncols(const std::string file_in_name, const std::string delimiter, const bool verbose_data){
/*
 * This function read an input file (file_in_name) which may or may not contain a header, labels, units indicator.
 * The number of column can be as small as 1
*/

    int cpt, Nrows;
    std::string line, line0, subline0; //token
    std::vector<std::string> header, labels, units, data_str;
    long double tmp_val;
    MatrixXd data;
    Data_Nd all_data_out; // The structure that encapsulate all the data, the header, the labels and units
    int data_Maxsize=1000000;
 
    std::ifstream file_in;
    std::cout << "Reading the Data File..." << std::endl;
    std::cout << "  Assumptions for the data file: " << std::endl;
    std::cout << "       - header lines appear on the top of the file and are indicated by a # for first character. The header is however optional (then no # are found)" << std::endl;
    std::cout << "       - labels appear in one single line that is just after the  header and is indicated by a ! for first character. Labels are optional (then no ! are found)" << std::endl;
    std::cout << "       - units appear in one single line that is just after the  labels or the header (if the labels are missing) and is indicated by a * for first character. Units are optional (then no * are found)" << std::endl;
    std::cout << "       - Maximum number of lines for the data: " << data_Maxsize << std::endl;
    file_in.open(file_in_name.c_str());
    if (file_in.is_open()) {
	std::cout << "...processing lines" << std::endl;

		// [1] Get the header
		cpt=0;
		std::getline(file_in, line0);
		line0=strtrim(line0); // remove any white space at the begining/end of the string
		subline0=strtrim(line0.substr(0, 1)); // pick the first character
		subline0=subline0.c_str();
		if (subline0 == "#"){
			while(subline0 == "#"){		
				header.push_back(strtrim(line0.substr(1, std::string::npos))); // add all characters except the first one (and any space at the begining)

				std::getline(file_in, line0);
				line0=strtrim(line0); // remove any white space at the begining/end of the string
				subline0=strtrim(line0.substr(0, 1)); // pick the first character
				subline0=subline0.c_str();
				
				cpt=cpt+1;
			}
			std::cout << "   [1] " << cpt << " header lines found..." << std::endl;
		} else{
			header.push_back("");
			std::cout << "   [1] Header not found. Header vector set to a blank vector<string> of size 1. Pursuing operations..." << std::endl;
		}

		// [2] Read the labels... these are expected just after the header... If not found, then labels is left blank of size 1
		// No need to read line... already read by the previous loop
		line0=strtrim(line0); // remove any white space at the begining/end of the string
		subline0=strtrim(line0.substr(0, 1)); // pick the first character
		subline0=subline0.c_str();
		if(subline0 == "!"){
			
			labels=strsplit(strtrim(line0.substr(1,std::string::npos)), " \t"); // remove either when you found a white space or a tabulation
			std::cout << "   [2] " << labels.size() << " labels found..." << std::endl;

			for(int i=0; i<labels.size();i++){
				labels[i]=strtrim(labels[i]); // remove spaces at begining/end of each substring
				std::cout << "         - " << labels[i] << std::endl;
			}

		} else {
			labels.push_back(""); // label for x
			labels.push_back(""); // label for y1
			labels.push_back(""); // label for y2
			labels.push_back(""); // label for y3
			labels.push_back(""); // label for y4
			std::cout << "   [2] No labels found. Label vector set to a blank vector<string> of size 5. Pursuing operations..." << std::endl;
		}
		
		// [3] Read the units... these are expected just after the labels... If not found, then units is left blank of size 1
		// Read a new line only if we found a label. Otherwise, it is not necessary...
		if (labels[0] != ""){std::getline(file_in, line0);}
		line0=strtrim(line0); 
		subline0=strtrim(line0.substr(0, 1)); 
		subline0=subline0.c_str();
		if(subline0 == "*"){
			units=strsplit(strtrim(line0.substr(1,std::string::npos)), " \t"); 
			std::cout << "   [3] " << units.size() << " units found..." << std::endl;
			for(int i=0; i<units.size();i++){
				units[i]=strtrim(units[i]);
				std::cout << "         - " << units[i] << std::endl;
			}
		} else {
 			units.push_back(""); // unit for x
			units.push_back(""); // label for y1
			units.push_back(""); // label for y2
			units.push_back(""); // label for y3
			units.push_back(""); // label for y4
			std::cout << "   [3] No units found. Units vector set to a blank vector<string> of size 5." << std::endl;
		}

		// [4] Read the data...
		std::cout <<  "   [4] Now processing the data..." << std::endl;
		if (labels[0] != "" || units[0] != ""){  // case where we need to read a new line before looping
			std::getline(file_in, line0);
		} 
		Nrows=0;
	   while(!file_in.eof()){
		data_str=strsplit(strtrim(line0), " \t"); 
		if (Nrows == 0) {
			data.resize(data_Maxsize, data_str.size());
		}
		for(int i=0; i<data_str.size();i++){
			if ( ! (std::istringstream(data_str[i]) >> tmp_val) ){tmp_val = nan("");} // If the number can be converted, then tmp_val=value. Otherwise tmp_val = NaN
			data(Nrows, i)=tmp_val;		
		}
		if (verbose_data == 1) {std::cout << data.row(Nrows) << std::endl;} // Show all entries only if requested
		std::getline(file_in, line0);
		Nrows=Nrows+1;
	    }
	file_in.close();
	data.conservativeResize(Nrows, data_str.size());
	std::cout << "         - Number of lines found: " << Nrows << std::endl;
	std::cout << "         - Number of columns found: " << data_str.size() << std::endl;
	std::cout << "      ----------------" << std::endl;
     } else {
	msg_handler(cfg_file, "openfile", "Config::read_data_ascii_Ncols()", "Could not open the configuration file!", 1);
     }

     all_data_out.data=data;
     all_data_out.header=header;
     all_data_out.labels=labels;
     all_data_out.units=units;

return all_data_out;
}

std::string Config::format_line(const std::string str){
	/*
	* This program is intended to prepare the string to be splitted
	* and properly formated.
	* [1] Detect if the line has a group indicator.
	*   If it is the case, return the group indicator with the
	*   group indicator symbol
	*   If not,
	* [2] Detect if the line is just a comment line. Comment
	* lines are indicated by a # symbol. 
	*   If not,
	*  	- Remove the comments at the end of a line (indicated 
	*  	  by the symbol ";"
	*  	- Remove blank at the begining/end of the lines
	*   If it is a comment line, returns an empty string

	* 	
	*/
	std::string str_out, str0, char0, terminator;
	size_t pos;

	terminator=";";

	str0=strtrim(str); // remove blanks
	char0=strtrim(str0.substr(0, 1)); 
	if(char0 != "!"){ // If the line has no group indicator
		char0=strtrim(str0.substr(0, 1)); 
		if(char0 != "#" && str0 != ""){ // If the line is not a comment line or an empty line
			pos = str0.find_first_of(terminator); // detect the end of the line
			if(pos != std::string::npos){	
				str_out=str0.substr(0, pos); // remove comments
				str_out=strtrim(str_out); // remove blanks
			} else {
				msg_handler("", "text_terminator", "Config_presets::format_line()", "", 1);
			}
		} else {
			str_out="";
		}
	} else { // If there is a group indicator
		str_out=str0;  // return the line without comments and any blank and without the ":" which terminate
		pos=str0.find_first_of(":"); 
		if(pos != std::string::npos){	
				str_out=str0.substr(0, pos); // remove comments
				str_out=strtrim(str_out); // remove blanks
		} else {
				msg_handler("", "text_terminator", "Config_presets::format_line()", "", 1);
			}
	}

return str_out;
}

void Config::write_cfg_file(std::string cfg_file_out){
/*
 * This function writes the configuration file by reading the configuration
 * within the class Config
*/	

    std::ofstream cfg_file_session;

    cfg_file_session.open(cfg_file_out.c_str());
    if (cfg_file_session.is_open()){
	std::cout << "Configuration File opened... writting lines...";

		cfg_file_session << "#---------------------------------------" << std::endl;
		cfg_file_session << "# Configuration File for the TAMCMC-C++ " << std::endl;
		cfg_file_session << "#   THIS FILE IS WRITTEN AUTOMATICALLY  " << std::endl;
		cfg_file_session << "#       BY Config::write_cfg_file()     " << std::endl;
		cfg_file_session << "#---------------------------------------" << std::endl;


		// The MALA structure
		cfg_file_session << "!MALA:" << std::endl;

		cfg_file_session << "target_acceptance=" << MALA.target_acceptance << ";" << std::endl;
		cfg_file_session << "c0="   << MALA.c0  << ";" << std::endl;
		cfg_file_session << "epsilon1=" << MALA.epsilon1 << ";" <<  std::endl;
		cfg_file_session << "epsilon2=" << MALA.epsi2 << ";" <<  std::endl;
		cfg_file_session << "A1=" << MALA.A1 << ";" << std::endl; 
		for(int i=0; i<MALA.Nt_learn.size(); i++){ cfg_file_session << "Nt_learn="  << MALA.Nt_learn[i] << ";" << std::endl; }
		for(int i=0; i<MALA.Nt_learn.size()-1; i++){ cfg_file_session << "periods_learn=" << MALA.periods_learn[i] << ";" << std::endl;} // MUST BE OF DIM = dim(Nt_learn) -1... 
		cfg_file_session << "use_drift="  << MALA.use_drift << ";" << std::endl; 
		cfg_file_session << "delta=" << MALA.delta << ";" << std::endl; 
		cfg_file_session << "delta_x=" << MALA.delta_x << ";" << std::endl; 
		cfg_file_session << "Nchains=" << MALA.Nchains << ";" << std::endl; 
		cfg_file_session << "dN_mixing=" << MALA.dN_mixing << ";" << std::endl; 
		cfg_file_session << "lambda_temp=" << MALA.lambda_temp << ";" << std::endl; 
		cfg_file_session << "proposal_type=" << MALA.proposal_type << ";" << std::endl; 

		// The modeling structure
		cfg_file_session << "!Modeling:" << std::endl;

		cfg_file_session << "model_fct_name=" << modeling.model_fct_name << ";" << std::endl; 
		cfg_file_session << "prior_fct_name=" << modeling.likelihood_fct_name << ";" << std::endl; 
		cfg_file_session << "likelihood_fct_name=" << modeling.prior_fct_name << ";" << std::endl; 
		cfg_file_session << "likelihood_params=" << modeling.likelihood_params << ";" << std::endl; 
		cfg_file_session << "cfg_model_file=" << modeling.cfg_model_file << ";" << std::endl; 

		// The data structure
		cfg_file_session << "!Data:" << std::endl;

		cfg_file_session << "verbose_data=" << data.verbose_data << ";" << std::endl; 
		cfg_file_session << "file_data=" << data.data_file << ";" << std::endl; 
		cfg_file_session << "x_col=" << data.x_col << ";" << std::endl; 
		cfg_file_session << "y_col=" << data.y_col << ";" << std::endl; 
		cfg_file_session << "ysig_col=" << data.ysig_col << ";" << std::endl; 

		cfg_file_session << "!Outputs:" << std::endl;

		cfg_file_session << "Nsamples=" << outputs.Nsamples << ";" << std::endl; 
		cfg_file_session << "Nbuffer=" << outputs.Nbuffer << ";" << std::endl; 
		cfg_file_session << "erase_old_files=" << outputs.erase_old_files << ";" << std::endl; 
		cfg_file_session << "file_format=" << outputs.file_out_format << ";" << std::endl; 
		cfg_file_session << "get_params=" << outputs.get_params << ";" << std::endl; 
		cfg_file_session << "get_statcriteria=" << outputs.get_statcriteria << ";" << std::endl; 
		cfg_file_session << "get_proposal_params=" << outputs.get_proposal_params << ";" << std::endl; 
		cfg_file_session << "get_parallel_tempering=" << outputs.get_parallel_tempering << ";" << std::endl; 
		cfg_file_session << "get_models=" << outputs.get_models << ";" << std::endl; 
		cfg_file_session << "output_root_name=" << outputs.output_root_name << ";" << std::endl; 
		cfg_file_session << "output_dir=" << outputs.dir_out << ";" << std::endl; 
		cfg_file_session << "params_txt_fileout=" << outputs.params_txt_fileout << ";" << std::endl; 
		cfg_file_session << "proposal_txt_fileout=" << outputs.proposal_txt_fileout << ";" << std::endl; 
		cfg_file_session << "parallel_tempering_txt_fileout=" << outputs.parallel_tempering_txt_fileout << ";" << std::endl; 
		cfg_file_session << "model_txt_fileout= " << outputs.model_txt_fileout << ";" << std::endl; 
		cfg_file_session << "stat_txt_fileout= " << outputs.stat_txt_fileout << ";" << std::endl; 
		cfg_file_session << "acceptance_txt_fileout= " << outputs.acceptance_txt_fileout << ";" << std::endl; 
		//cfg_file_session << "hdf5_fileout=" << outputs.hdf5_fileout << ";" << std::endl; 
		cfg_file_session << "do_restore_last_index=" << outputs.do_restore_last_index << ";" << std::endl; 
		cfg_file_session << "do_restore_proposal=" << outputs.do_restore_proposal << ";" << std::endl; 
		cfg_file_session << "do_restore_variables=" << outputs.do_restore_variables << ";" << std::endl; 
		cfg_file_session << "restore_dir=" << outputs.restore_dir << ";" << std::endl; 
		cfg_file_session << "restore_file_in=" << outputs.restore_file_in << ";" << std::endl; 
		cfg_file_session << "restore_file_out=" << outputs.restore_file_out << ";" << std::endl;

		// The diagnostics structure
		cfg_file_session << "!Diagnostics:" << std::endl;

		cfg_file_session << "chains_diags=" << diags.chains_diags << ";" << std::endl; 
		cfg_file_session << "evidence_diags=" << diags.evidence_diags << ";" << std::endl; 
		cfg_file_session << "pdfs_diags=" << diags.pdfs_diags << ";" << std::endl; 
		cfg_file_session << "model_initial_diags=" << diags.model_initial_diags << ";" << std::endl; 
		cfg_file_session << "model_buffer_diags=" << diags.model_buffer_diags << ";" << std::endl; 
		cfg_file_session << "model_final_diags=" << diags.model_final_diags << ";" << std::endl; 
		cfg_file_session << "output_root_name=" << diags.output_root_name << ";" << std::endl; 
		cfg_file_session << "output_dir=" << diags.output_dir << ";" << std::endl; 
		cfg_file_session << "file_chains_diags=" << diags.file_chains_diags << ";" << std::endl; 
		cfg_file_session << "file_evidence_diags=" << diags.file_evidence_diags << ";" << std::endl; 
		cfg_file_session << "file_pdfs_diags=" << diags.file_pdfs_diags << ";" << std::endl; 
		cfg_file_session << "file_model_init_diags=" << diags.file_model_init_diags << ";" << std::endl; 
		cfg_file_session << "file_model_buffer_diags=" << diags.file_model_buffer_diags << ";" << std::endl; 
		cfg_file_session << "file_model_final_diags=" << diags.file_model_final_diags << ";" << std::endl; 
		cfg_file_session << "data_scoef1=" << diags.data_scoef1 << ";" << std::endl; 
		cfg_file_session << "data_scoef2=" << diags.data_scoef2 << ";" << std::endl; 
		cfg_file_session << "Nclasses=" << diags.Nclasses << ";" << std::endl; 

	std::cout << "Done" << std::endl;
	}
 	cfg_file_session.close();

}

void Config::read_cfg_file(bool verbose){
/*
 * This function reads the configuration file and set all known keywords
 * into the class structures. If a keyword is not recognized, A warning is 
 * issued and the program stops
 * If verbose is 1, then show all the configuration lines. Otherwise shows only the warnings.
*/

    bool keyword_found=0, group_found=0;
    int cpt;
    std::string line0, char0, char_grp; // token
    std::vector<std::string> word;
 
    std::ifstream cfg_file_session;
    cfg_file_session.open(cfg_file.c_str());
    if (cfg_file_session.is_open()) {
	if(verbose ==1){std::cout << "Configuration File opened... processing lines" << std::endl;}

		cpt=0;
		std::getline(cfg_file_session, line0);
		while(!cfg_file_session.eof()){
			line0=format_line(line0);
			if (line0 == "") { 
				std::getline(cfg_file_session, line0); // skip the line
			} 
			if (line0 == "!MALA"){
				if(verbose ==1){std::cout << line0 << std::endl;}
				cpt=cpt+1;
				if(verbose ==1){std::cout << "   [" << cpt << "] " << " Processing the MALA group..." << std::endl;}
				std::getline(cfg_file_session, line0); // get a new line
				char0=strtrim(line0.substr(0, 1)); 
				while(char0 != "!" && line0 != "/END" && !cfg_file_session.eof()){ // Extract parameters for the MALA group
					keyword_found=0;
					line0=format_line(line0);
					word=strsplit(line0, "=");
					if (line0 == ""){
						keyword_found=1; // empty lines are skipped
					}
					if (word[0] == "target_acceptance"){ 
						keyword_found=1;
						MALA.target_acceptance=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      target_acceptance= " << MALA.target_acceptance << std::endl;}
					}
					if (word[0] == "c0"){ 
						keyword_found=1;
						MALA.c0=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      c0= " << MALA.c0 << std::endl;}
					}
					if (word[0] == "epsilon1"){ 
						keyword_found=1;
						MALA.epsilon1=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      epsilon1= " << MALA.epsilon1 << std::endl;}
					}
					if (word[0] == "epsilon2"){ 
						keyword_found=1;
						MALA.epsi2=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      epsilon2= " << MALA.epsi2 << std::endl;}
					}
					if (word[0] == "A1"){ 
						keyword_found=1;
						MALA.A1=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      A1= " << MALA.A1 << std::endl;}
					}
					if (word[0] == "Nt_learn"){ 
						keyword_found=1;
						MALA.Nt_learn=str_to_arrint(word[1], ",");
						if(verbose ==1){
							for (int i=0; i<MALA.Nt_learn.size();i++){
								std::cout << "      Nt_learn[" << i << "]= "<< MALA.Nt_learn[i] << std::endl;
							}
						}
					}
					if (word[0] == "periods_learn"){ 
						keyword_found=1;
						MALA.periods_learn=str_to_arrint(word[1], " ,");
						if(verbose ==1){
							for (int i=0; i<MALA.periods_learn.size();i++){
								std::cout << "      periods_learn[" << i << "]= "<< MALA.periods_learn[i] << std::endl;
							}
						}
					}
					if (word[0] == "use_drift"){ 
						keyword_found=1;
						MALA.use_drift=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      use_drift= " << MALA.use_drift << std::endl;}
					}
					if (word[0] == "delta"){ 
						keyword_found=1;
						MALA.delta=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      delta= " << MALA.delta << std::endl;}
					}
					if (word[0] == "delta_x"){ 
						keyword_found=1;
						MALA.delta_x=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      delta_x= " << MALA.delta_x << std::endl;}
					}
					if (word[0] == "Nchains"){ 
						keyword_found=1;
						MALA.Nchains=str_to_long(word[1]);
						if(verbose ==1){std::cout << "      Nchains= "  << MALA.Nchains << std::endl;}
					}
					if (word[0] == "dN_mixing"){ 
						keyword_found=1;
						MALA.dN_mixing=str_to_int(word[1]);
						if(verbose ==1){std::cout << "      dN_mixing= " << MALA.dN_mixing << std::endl;}
					}
					if (word[0] == "lambda_temp"){ 	
						keyword_found=1;
						MALA.lambda_temp=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      lambda_temp= " << MALA.lambda_temp << std::endl;}
					}
					if (word[0] == "proposal_type"){ 
						keyword_found=1;
						MALA.proposal_type=word[1];
						if(verbose ==1){std::cout << "      proposal_type= " << MALA.proposal_type << std::endl;}
					}
					// After all keywords see whether we detected a known keyword
					if (keyword_found == 0 && char0 != "!" && line0 != "/END" ){
						std::cout << "char0.c_str()= " << char0.c_str() << std::endl;
						std::cout << "Warning: The keyword " << word[0] << " is not a known keyword!" << std::endl;
						std::cout << "Check the syntax of the configuration file " << std::endl;
						std::cout << "The program will exit now" << std::endl;
						exit(EXIT_FAILURE);
					}
					std::getline(cfg_file_session, line0); // get a new line
					char0=strtrim(line0.substr(0, 1)); 
				} // endwhile
			} // endif !MALA
			if (line0 == "!Modeling"){
				cpt=cpt+1;
				if(verbose ==1){std::cout << "   [" << cpt << "] " << " Processing the Modeling group..." << std::endl;}
				std::getline(cfg_file_session, line0); // get a new line
				char0=strtrim(line0.substr(0, 1)); 
				while(char0 != "!" && line0 != "/END" && !cfg_file_session.eof()){ // Extract parameters for the MALA group
					keyword_found=0;
					line0=format_line(line0);
					word=strsplit(line0, "=");
					if (line0 == ""){
						keyword_found=1; // empty lines are skipped
					}
					if (word[0] == "model_fct_name"){ 
							keyword_found=1;
							modeling.model_fct_name=word[1];
							if(verbose ==1){std::cout << "      model_fct_name= " << modeling.model_fct_name << std::endl;}
						}
						if (word[0] == "prior_fct_name"){ 
							keyword_found=1;
							modeling.prior_fct_name=word[1];
							if(verbose ==1){std::cout << "      prior_fct_name= " << modeling.prior_fct_name << std::endl;}
						}

					if (word[0] == "likelihood_fct_name"){ 
						keyword_found=1;
						modeling.likelihood_fct_name=word[1];
						if(verbose ==1){std::cout << "      likelihood_fct_name= " << modeling.likelihood_fct_name << std::endl;}
					}
					if (word[0] == "likelihood_params"){ 
						keyword_found=1;
						modeling.likelihood_params=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      likelihood_params= " << modeling.likelihood_params << std::endl;}
					}
					if (word[0] == "cfg_model_file"){ 
						keyword_found=1;
						modeling.cfg_model_file=word[1];
						if(verbose ==1){std::cout << "      cfg_model_file= " << modeling. cfg_model_file << std::endl;}
					}
					// After all keywords see whether we detected a known keyword
					if (keyword_found == 0 && char0 != "!" && line0 != "/END" ){
						std::cout << "Warning: The keyword " << word[0] << " is not a known keyword!" << std::endl;
						std::cout << "Check the syntax of the configuration file " << std::endl;
						std::cout << "The program will exit now" << std::endl;
						exit(EXIT_FAILURE);
					}
					std::getline(cfg_file_session, line0); // get a new line
					char0=strtrim(line0.substr(0, 1)); 
				} // endwhile
			} // endif !Modeling
			if (line0 == "!Data"){
				cpt=cpt+1;
				if(verbose ==1){std::cout << "   [" << cpt << "] " << " Processing the Data group..." << std::endl;}
				std::getline(cfg_file_session, line0); // get a new line
				char0=strtrim(line0.substr(0, 1)); 
				while(char0 != "!" && line0 != "/END" && !cfg_file_session.eof()){ // Extract parameters for the MALA group
					keyword_found=0;
					line0=format_line(line0);
					word=strsplit(line0, "=");
					if (line0 == ""){
						keyword_found=1; // empty lines are skipped
					}
					if (word[0] == "verbose_data"){ 
						keyword_found=1;
						data.verbose_data=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      verbose_data= " << data.verbose_data << std::endl;}
					}
					if (word[0] == "file_data"){ 
						keyword_found=1;
						data.data_file=word[1];
						if(verbose ==1){std::cout << "      data_file= "<< data.data_file << std::endl;}
					}
					if (word[0] == "x_col"){ 
						keyword_found=1;
						data.x_col=str_to_int(word[1]);
						if(verbose ==1){std::cout << "      x_col= "<< data.x_col << std::endl;}
					}
					if (word[0] == "y_col"){ 
						keyword_found=1;
						data.y_col=str_to_int(word[1]);
						if(verbose ==1){std::cout << "      y_col= "<< data.y_col << std::endl;}
					}
					if (word[0] == "ysig_col"){ 
						keyword_found=1;
						data.ysig_col=str_to_int(word[1]);
						if(verbose ==1){std::cout << "      ysig_col= "<< data.ysig_col << std::endl;}
					}

					// After all keywords see whether we detected a known keyword
					if (keyword_found == 0 && char0 != "!" && line0 != "/END" ){
						std::cout << "Warning: The keyword " << word[0] << " is not a known keyword!" << std::endl;
						std::cout << "Check the syntax of the configuration file " << std::endl;
						std::cout << "The program will exit now" << std::endl;
						exit(EXIT_FAILURE);
					}
					std::getline(cfg_file_session, line0); // get a new line
					char0=strtrim(line0.substr(0, 1)); 
				} // endwhile
			} // endif !Data			
			if (line0 == "!Outputs"){
				cpt=cpt+1;
				if(verbose ==1){std::cout << "   [" << cpt << "] " << " Processing the Outputs group..." << std::endl;}
				std::getline(cfg_file_session, line0); // get a new line
				char0=strtrim(line0.substr(0, 1)); 
				while(char0 != "!" && line0 != "/END"  && !cfg_file_session.eof()){ // Extract parameters for the MALA group
					keyword_found=0;
					line0=format_line(line0);
					word=strsplit(line0, "=");
					if (line0 == ""){
						keyword_found=1; // empty lines are skipped
					}
					if (word[0] == "Nsamples"){ 
						keyword_found=1;
						outputs.Nsamples=str_to_long(word[1]);
						if(verbose ==1){std::cout << "      Nsamples= " << outputs.Nsamples << std::endl;}
					}
					if (word[0] == "Nbuffer"){ 
						keyword_found=1;
						outputs.Nbuffer=str_to_long(word[1]);
						if(verbose ==1){std::cout << "      Nbuffer= " << outputs.Nbuffer << std::endl;}
					}
					if (word[0] == "erase_old_files"){ 
						keyword_found=1;
						outputs.erase_old_files=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      erase_old_files= " << outputs.erase_old_files << std::endl;}
					}
					if (word[0] == "file_format"){ 
						keyword_found=1;
						outputs.file_out_format=strtrim(word[1]);
						if(verbose ==1){std::cout << "      file_format= " << outputs.file_out_format << std::endl;}
					}
					if (word[0] == "get_params"){ 
						keyword_found=1;
						outputs.get_params=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      get_params= " << outputs.get_params << std::endl;}
					}
					if (word[0] == "get_statcriteria"){ 
						keyword_found=1;
						outputs.get_statcriteria=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      statcriteria=" << outputs.get_statcriteria << std::endl;}
					}
					if (word[0] == "get_proposal_params"){ 
						keyword_found=1;
						outputs.get_proposal_params=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      get_proposal_params= " << outputs.get_proposal_params << std::endl;}
					}
					if (word[0] == "get_parallel_tempering"){ 
						keyword_found=1;
						outputs.get_parallel_tempering=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      get_parallel_tempering= " << outputs.get_parallel_tempering << std::endl;}
					}
					if (word[0] == "get_models"){ 
						keyword_found=1;
						outputs.get_models=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      get_models= " << outputs.get_models << std::endl;}
					}
					if (word[0] == "output_root_name"){ 
						keyword_found=1;
						outputs.output_root_name=strtrim(word[1]);
						if(verbose ==1){std::cout << "      output_root_name= " << outputs.output_root_name << std::endl;}
					}
					if (word[0] == "output_dir"){ 
						keyword_found=1;
						outputs.dir_out=strtrim(word[1]);
						if(verbose ==1){std::cout << "      output_dir= " << outputs.dir_out << std::endl;}
					}
					if (word[0] == "params_txt_fileout"){ 
						keyword_found=1;
						outputs.params_txt_fileout=strtrim(word[1]);
						if(verbose ==1){std::cout << "      params_txt_fileout= " << outputs.params_txt_fileout << std::endl;}
					}
					if (word[0] == "proposal_txt_fileout"){ 
						keyword_found=1;
						outputs.proposal_txt_fileout=strtrim(word[1]);
						if(verbose ==1){std::cout << "      proposal_txt_fileout= " << outputs.proposal_txt_fileout << std::endl;}
					}
					if (word[0] == "parallel_tempering_txt_fileout"){ 
						keyword_found=1;
						outputs.parallel_tempering_txt_fileout=strtrim(word[1]);
						if(verbose ==1){std::cout << "      parallel_tempering_txt_fileout= " << outputs.parallel_tempering_txt_fileout << std::endl;}
					}
					if (word[0] == "model_txt_fileout"){ 
						keyword_found=1;
						outputs.model_txt_fileout=strtrim(word[1]);
						if(verbose ==1){std::cout << "      model_txt_fileout= " << outputs.model_txt_fileout<< std::endl;}
					}
					if (word[0] == "stat_txt_fileout"){ 
						keyword_found=1;
						outputs.stat_txt_fileout=strtrim(word[1]);
						if(verbose ==1){std::cout << "      stat_txt_fileout= " << outputs.stat_txt_fileout << std::endl;}
					}
					if (word[0] == "acceptance_txt_fileout"){
						keyword_found=1;
						outputs.acceptance_txt_fileout=strtrim(word[1]);
						if(verbose ==1){std::cout << "      acceptance_txt_fileout= " << outputs.acceptance_txt_fileout << std::endl;}
					}
					/*if (word[0] == "hdf5_fileout"){  // REMOVED ON 23/11/2016 ==> no plan to incorporate hdf5 at short term and probably useless even then
						keyword_found=1;
						outputs.hdf5_fileout=word[1];
						if(verbose ==1){std::cout << "      hdf5_fileout= " << outputs.hdf5_fileout << std::endl;}
					}
					*/
					if (word[0] == "do_restore_proposal"){ 
						keyword_found=1;
						outputs.do_restore_proposal=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      do_restore_proposal= " << outputs.do_restore_proposal << std::endl;}
					}
					if (word[0] == "do_restore_proposal_mean"){ 
						keyword_found=1;
						outputs.do_restore_proposal_mean=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      do_restore_proposal_mean= " << outputs.do_restore_proposal_mean << std::endl;}
					}
					if (word[0] == "do_restore_variables"){ 
						keyword_found=1;
						outputs.do_restore_variables=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      do_restore_variables= " << outputs.do_restore_variables << std::endl;}
					}
					if (word[0] == "do_restore_last_index"){ 
						keyword_found=1;
						outputs.do_restore_last_index=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      do_restore_last_index= " << outputs.do_restore_last_index << std::endl;}
					}
					if (word[0] == "restore_dir"){ 
						keyword_found=1;
						outputs.restore_dir=strtrim(word[1]);
						if(verbose ==1){std::cout << "      restore_dir= " << outputs.restore_dir << std::endl;}
					}
					if (word[0] == "restore_file_in"){ 
						keyword_found=1;
						outputs.restore_file_in=strtrim(word[1]);
						if(verbose ==1){std::cout << "      restore_file_in= " << outputs.restore_file_in << std::endl;}
					}
					if (word[0] == "restore_file_out"){ 
						keyword_found=1;
						outputs.restore_file_out=strtrim(word[1]);
						if(verbose ==1){std::cout << "      restore_file_out= " << outputs.restore_file_out << std::endl;}
					}
					if (word[0] == "do_backup_input_file"){ 
						keyword_found=1;
						outputs.do_backup_input_file=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      do_backup_input_file= " << outputs.do_backup_input_file << std::endl;}
					}
					if (word[0] == "do_backup_cfg_files"){ 
						keyword_found=1;
						outputs.do_backup_cfg_files=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      do_backup_cfg_files= " << outputs.do_backup_cfg_files << std::endl;}
					}
					// After all keywords see whether we detected a known keyword
					if (keyword_found == 0 && char0 != "!" && line0 != "/END" ){
						std::cout << "Warning: The keyword " << word[0] << " is not a known keyword!" << std::endl;
						std::cout << "Check the syntax of the configuration file " << std::endl;
						std::cout << "The program will exit now" << std::endl;
						exit(EXIT_FAILURE);
					}
					std::getline(cfg_file_session, line0); // get a new line
					char0=strtrim(line0.substr(0, 1)); 
				} // endwhile
			} // endif !Ouputs
			if (line0 == "!Diagnostics"){
				cpt=cpt+1;
				if(verbose ==1){std::cout << "   [" << cpt << "] " << " Processing the Diagnostics group..." << std::endl;}
				std::getline(cfg_file_session, line0); // get a new line
				char0=strtrim(line0.substr(0, 1)); 
				while(char0 != "!" && line0 != "/END"  && !cfg_file_session.eof()){ // Extract parameters for the MALA group
					keyword_found=0;
					line0=format_line(line0);
					word=strsplit(line0, "=");
					if (line0 == ""){
						keyword_found=1; // empty lines are skipped
					}
					if (word[0] == "chains_diags"){ 
						keyword_found=1;
						diags.chains_diags=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      chains_diags= " << diags.chains_diags << std::endl;}
					}
					if (word[0] == "pdfs_diags"){ 
						keyword_found=1;
						diags.pdfs_diags=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      pdfs_diags= " << diags.pdfs_diags << std::endl;}
					}
					if (word[0] == "evidence_diags"){ 
						keyword_found=1;
						diags.evidence_diags=str_to_bool(word[1]);
						if(verbose ==1){std::cout << "      evidence_diags= " << diags.evidence_diags << std::endl;}
					}
					if (word[0] == "output_root_name"){ 
						keyword_found=1;
						diags.output_root_name=strtrim(word[1]);
						if(verbose ==1){std::cout << "      output_root_name= " << diags.output_root_name << std::endl;}
					}
					if (word[0] == "output_dir"){ 
						keyword_found=1;
						diags.output_dir=strtrim(word[1]);
						if(verbose ==1){std::cout << "      output_dir= " << diags.output_dir << std::endl;}
					}
					if (word[0] == "file_chains_diags"){ 
						keyword_found=1;
						diags.file_chains_diags=strtrim(word[1]);
						if(verbose ==1){std::cout << "      file_chains_diags= " << diags.file_chains_diags << std::endl;}
					}
					if (word[0] == "file_evidence_diags"){ 
						keyword_found=1;
						diags.file_evidence_diags=strtrim(word[1]);
						if(verbose ==1){std::cout << "      file_evidence_diags= " << diags.file_evidence_diags << std::endl;}
					}
					if (word[0] == "file_pdfs_diags"){ 
						keyword_found=1;
						diags.file_pdfs_diags=strtrim(word[1]);
						if(verbose ==1){std::cout << "      file_pdfs_diags= " << diags.file_pdfs_diags << std::endl;}
					}
					if (word[0] == "model_initial_diags"){ 
						keyword_found=1;
						diags.model_initial_diags=str_to_bool(word[1]);;
						if(verbose ==1){std::cout << "      model_initial_diags= " << diags.model_initial_diags << std::endl;}
					}
					if (word[0] == "model_buffer_diags"){ 
						keyword_found=1;
						diags.model_buffer_diags=str_to_bool(word[1]);;
						if(verbose ==1){std::cout << "      model_buffer_diags= " << diags.model_buffer_diags << std::endl;}
					}
					if (word[0] == "model_final_diags"){ 
						keyword_found=1;
						diags.model_final_diags=str_to_bool(word[1]);;
						if(verbose ==1){std::cout << "      model_final_diags= " << diags.model_final_diags << std::endl;}
					}
					if (word[0] == "file_model_init_diags"){ 
						keyword_found=1;
						diags.file_model_init_diags=word[1];
						if(verbose ==1){std::cout << "      file_model_init_diags= " << diags.file_model_init_diags << std::endl;}
					}
					if (word[0] == "file_model_buffer_diags"){ 
						keyword_found=1;
						diags.file_model_buffer_diags=strtrim(word[1]);
						if(verbose ==1){std::cout << "      file_model_buffer_diags= " << diags.file_model_buffer_diags << std::endl;}
					}
					if (word[0] == "file_model_final_diags"){ 
						keyword_found=1;
						diags.file_model_final_diags=strtrim(word[1]);
						if(verbose ==1){std::cout << "      file_model_final_diags= " << diags.file_model_final_diags << std::endl;}
					}
					if (word[0] == "data_scoef1"){ 
						keyword_found=1;
						diags.data_scoef1=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      data_scoef1= " << diags.data_scoef1 << std::endl;}
					}
					if (word[0] == "data_scoef2"){ 
						keyword_found=1;
						diags.data_scoef2=str_to_dbl(word[1]);
						if(verbose ==1){std::cout << "      data_scoef2= " << diags.data_scoef2 << std::endl;}
					}
					if (word[0] == "show_original_data"){ 
						keyword_found=1;
						diags.show_original_data=str_to_bool(word[1]);;
						if(verbose ==1){std::cout << "      show_original_data= " << diags.show_original_data << std::endl;}
					}
					if (word[0] == "Nclasses"){ 
						keyword_found=1;
						diags.Nclasses=str_to_int(word[1]);
						if(verbose ==1){std::cout << "      Nclasses= " << diags.Nclasses << std::endl;}
					}
					if (word[0] == "evidence_interpolation_factor"){ 
						keyword_found=1;
						diags.evidence_interpolation_factor=str_to_int(word[1]);
						if(verbose ==1){std::cout << "      evidence_interpolation_factor= " << diags.evidence_interpolation_factor << std::endl;}
					}
					// After all keywords see whether we detected a known keyword
					if (keyword_found == 0 && char0 != "!" && line0 != "/END" ){
						std::cout << "Warning: The keyword " << word[0] << " is not a known keyword!" << std::endl;
						std::cout << "line0= " << line0 << std::endl;
						std::cout << "Check the syntax of the configuration file " << std::endl;
						std::cout << "The program will exit now" << std::endl;
						exit(EXIT_FAILURE);
					}
					std::getline(cfg_file_session, line0); // get a new line
					char0=strtrim(line0.substr(0, 1)); 
				} // endwhile
			} // endif !Diagnostics
		}
		if(verbose ==1){std::cout << " Configuration file read successfully..." << std::endl << "    ------------- " << std::endl;}
     } else {
	msg_handler(cfg_file, "openfile", "Config::read_cfg_file()", "Could not open the configuration file!", 1);
     }
}

void Config::read_restore_files(){
/* 
 * Fonction used to read the files in the restore directory
*/

    bool verbose_all=0;
    int keyword_found=0;
    int cpt, i;
    std::string line0, char0;
    std::vector<std::string> word;

    std::string filein_1, filein_2, filein_3;
    std::ifstream restore_file_session1, restore_file_session2, restore_file_session3;

    filein_1=outputs.restore_dir + outputs.restore_file_in + "1.dat";
    filein_2=outputs.restore_dir + outputs.restore_file_in + "2.dat";
    filein_3=outputs.restore_dir + outputs.restore_file_in + "3.dat";

    // -------- Reading file 1 (parameters) ----------
    restore_file_session1.open(filein_1.c_str());
    restore_file_session2.open(filein_2.c_str());
    restore_file_session3.open(filein_3.c_str());

   if (restore_file_session1.is_open() && restore_file_session2.is_open() && restore_file_session3.is_open()) {
		std::cout << "Opening Restoring Files..." << std::endl;
	
		cpt=1;
		std::cout << "   [" << cpt << "] " << " Processing File " << filein_1 << " (do_restore_variables =1)..." << std::endl;	
		keyword_found=0;
		while(!restore_file_session1.eof()){
			std::getline(restore_file_session1, line0);
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1)); 
			if(char0 == "!" && !restore_file_session1.eof()){ 
				line0=strtrim(line0);
				word=strsplit(line0, "=");
				if (word[0] == "! Nvars"){ 
						keyword_found=keyword_found + 1;
						restored_vals.Nvars=str_to_int(word[1]);
						std::cout << "      Nvars= " << restored_vals.Nvars << std::endl;
				}
				if (word[0] == "! Nchains"){ 
					keyword_found=keyword_found + 1;
					restored_vals.Nchains=str_to_int(word[1]);
					std::cout << "      Nchains= " << restored_vals.Nchains << std::endl;
				}
			    if (outputs.do_restore_last_index ==1){
					if (word[0] == "! iteration"){ 
						keyword_found=keyword_found + 1;
						restored_vals.iteration=str_to_int(word[1]);
						std::cout << "      iteration= " << restored_vals.iteration << std::endl;
					}
				} else {
			 		if (word[0] == "! iteration"){
						keyword_found=keyword_found + 1;
						restored_vals.iteration=0;
			 			std::cout << "      iteration= " << restored_vals.iteration  << " (forced to 0 because do_restore_last_index = 0) " << std::endl;
			 		}
		     	}
			if (word[0] == "! variable_names"){ 
				//std::cout << line0 << std::endl;			
		  		//for (int ijk=0; word.size(); ijk++){
		  		//	std::cout << word[ijk] << std::endl;
		  		//}
		  		//exit(EXIT_SUCCESS);
				keyword_found=keyword_found + 1;
				restored_vals.variable_names=strsplit(word[1], " ");
				std::cout << "      variable_names= " << restored_vals.variable_names[0] << " ";
				for(i=1; i<restored_vals.Nvars;i++){
					std::cout << restored_vals.variable_names[i] << " ";
				}
				std::cout << std::endl;
			}
			if (word[0] == "! vars"){ 
					keyword_found=keyword_found + 1;
					restored_vals.vars.resize(restored_vals.Nchains, restored_vals.Nvars); // This assumes that MALA.Nchains and restored_vals.Nvars was set
                    std::cout << restored_vals.Nchains << "   " << restored_vals.Nvars << std::endl;
                    std::cout << "      vars=" << std::endl;
					for(int chain=0; chain<restored_vals.Nchains; chain++){ 
						std::getline(restore_file_session1, line0);
						restored_vals.vars.row(chain)=str_to_Xdarr(line0, " \t");
						std::cout << "         [" << chain << "]   " << restored_vals.vars.row(chain) << std::endl;
					}
			}
			if (word[0] == "! vars_mean"){ 
					keyword_found=keyword_found + 1;
					restored_vals.vars_mean.resize(restored_vals.Nchains, restored_vals.Nvars); // This assumes that MALA.Nchains and restored_vals.Nvars was set
					std::cout << "      vars_mean=" << std::endl;
					for(int chain=0; chain<restored_vals.Nchains; chain++){ 
						std::getline(restore_file_session1, line0);
						restored_vals.vars_mean.row(chain)=str_to_Xdarr(line0, " ");
						std::cout << "         [" << chain << "]   " << restored_vals.vars_mean.row(chain) << std::endl;
					}
			}
			
		}
	}
	restore_file_session1.close();
	if (keyword_found != 6){ 
		msg_handler("", "text_missing_keywords", "Config::read_restore_files()", "", 1);
	}
	
	cpt=2;
	std::cout << std::endl << "   [" << cpt << "] " << " Processing File " << filein_2 << " (do_restore_proposal =1)..." << std::endl;
	keyword_found=0;	
	while(!restore_file_session2.eof()){
             std::getline(restore_file_session2, line0);
	     line0=strtrim(line0);
	     char0=strtrim(line0.substr(0, 1)); 
	     if(outputs.do_restore_proposal == 1){
		if(char0 == "!" && !restore_file_session2.eof()){ 
			line0=strtrim(line0);
			word=strsplit(line0, "=");
			if (word[0] == "! sigmas"){ 
					keyword_found=keyword_found + 1;
					restored_vals.sigma=str_to_Xdarr(word[1], " ");
					if(verbose_all == 1){ std::cout << "      sigma= " << restored_vals.sigma.transpose() << std::endl;}
			}
			if (word[0] == "! mus"){ 
					keyword_found=keyword_found + 1;
					restored_vals.mu.resize(restored_vals.Nchains, restored_vals.Nvars); // This assumes that MALA.Nchains and restored_vals.Nvars was set
					if(verbose_all == 1){ std::cout << "      mus=" << std::endl;}
					for(int chain=0; chain<restored_vals.Nchains; chain++){ 
						std::getline(restore_file_session2, line0);
						restored_vals.mu.row(chain)=str_to_Xdarr(line0, " ");
						if(verbose_all == 1){ std::cout << "         [" << chain << "]   " << restored_vals.mu.row(chain) << std::endl;}
					}
			}
			if (word[0] == "! sigmas_mean"){ 
					keyword_found=keyword_found + 1;
					restored_vals.sigma_mean=str_to_Xdarr(word[1], " ");
					if(verbose_all == 1){ std::cout << "      sigma_mean= " << restored_vals.sigma_mean.transpose() << std::endl;}
			}
			if (word[0] == "! mus_mean"){ 
					keyword_found=keyword_found + 1;
					restored_vals.mu_mean.resize(restored_vals.Nchains, restored_vals.Nvars); // This assumes that MALA.Nchains and restored_vals.Nvars was set
					if(verbose_all == 1){ std::cout << "      mus_mean=" << std::endl;}
					for(int chain=0; chain<restored_vals.Nchains; chain++){ 
						std::getline(restore_file_session2, line0);
						restored_vals.mu_mean.row(chain)=str_to_Xdarr(line0, " ");
						if(verbose_all == 1){ std::cout << "         [" << chain << "]   " << restored_vals.mu_mean.row(chain) << std::endl;}
					}
			}
			

		}
	   } // end if do_restore_proposal
	   else{
		keyword_found=4;
	   }
	}
	if(outputs.do_restore_proposal == 0){ 
		std::cout << "do_restore_proposal = 0 ==> The saved sigmas and mus are not loaded" << std::endl;
	}

        restore_file_session2.close();
	if (keyword_found != 4){ 
		msg_handler("", "text_missing_keywords", "Config::read_restore_files()", "", 1);
	}

	cpt=3;
	std::cout << std::endl << "   [" << cpt << "] " << " Processing File " << filein_3 << " (do_restore_proposal =1)..." << std::endl;
	keyword_found=0;	
	while(!restore_file_session3.eof()){
	    std::getline(restore_file_session3, line0);
	    line0=strtrim(line0);
	    char0=strtrim(line0.substr(0, 1)); 
	    if(outputs.do_restore_proposal == 1){
		if(char0 == "!" && !restore_file_session3.eof()){ 
			line0=strtrim(line0);
			word=strsplit(line0, "=");
			if (word[0] == "! covarmats"){ 
					keyword_found=keyword_found + 1;
					restored_vals.covarmat=initialize_3dMatrix(restored_vals.Nchains, restored_vals.Nvars, restored_vals.Nvars);
					
					std::getline(restore_file_session3, line0); // line with indicator of begining of a new matrix (*0)
					line0=strtrim(line0);
					char0=strtrim(line0.substr(0, 2)); 
					if(char0 == "*0"){
						for(int chain=0; chain<restored_vals.Nchains; chain++){ 
							std::getline(restore_file_session3, line0);
							line0=strtrim(line0);
							char0=strtrim(line0.substr(0, 1));
							i=0;
							if(verbose_all == 1){ std::cout << "      covarmat[" << chain << "]" << std::endl;}
							while(char0 != "*" && char0 != "!" && !restore_file_session3.eof()){
								restored_vals.covarmat[chain]->row(i)=str_to_Xdarr(line0, " ");
								if(verbose_all == 1){ std::cout << "          [" << i << "]   " << restored_vals.covarmat[chain]->row(i) << std::endl;}
							
								std::getline(restore_file_session3, line0);
								line0=strtrim(line0);
								char0=strtrim(line0.substr(0, 1));
								i=i+1;
							}
							if(i > restored_vals.Nvars){ 
								msg_handler("", "readfile", "Config::read_restore_files()", "Syntax error while reading the covarmat: The number of lines do not match the number of variables", 1);
							}
							if(char0 == "!" && !restore_file_session3.eof()){
								word=strsplit(line0, "=");
							}
						}
					} else {
						msg_handler("", "readfile", "Config::read_restore_files()", "Syntax error while reading the covarmat: indicator of first matrix not found", 1);
					}
			}
			if (word[0] == "! covarmats_mean"){ 
					keyword_found=keyword_found + 1;
					restored_vals.covarmat_mean=initialize_3dMatrix(restored_vals.Nchains, restored_vals.Nvars, restored_vals.Nvars);
					
					std::getline(restore_file_session3, line0); // line with indicator of begining of a new matrix (*0)
					line0=strtrim(line0);
					char0=strtrim(line0.substr(0, 2)); 
					if(char0 == "*0"){
						for(int chain=0; chain<restored_vals.Nchains; chain++){ 
							std::getline(restore_file_session3, line0);
							line0=strtrim(line0);
							char0=strtrim(line0.substr(0, 1));
							i=0;
							if(verbose_all == 1){ std::cout << "      covarmat_mean[" << chain << "]" << std::endl;}
							while(char0 != "*" && char0 != "!" && !restore_file_session3.eof()){
								restored_vals.covarmat_mean[chain]->row(i)=str_to_Xdarr(line0, " ");
								if(verbose_all == 1){ std::cout << "          [" << i << "]   " << restored_vals.covarmat_mean[chain]->row(i) << std::endl;}
							
								std::getline(restore_file_session3, line0);
								line0=strtrim(line0);
								char0=strtrim(line0.substr(0, 1));
								i=i+1;
							}
							if(i != restored_vals.Nvars){ 
								msg_handler("", "readfile", "Config::read_restore_files()", "Syntax error while reading the covarmat: The number of lines do not match the number of variables", 1);
							}
							if(char0 == "!" && !restore_file_session3.eof()){
								word=strsplit(line0, "=");
							}
						}
					} else {
						msg_handler("", "readfile", "Config::read_restore_files()", "Syntax error while reading the covarmat: indicator of first matrix not found", 1);
					}			
			}
			
	   }
	   }  // end if do_restore_proposal
	   else{
		keyword_found=2;
	   }
	}
	if (outputs.do_restore_proposal == 0){
		std::cout << "do_restore_proposal = 0 ==> The saved covarmat is not loaded" << std::endl;
	}

	restore_file_session3.close();
	if (keyword_found != 2){ 
		msg_handler("", "text_missing_keywords", "Config::read_restore_files()", "", 1);
	}
   } else {
	msg_handler(outputs.restore_dir + outputs.restore_file_in +"*", "openfile", "Config::read_restore_files()", "At least one file is missing. Please check that the 3 restore files are in the Restore directory", 1);
 }

}


void Config::read_inputs_files(){

    bool passed=0; 

	if(modeling.prior_fct_name == "prior_Test_Gaussian"){ // The structure of such a file is quite simple: Comments (#), Params names (!), Params Inputs, Priors names (!), Priors Inputs
		read_inputs_prior_Simple_Matrix();
		passed=1;
	}
	if(modeling.prior_fct_name == "priors_Harvey_Gaussian"){ // The structure of such a file is quite simple: Comments (#), Params names (!), Params Inputs, Priors names (!), Priors Inputs
		read_inputs_prior_Simple_Matrix();
		passed=1;
	}
	if((modeling.prior_fct_name == "io_MS_Global")){
		read_inputs_priors_MS_Global();
		passed=1;
	}
	if((modeling.prior_fct_name == "io_local")){
		read_inputs_priors_local();
		passed=1;
	}
	if((modeling.prior_fct_name == "io_asymptotic")){
		read_inputs_priors_asymptotic();
		passed=1;
	}
	if(passed == 0){
		msg_handler("", "prior_name", "Config::read_inputs_files()", modeling.prior_fct_name, 1);
	}

}

void Config::read_defautlerrors(bool verbose){
/*
 * This function read a 3 column file that contains the setup for the initial covariance matrix.
 * The diagonal terms of the covariance matrix are given by 1/err with err = Col2 * var + Col3
 * 	- Column 1: Name of the variable
 *	- Column 2: Fraction of the input that defines the err vector (value Col2)
 *	- Column 3: Offset value that defines the err vector (value Col3)
*/

	int i;
	double dbl_out;
	std::vector<int> pos;
	std::string line0, char0;
	std::vector<std::string> tmp;
	std::ifstream cfg_session;
	std::vector<double> frac_errors, abs_errors;

    cfg_session.open(errordefault_file.c_str());
    std::cout << " In read default errors " << std::endl;
    if (cfg_session.is_open()) {

	char0="#"; 
	// ------ Skip header lines ------
	while(char0 == "#" && !cfg_session.eof()){ 
		std::getline(cfg_session, line0);
		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1)); 
	}

	while(!cfg_session.eof()){

		tmp=strsplit(line0, " \t");  // split the line assuming spaces (" ") OR tabulations ("\t") for separator
		MALA.var_names_errors.push_back(strtrim(tmp[0])); // take the var_names_errors
		std::stringstream(strtrim(tmp[1])) >> dbl_out;
		frac_errors.push_back( dbl_out );  // take the fractional error
		std::stringstream(strtrim(tmp[2])) >> dbl_out;
		abs_errors.push_back( dbl_out );  // take the offset error

		std::getline(cfg_session, line0);
		line0=strtrim(line0);
	}
	cfg_session.close();
    } else {
	msg_handler(errordefault_file, "openfile", "Config::read_defautlerrors()", "Could not open the input configuration file! Check if the config file is in the Config directory", 1);
    }

    // ------ Convert the std::vector into VectroXd -------
    MALA.fraction_errors.resize(MALA.var_names_errors.size());
    MALA.offset_errors.resize(MALA.var_names_errors.size());

    for(int i=0; i<MALA.var_names_errors.size(); i++){
	MALA.fraction_errors[i]=frac_errors[i];
	MALA.offset_errors[i]=abs_errors[i];
    }

    // ------ Verbose if requested -------
    if (verbose == 1){
	std::cout << " ------------ Setup for the initial covariance matrix  ---------------------" << std::endl;
	std::cout << "  Configuration File: " << errordefault_file << std::endl;
	std::cout << " ----------------------------  Values   ------------------------------------" << std::endl;
	std::cout << "   Variable name		A			B" << std::endl;
	for(int i=0; i<MALA.var_names_errors.size(); i++){
		std::cout << "  " << MALA.var_names_errors[i]  << "	" << MALA.fraction_errors[i] << "		" << MALA.offset_errors[i] << std::endl;
	}
     }

}

void Config::read_inputs_prior_Simple_Matrix(){
/*
 * Contains the internal structure of an input file when inputs and priors are 
 * Simply organized in a Matrix of parameters. Useful format for simple models
 * ie, when the number of parameters is small and without hyper-parameters. 
 * Example of configuration files that can be read with this function: 
 *         - model_Test_Gaussian
 *         - model_Harvey_Gaussian
*/

	int i;
	std::vector<int> pos;
	std::string line0, char0;
	std::vector<std::string> word, tmp;
	std::ifstream cfg_session;

    cfg_session.open(modeling.cfg_model_file.c_str());
    if (cfg_session.is_open()) {

	char0="#"; 
	// ------ Skip header lines and the range information (because a simple fit does the full spectrum) ------
	while((char0 == "#" || char0 == "*") && !cfg_session.eof()){ 
		std::getline(cfg_session, line0);
		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1)); 
	}
	while(!cfg_session.eof()){
		// ---------- Dealing with the inputs names -----------
		if(char0 == "!"){
			modeling.inputs.inputs_names=strsplit(line0.substr(1), " \t");  // Remove the "!" and split what follows assuming spaces (" ") OR tabulations ("\t") for separator
		} else{
			msg_handler(modeling.cfg_model_file, "text_indicator", "Config::read_inputs_prior_Simple_Matrix()", "", 1);
		}

		// ----------- Dealing with the inputs values ------------
		std::getline(cfg_session, line0);
		line0=strtrim(line0);
		modeling.inputs.inputs=arrstr_to_Xdarrdbl(strsplit(line0, " \t"));

		// ---------- Dealing with relax ---------
		std::getline(cfg_session, line0);
		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1)); 
		if(char0 == "!"){
			std::getline(cfg_session,line0);
			modeling.inputs.relax=arrstr_to_Xiarr(strsplit(line0.substr(1), " \t"));  // Remove the "!" and split what follows assuming spaces (" ") OR tabulations ("\t") for separator
		} else{
			msg_handler(modeling.cfg_model_file, "text_indicator", "Config::read_inputs_prior_Simple_Matrix()", "", 1);
		}

		// ---------- Dealing with the priors names ---------
		std::getline(cfg_session, line0); // Normally, we should get the priors names in line0 now
		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1)); 
		if(char0 == "!"){
			modeling.inputs.priors_names=strsplit(line0.substr(1), " \t");  // Remove the "!" and split what follows assuming spaces (" ") OR tabulations ("\t") for separator
		} else{
			msg_handler(modeling.cfg_model_file, "text_indicator", "Config::read_inputs_prior_Simple_Matrix()", "", 1);
		}

		//  ------- Dealing with the priors values -------
		modeling.inputs.priors.resize(4, modeling.inputs.priors_names.size()); // The maximum number of input parameters for the priors is 4
		modeling.inputs.priors.setConstant(-9999); // Put a recognizable values to indicate empty slots
		i=0;
		
		while(!cfg_session.eof()){
			std::getline(cfg_session, line0);
			line0=strtrim(line0);
			if(line0.size() != 0){
				modeling.inputs.priors.row(i)=arrstr_to_Xdarrdbl(strsplit(line0, " \t"));	
				i=i+1;	
			}		
		}

	}
	cfg_session.close();

	// -------- Determining plength using the parameters names ------
	std::cout << "inputs_names = ";
	for(i=0; i<modeling.inputs.inputs_names.size(); i++){
		std::cout << modeling.inputs.inputs_names[i] << "  ";
	}
	std::cout << std::endl << "     -----------" << std::endl;

	tmp=modeling.inputs.inputs_names;
	i=0;
	while(tmp.size() > 0){
		pos=where_str(tmp, tmp[0]);
		for(int j=0; j<pos.size(); j++){ // Remove the detected values to avoid to read them again
			tmp.erase(tmp.begin() + pos[j]);
		}
		modeling.inputs.plength.conservativeResize(i+1);
		modeling.inputs.plength[i]=pos.size();
		i=i+1;
	}

   } else {
	msg_handler(modeling.cfg_model_file, "openfile", "Config::read_inputs_prior_Simple_Matrix()", "Please check that the config file is in the Config directory", 1);
   }
}


struct Data_Basic Config::read_listfiles(const std::string file, const bool verbose){
/*
 Simple function that read a two-column file that contains the case index in the first column and the name of a function in the second
 Comment lines can be put as header

*/

	std::vector<std::string> tmp, case_val_str, strarr;
	std::string line0, char0;
	std::ifstream cfg_session;
	Data_Basic outputs;
	
    cfg_session.open(file.c_str());
    if (cfg_session.is_open()) {

		char0="#"; 
		// ------ Skip header lines ------
		while(char0 == "#" && !cfg_session.eof()){ 
			std::getline(cfg_session, line0);
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1)); 
		}

		while(!cfg_session.eof()){
			// ---------- Dealing with the inputs strings -----------
			tmp=strsplit(line0.substr(0), " \t");
			if (line0.empty() == 0){ // Increment only lines that are not empty
				case_val_str.push_back(tmp[0]); 
				strarr.push_back(tmp[1]);  // Split the line assuming spaces (" ") OR tabulations ("\t") for separator
			}
				std::getline(cfg_session, line0);
				line0=strtrim(line0);
		}
		cfg_session.close();

	// Process the last line		
		tmp=strsplit(line0.substr(0), " \t");
		if (line0.empty() == 0){ // Increment only lines that are not empty
			case_val_str.push_back(tmp[0]); 
			strarr.push_back(tmp[1]);  // Split the line assuming spaces (" ") OR tabulations ("\t") for separator
		}
		
	// Agregate data into the structure
		outputs.vecXi=arrstr_to_Xiarr(case_val_str);
		outputs.strarr=strarr;
		
	} else{
		msg_handler(file, "openfile", "Config::read_listfiles()", "", 1);
	}
	if (verbose == 1){
		std::cout << " ------ Read parameters ----" << std::endl;
		std::cout << " Model name / Case number" << std::endl;
		for (int i=0; i<outputs.vecXi.size(); i++){
			std::cout << outputs.strarr[i] << "     "  << outputs.vecXi[i] << std::endl;
		}
	}
	return outputs;
}


/*
 * Function that show an error message depending on error_type
 * It receives the filename, fonction name specific message values (arguments) and returns an error code
 * A negative error code correspond to a fatal error. A positive
 * code correspond to a warning.
 * Exit from the program is controlled by the boolean fatal
*/
int Config::msg_handler(const std::string file, const std::string error_type, const std::string fct_name, const std::string arguments, const short int fatal){

	bool err_msg=0;

	if (fatal == 1){
		std::cout << "Warning Fatal error :" << std::endl;
	}
	if (file != ""){ std::cout << "          Accessing file: " << file << std::endl;}
	if (fct_name !=""){ std::cout << "          Fonction name: " << fct_name << std::endl;}

	if(error_type == "text_indicator"){
		std::cout << "Expected the '!' indicator for the relax vector, but not found" << std::endl;
		err_msg=1;
	}
	if(error_type == "text_missing_keywords"){
		std::cout << "Syntax error: At least one of the expected keywords were not found " << std::endl;
		std::cout << "Check the restore files" << std::endl;
		err_msg=1;
	}
	if(error_type == "text_incorrect_keywords"){
		std::cout << "           Incorrect number of keywords. Expected number of keywords: " << arguments << std::endl;
		std::cout << "           Probable reason: mispelt keywords ==> unrecognized. Check the syntax of the configuration file " << std::endl;
		err_msg=1;
	}
	if(error_type == "text_terminator"){
		std::cout << "Terminator not found!" << std::endl;
		std::cout << "   - Each uncommented line that not indicate a group must finish by a ; symbol" << std::endl;
		std:: cout << "  - Each group name must finish by a : symbol" << std::endl;
		std::cout << "Comments may follow the terminator symbol in both cases" << std::endl;
		err_msg=1;
	}
	if(error_type =="restore"){
		std::cout << " restore[" << arguments << "] is not allowed" << std::endl;
		std::cout << " restore[i] must be a number not greater than 3" << std::endl;
		err_msg=1;
	}
        if(error_type == "openfile"){
        	std::cout << arguments << std::endl;
        	std::cout << "Check that the destination path exists" << std::endl;
		err_msg=1;
    	}
	if(error_type == "readfile"){
		std::cout << arguments << std::endl;
		std::cout << "Check the syntax of you restore file" << std::endl;
	}
	if(error_type == "model_name"){
		std::cout << "         Problem in the definition of the model name" << std::endl;
		std::cout << "         Unknown model name detected: " << arguments << std::endl;
		std::cout << "         The name of priors MUST BE those on the following list:" << std::endl;
		for(int i=0; i< modeling.models_case_list_ctrl.size();i++){
			std::cout <<  "              - " << modeling.models_list_ctrl[i] << std::endl;
		}
		std::cout << "         Check that this the names that you use in the configuration file" << std::endl;
	}
	if(error_type == "prior_fctname"){
			std::cout << "         Problem in the definition of the prior names" << std::endl;
			std::cout << "         Unknown prior name detected: " << arguments << std::endl;
			std::cout << "         The name of priors MUST BE those on the following list:" << std::endl;
				for(int i=0; i< modeling.primepriors_case_list_ctrl.size();i++){
					std::cout <<  "              - " << modeling.primepriors_list_ctrl[i] << std::endl;
				}
	}
	if(error_type == "prior_name"){
		std::cout << "         Problem in the definition of the prior name" << std::endl;
		std::cout << "         Unknown model name detected: " << arguments << std::endl;
		std::cout << "         The name of priors MUST BE those on the following list:" << std::endl;
		for(int i=0; i< modeling.priors_case_list_ctrl.size();i++){
			std::cout <<  "              - " << modeling.priors_list_ctrl[i] << std::endl;
		}
		std::cout << "         Check that this the names that you use in the configuration file" << std::endl;
	}
    	if(error_type == "" || err_msg == 0){
    	    std::cout << "Unknown Fatal Error in " << fct_name << std::endl;
    	    std::cout << "Could not open the file: " << file << std::endl;
    	    std::cout << "Neeed debuging..." << std::endl;
    	}
	if(error_type =="likelihood_name"){
		std::cout << "Warning: Error in config.cpp" << std::endl;
		std::cout << "Warning: Problem in the definition of the likelihood name" << std::endl;
		std::cout << "         Unknown model name detected: " << arguments << std::endl;
		std::cout << "         The name of priors MUST BE those on the following list:" << std::endl;
		for(int i=0; i< modeling.likelihoods_case_list_ctrl.size();i++){
			std::cout <<  "              - " << modeling.likelihoods_list_ctrl[i] << std::endl;
		}
		std::cout << "         Check that this the names that you use in the configuration file" << std::endl;
	}
	if(fatal == 1){
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
		return -1;
	} else{
		if(fatal == 2){
		std::cout << "End" << std::endl;
		exit(EXIT_SUCCESS);
		} else{
			return 1;
		}
	}
	
}

