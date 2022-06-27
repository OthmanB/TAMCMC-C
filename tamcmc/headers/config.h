/*
 * Config.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "data.h" // contains the structure Data
#include "matrices.h"
//#include "string_handler.h"
#include "io_ms_global.h"
#include "io_local.h"
#include "io_asymptotic.h"
#include "io_ajfit.h"

//using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

class Config{
/* 
 * The full configuration setup should be handled here. 
*/
	std::string working_dir; 

	public:
		std::string cfg_file; // The configuration file for everything
		std::string errordefault_file; // The configuration file that contains the error values used to initialize the covariance matrix
		std::string cfg_models_ctrl_file;
		std::string cfg_priors_ctrl_file;
		std::string cfg_likelihoods_ctrl_file;
		std::string cfg_primepriors_ctrl_file;
		
		struct MALA_cfg{
			// ----------- Controls for the initial covariance matrix ---------
			std::vector< std::string > var_names_errors; //List of variable names that are recognized in order to set the initial covariance matrix
			VectorXd fraction_errors; // percent of the input value used to initialized the covariance matrix. error[i]=fraction_errors[i]*vars[i] + abs_errors[i]
			VectorXd offset_errors; // Zero offset value used to initialized the covariance matrix. error[i]=fraction_errors[i]*vars[i] + offset_errors[i]
			// ----------- Controls for the Adaptive process -----------
			double target_acceptance;
			double c0;
			long double epsilon1;
			long double epsi2;
			//VectorXd epsilon2; // same as epsi2 but in Vector form
			long double A1;
			std::vector<int> Nt_learn;
			std::vector<int> periods_learn;
			// ---------------- Controls for the Langevin scheme  -----------------
			bool use_drift; // Not operational
			double delta; // Not operational
			long double delta_x; // Not operational
			// ----------- Controls for the parallel tempering ---------
			long Nchains;
			int dN_mixing;
			double lambda_temp;			
			// --------- Experimental extra parameters ----------
			std::string proposal_type;
		};
		struct Models_cfg{
			// --------------- Model parameters --------------
			VectorXi models_case_list_ctrl;
			VectorXi priors_case_list_ctrl, primepriors_case_list_ctrl;
			VectorXi likelihoods_case_list_ctrl;
			
			std::vector<std::string> models_list_ctrl; // List of available models
			std::vector<std::string> priors_list_ctrl, primepriors_list_ctrl; // List of available priors
			std::vector<std::string> likelihoods_list_ctrl; // List of available likelihoods
			
			std::string model_fct_name;
			std::string likelihood_fct_name;
			std::string prior_fct_name;

			int model_fct_name_switch; // the names are encoded into integers so that we can use the switch statement. See Config::convert_model_fct_name_to_switch()
			int likelihood_fct_name_switch; // See Config::convert_likelihood_fct_name_to_switch() for further details
			int prior_fct_name_switch; // See Config::convert_prior_fct_name_to_switch()

			double likelihood_params; // Parameters that could define the likelihood. e.g., in the case of a chi(2,2p) statistics (chi22p function), we need p
			std::string cfg_model_file; // Contains the initial guesses and the priors in an ASCII format. To each model_fct_name, a given format is expected
			int slice_ind; // Index of the subset of data that should be analysed (setting e.g. the frequency range when multiple ranges are requested)
			Input_Data inputs;
		};
		struct Data_cfg{ // all the configuration required for the Data
			bool verbose_data; // Should we show them on screen before proceeding?
			std::string data_file; // The file containing all the data
			std::string type_data;
			int x_col;
			int y_col;
			int ysig_col;
			Data_Nd data_all; // this is the raw outputs of the input data file
			Data data; // this is a formated outputs suitable for the MCMC analysis
		};

		struct Outputs_cfg{
			long Nsamples;
			long Nbuffer;

			bool erase_old_files;

			// --------- Variables that decide which outputs we do write ------------
			bool get_params;  //if set to 1, save the parameters of the fitted model. This option might always be set to 1
			bool get_statcriteria; //if set to 1, save the statistical parameters (logLikelihood, logPrior, logPosterior, ...). Recommended to always set it to 1	
			bool get_proposal_params;  //if set to 1, save the proposal parameters (sigmas, mus, covarmats) of the parallel chains
			bool get_parallel_tempering; //if set to 1, save the parallel tempering information (swaped chain, probability, ...)
			bool get_models;  //if set to 1, save the models for all chains. WARNING: Might greatly slow-down the process and will take a HUGE amount of space

			//  -----------  Variables defining the output files. Do not write the filename full path and extension -----------
			std::string output_root_name;  // The root of the names of the outputs file(e.g. the star identifier such as the KIC number)
			std::string dir_out;  // Output directory for all outputs
			std::string file_out_format; // PLAIN ASCII or BINARY or HDF5 (NOT IMPLEMENTED YET)

			std::string params_txt_fileout;  // extension name for the txt file in output for the parameters
			std::string proposal_txt_fileout; // extension name for the txt file in output for the proposal parameters
			std::string parallel_tempering_txt_fileout; // extension name for the txt file in output for the parallel tempering parameters
			std::string model_txt_fileout; // extension name for the txt file in output for the models
			std::string stat_txt_fileout; // extension name for the txt file in output for the statistical info such as the logLikelihood, logPrior and logPosterior
			std::string acceptance_txt_fileout; // extension name for the txt file with the acceptance rates

			//std::string hdf5_fileout; // extension name for the hdf5 file. All data are put in one single HDF5 if this option is chosen
		
			// ------------ Information required for the completion of an aborted run OR for pursuing a complete run  -------------
			bool do_restore_proposal;
			bool do_restore_variables;
			bool do_restore_last_index;
			bool do_restore_proposal_mean; // Instead of restarting using last position, use the mean ==> Use in principle only for the Acquire phase
			
			std::string restore_dir;
			std::string restore_file_in;
			std::string restore_file_out;

			// ----- Extra backup options ------
			bool do_backup_cfg_files;
			bool do_backup_input_file;
		};
		struct Diagnostics_cfg{ // all the configuration required for the Diagnostics variables
			// --------- Variables that decide which outputs we do plot ------------
			bool chains_diags;
			bool evidence_diags; // If set to 1, plot the evidence function in third page of the acceptance plot and return a file with the values of the evidence at each iteration
			bool pdfs_diags;
	
			std::string output_dir;
			std::string output_root_name;
			std::string file_chains_diags;
			std::string file_evidence_diags;
			std::string file_pdfs_diags;

			bool model_initial_diags;
			bool model_buffer_diags;
			bool model_final_diags;

			std::string file_model_init_diags;
			std::string file_model_buffer_diags;
			std::string file_model_final_diags;

			bool show_original_data;
			int evidence_interpolation_factor; // By how much we increase the number of data points when computing the L_beta
			double data_scoef1; 
			double data_scoef2;
			int Nclasses;
		};
		struct Restored_values{ // Contains the values that were restored but that do not belong to this class (variables, sigma, mu and the covarmat)
			std::vector<std::string> variable_names;
			long iteration; // by default it is 0. If we restore from a previous run, this is changed to the restore file value
			int Nchains;
			int Nvars;
			MatrixXd vars;
			VectorXd sigma;
			MatrixXd mu;
			MatrixXd **covarmat;
			
			VectorXd sigma_mean;
			MatrixXd mu_mean;
			MatrixXd **covarmat_mean;
			MatrixXd vars_mean;
		};

		Data_cfg data;
		MALA_cfg MALA;
		Outputs_cfg outputs;
		Models_cfg modeling;
		Diagnostics_cfg diags;
		Restored_values restored_vals;

		Config(std::string current_path, std::string cfg_file_in, std::string cfg_file_errors,
			   std::string cfg_models_ctrl_file_in, std::string cfg_priors_ctrl_file_in, std::string cfg_likelihoods_ctrl_file_in,
			   std::string cfg_primepriors_ctrl_file_in); // The constructor

		Config(); // Empty constructor if one wants to only use the internal functions
		void init(std::string current_path, std::string cfg_file_in); // An initialization method
		void setup(const int slice_ind);
		void reset(); // reset to 0-sized containers

		void read_defautlerrors(bool verbose); // Function in charge of reading the default error configuration file
		void read_cfg_file(bool verbose);
		void write_cfg_file(std::string cfg_file_out);
		struct Data_Basic read_listfiles(const std::string file, const bool verbose);	// NEW
		VectorXi convert_priors_names_to_switch(const std::vector<std::string> p_names);
		int convert_model_fct_name_to_switch(const std::string model_name);
		int convert_model_fct_name_to_switch(const std::string model_name, const Data_Basic models_ctrl); // designed for stand-alone execution (e.g. for the getmodel tool)
		int convert_prior_fct_name_to_switch(const std::string prior_name);
		int convert_likelihood_fct_name_to_switch(const std::string likelihood_name);

		std::string get_model_fct_name_to_switch(const std::string model_name); // same as convert_model_fct_name_to_switch() but return a string. Can list all models as well (used by getmodel tool)

		inline bool file_exists(const std::string& name); 
		Data_Nd read_data_ascii_Ncols(const std::string file_in_name, const std::string delimiter, const bool verbose_data); // The main function to read ASCII files
		std::string format_line(const std::string str);
		void read_restore_files();
		void read_inputs_files(); // Function that is in charge of handling the input files. These contain the initial values for the parameters as well as the priors
		void read_inputs_prior_Simple_Matrix(); // Procedure that reads simple configuration file, organized as a Matrix of information. See Config.cpp for further details
		void read_inputs_priors_MS_Global(); // For reading a MCMC and setting the modeling.inputs structure of global models
		void read_inputs_priors_local(); // For reading a MCMC and setting the modeling inputs structure of local models
		void read_inputs_priors_asymptotic(); // For read a MCMC and setting the modeling inputs structure for asymptotic models (mixed modes)
		void read_inputs_ajfit(); // Read the inputs in the ajfit format
		int msg_handler(const std::string file, const std::string error_type, const std::string fct_name, const std::string arguments, const short int fatal);
};

