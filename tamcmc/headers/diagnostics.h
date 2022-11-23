/*
 *  diagnostics.h
 *
 *  Contains the header to handle graphics (using gnuplot) and other diagnostics outputs
 * 
 *  Created on: 20 Apr 2016
 *      Author: obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>
#include "gnuplot-iostream.h"
#ifdef TAMCMC_WITH_GSL
    #include <gsl/gsl_histogram.h>
    #include <gsl/gsl_rng.h>
#endif

#include "config.h"
#include "interpol.h"

using Eigen::VectorXd;

struct Ptempering_out{
	/*
	 * This is a structure that contains most of values that are saved
	 * for the parallel tempering parameters. Tcoef is available in the header (separated from binary files)
	*/
	std::vector<bool> attempt_mixing; // Tells us if we actually tried to mix the chains
	VectorXi chain0s; // Index of the attempted switched chains. by definition chain1 = chain0 + 1	
	VectorXd Pswitchs; // Probability of switching
	std::vector<bool> switcheds; // A vector of booleans: Switched or not?
};

struct Moves_out{
	MatrixXd Pmoves;
	std::vector< std::vector<bool> > moveds; // A matrix of booleans
};

struct Evidence_out{
		VectorXd beta; // Value of beta=1/Tcoefs used to compute the evidence (x-axis) BEFORE INTERPOLATION ==> FROM RAW DATA
		VectorXd L_beta; // Value of Global Likelihoods at each beta (<P(D|M,I)>_beta)  BEFORE INTERPOLATION ==> FROM RAW DATA
		VectorXd beta_interp; // Value of beta=1/Tcoefs used to compute the evidence (x-axis) AFTER QUADRATIC INTERPOLATION
		VectorXd L_beta_interp; // Value of Global Likelihoods at each beta (<P(D|M,I)>_beta) AFTER QUADRATIC INTERPOLATION
		int interpolation_factor; // multiplicative factor that defines the number of points in the interpolated functions
		long double evidence; // Integral of the L_beta which gives P(M|D, I), USES OF THE INTERPOLATED VALUES.
};

class Diagnostics{
	private:
		bool ps;
		int Nchains;
		int Nbuffer;
		int Nsamples;
		int Nclasses; // number of classes for the histograms
                int firstpass; // Set to 1 initially, when the diagnostic is written, becomes 0
    
		// -- Parameter used for computing the evidence --
		int evidence_interpolation_factor;
		// ------------------------------------------------

		std::string file_out_format;
		std::string filename_likelihood;
		std::string filename_acceptance;
		std::string filename_params; // contains the params for each chain
		std::string filename_data; // contains the data
		std::string filename_likelihood_hdr; // Header of binary files
		std::string filename_params_hdr; // Header if binary files

		bool chains_diags;
		bool evidence_diags;
		bool pdfs_diags;
	
		std::string output_dir;
		std::string output_root_name;
		std::string file_chains_diags;
		std::string file_evidence_diags; 
		std::string file_pdfs_diags; // Output file for the eps file with the pdfs

		bool model_initial_diags; 
		bool model_buffer_diags; 
		bool model_final_diags; 
		std::string file_model_init_diags; 
		std::string file_model_buffer_diags; 
		std::string file_model_final_diags; 

		double data_scoef1; 
		double data_scoef2;
		bool show_original_data; // If set to 1, shows the data with the smooth
		std::string file_data_tmp; // a temporary file used to plot the data with gnuplot
		std::string file_data_likelihood_tmp; // a temporary file used to plot the data with gnuplot
		std::string file_hist_tmp; // a temporary file used to plot the data with gnuplot
	public:
		Diagnostics(Config *config);
		Diagnostics();
		~Diagnostics();

		// --- String Handling procedures ----
		//std::string int_to_str(const int value);
		//std::string dbl_to_str(const double ind);
		//std::string strtrim(const std::string& str);
		//std::vector<std::string> strsplit(const std::string str, const std::string delimiters);
		//Eigen::VectorXd str_to_Xdarr(const std::string str, const std::string delimiters);
		//std::vector<double> str_to_dblarr(const std::string str, const std::string delimiters);
		std::string formated_int_to_str(const int ind_param);
		std::vector<double> vectXd_to_vec(const VectorXd& vecXd_in);

		bool get_model_buffer_diags();
		bool get_model_final_diags();

		// ---- Histograms handling -----
		VectorXd smooth(const VectorXd& xin, const VectorXd& in, const double scoef);
		//VectorXd read_txt_output_params(std::string file); // used to read the output parameters, stored in the txt files
		void write_data_tmp(const VectorXd& x, const VectorXd& y, const VectorXd& z); // used as a temporary storage for the plots of the data with gnuplot
		void fwrite_histogram(const Data_Nd samples, const std::string file_out, const size_t Nclasses, const int ind_param);
		Data_Nd read_txt_output_params(const std::string file, const int Nlines, const bool verbose);

		// ---- Plots handling -------
		void gnuplt_chains_diags(const int i, const VectorXd& Tcoefs);
		void gnuplt_model_diags(Data *data_struc, const VectorXd& model, const std::string phase);
		void gnuplt_model_core(Data *data_struc, const VectorXd& model, const std::string phase);

		void gnuplt_pdfs_diags_main(const int i);
		void gnuplt_pdfs_diags(const std::string dir_pdfs_files, const std::vector<std::string> varnames);

		//Ptempering_out ptemp;
	
		// ---- Handling binary outputs -----
		// Function to read *stat_criteria.bin, *params_chain-*.bin, *model*.bin *sigmas.bin and *mus.bin files. Basically anything that has only double in it
		Eigen::MatrixXd read_bin_matrix_dbl(const std::string binfile, const long Ncols, const long Nrows, const std::string type); 
		// Reading the parallel_tempering outputs
		Ptempering_out read_bin_parallel_temp_params(const std::string binfile, const long Nrows);
		// Reading the Pmoves and moveds variables
		Moves_out read_bin_moves(const std::string binfile, const long Nchains, const long Nrows);

		// Reading the Header of the parameters
		Params_hdr read_params_header(const std::string file);
		// Write a Matrix into an ASCII file
		void write_txt_matrix_dbl(const MatrixXd& Mat, const std::string file_out);

		// Handling the Evidence
		Evidence_out evidence_calc(const VectorXd& Tcoefs, const MatrixXd& Likelihoods);
		void write_evidence(const Evidence_out& evidence, const int i);

		// ------ Execute Shell commands ------
		std::string shell_exec(const std::string cmd);

        	// ------ Error Handler ------
        	int file_error(const std::string file, const std::string error_type, const std::string fct_name);
};
