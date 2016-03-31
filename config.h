/*
 * Config.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "data.h" // contains the structure Data

using Eigen::VectorXd;
using Eigen::VectorXi;

class Config{
/* 
 * The full configuration setup should be handled here. FOR THE MOMENT I DO NOT USE THIS BUT MUST BE DONE!
*/
	public:
		struct Data_cfg{ // all the configuration required for the Data
			std::string data_file; // The file containing all the data
			std::string type_data;
			Data data;
		};
		struct MALA_cfg{
			string MALA_cfg_file; // The configuration file that contains all the variables for the MALA
			std::string proposal_type;
			std::string model_fct_name;
			std::string likelihood_fct_name;
			std::string prior_fct_name;
			VectorXi relax;
			VectorXi plength;
			VectorXd params_init;
			VectorXd priors_params;
			std::vector<std::string> params_names;
			std::vector<std::string> priors_params_names;
			long double espilon1;
			VectorXd epsilon2;
			long double A1;
			bool use_drift;
			double delta;
			long double delta_x;
			double c0;
			double lambda_temp;
			std::vector<int> Nt_learn;
			std::vector<int> periods_learn;
			int dN_mixing;
			long Nchains;
		};
		struct Outputs_cfg{
			string output_cfg_file; // All the configuration for the outputs
			long Nsamples;
			string output_root_name;  // The root of the names of the outputs file(e.g. the star identifier such as the KIC number)
			string dir_out;  // Output directory for all outputs
			string file_out_format; // ASCII or HDF5
		};
		Data_cfg data;
		MALA_cfg MALA;
		Outputs_cfg outputs;
		Config(string data_file_in, string MALA_cfg_file_in, string output_cfg_file_in); // The constructor
		inline bool file_exists(const std::string& name); 
		std::string strtrim(const std::string& str);
		vector<std::string> strsplit(const std::string str, const std::string delimiter);
		Data_Nd read_data_ascii_Ncols(const string file_in_name, const string delimiter, const bool verbose_data); // The main function to read ASCII files
};

