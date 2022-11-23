/*
 * io_models.h
 *
 *  header for the core instructions that allows you
 *  to generate the correct input vectors for the models
 *
 *  Created on: 16 May 2018
 */
 
#pragma once
//#include <string>
//#include <vector>
#include <Eigen/Dense>
#include "data.h" // contains the structure Data
//#include "matrices.h"
//#include "string_handler.h"

//using Eigen::MatrixXd;
using Eigen::VectorXd;
//using Eigen::VectorXi;

class IO_models{		

	public:
		IO_models();
		~IO_models();
		short int initialise_param(Input_Data * data, const int size_vec, const int Nrows, const int Nplength, const int Nextra_priors);
        short int initialise_param(Input_Data *data, const int size_vec, const int Nrows, const VectorXi& plength, const int Nextra_priors);
        short int initialise_param(Input_Data *data, const int size_vec, const int Nrows, const int Nplength, const VectorXd& extra_priors);
        short int initialise_param(Input_Data *data, const int size_vec, const int Nrows, const VectorXi& plength, const VectorXd& extra_priors);
		//Input_Data create_param();
        short int fill_param(Input_Data *data, const std::string input_name, const std::string prior_name, const double in_val, const VectorXd& prior_vals, const int pos, const int i0_prior);
        //short int fill_param(Input_Data *data, const std::string input_name, const std::string prior_name, const double in_val, const double prior_val, const int pos);
		short int fill_param_vect(Input_Data *data, const std::vector<double> vec_inputs, 
				const std::vector<bool> vec_relax, const std::string input_name, const std::string prior_name, 
				const VectorXd& prior_vals, const int pos, const int i0_prior_if, const int i0_prior_else); 	
		short int fill_param_vect2(Input_Data *data, const std::vector<double> vec_inputs, 
				const std::vector<bool> vec_relax, const std::string input_name, const std::string prior_name, 
				const MatrixXd& prior_vals, const int pos, const int i0_prior_if, const int i0_prior_else); // same as fill_param_vect but prior_vals is a Matrix (overload did not work due to ambiguity between VectorXd and MatrixXd)
		short int add_param(Input_Data *all_in, Input_Data *data_param, const int pos0); // Get a structure of parameters (data_param) that we wish to add into the main structure
        short int append_param(Input_Data *all_in, Input_Data *data_param, const int Nparams);
		short int add_extra_priors(Input_Data *all_in, const VectorXd& extra, const int pos0); //
    
    	short int show_param(const Input_Data& data, const bool show_metadata);
    	short int show_param(const Input_Data& data);
    	
		void msg_handler(const std::string msg, const short int severity);
		
};
