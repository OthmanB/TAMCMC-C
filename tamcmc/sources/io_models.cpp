/*
 * io_models.cpp
 *
 *  This file contains the core instructions that allows you
 *  to generate the correct input vectors for the models
 *
 *  Created on: 16 May 2018
 */

#include <Eigen/Dense>
#include "data.h" // contains the structure Data
//#include "string_handler.h"
#include "io_models.h"
 
IO_models::IO_models(){ // Empty constructor

}

IO_models::~IO_models(){ // Destructor

}

//int IO_models::size(*data){
//	return (*data).priors.rows()
//}

short int IO_models::fill_param(Input_Data *data, const std::string input_name, const std::string prior_name, const double in_val, const VectorXd& prior_vals, const int pos, const int i0_prior){
/*
 * 
 * This function stores the minimum-required values of an Input_Data type structure
 * Note that the data block must have been initialised e.g. using IO_models::initialise_param()
 * No internal checks are performed in this function such that it is up to the user to be careful and use it properly
 * The minimal set of values is:
 *         input_name: Name of the input variable
 *         prior_name: Name of the prior for the variable
 *             in_val: Value of the input variable
 *         prior_vals: Values of the parameters of the prior
 *                pos: The position at which the data will be written into the data structure 
 *			i0_prior:  The zero-index for the values of prior_vals. Usually 1 or 0          
 *
*/
	const int Nmax_prior_params=(*data).priors.rows();
	
	(*data).inputs_names[pos]=input_name; 
	(*data).priors_names[pos]=prior_name;
	
	(*data).inputs[pos]=in_val; // The input value is always in the first element.
	if(prior_name == "Fix"){
		(*data).relax[pos]=0;
		for(int k=0; k<Nmax_prior_params; k++){
			(*data).priors(k,pos)=-9999;//prior_vals(k+1);  // The row becomes a col and this is normal
		}
	} else{
		//std::cout << "taking care of relax...";
		(*data).relax[pos]=1;
		//std::cout << "putting priors...";
		for(int k=0; k<Nmax_prior_params; k++){
			//std::cout << k << "...";
			(*data).priors(k,pos)=prior_vals(k+i0_prior);  // The row becomes a col and this is normal
		}
	}
	return 0;
}

short int IO_models::fill_param_vect(Input_Data *data, const std::vector<double> vec_inputs, 
				const std::vector<bool> vec_relax, const std::string input_name, const std::string prior_name, 
				const VectorXd& prior_vals, const int pos, const int i0_prior_if, const int i0_prior_else){ 	
/* 
	This function takes a vector of inputs and fill an initialized Input_Data structure starting from pos
	All fields are filled with the same parameter name and the same prior definition
	            
*/
   int i;
    for(i=0; i<vec_inputs.size(); i++){
 		if(vec_relax[i]){
			fill_param(data, input_name, prior_name, vec_inputs[i],  prior_vals, i+pos, i0_prior_if);
		} else{
			fill_param(data, input_name,  "Fix", vec_inputs[i],  prior_vals, i+pos, i0_prior_else);		
		}
	}
	//show_param(*data);
	return 0;
}

short int IO_models::fill_param_vect2(Input_Data *data, const std::vector<double> vec_inputs, 
				const std::vector<bool> vec_relax, const std::string input_name, const std::string prior_name, 
				const MatrixXd& prior_vals, const int pos, const int i0_prior_if, const int i0_prior_else){ 	
/* 
	This function takes a vector of inputs and fill an initialized Input_Data structure starting from pos
	All fields are filled with the same parameter name and the same prior definition
	            
*/
   int i;
    for(i=0; i<vec_inputs.size(); i++){
    	std::cout << "      [" << i << "]" << std::endl;
    	std::cout << "       " << prior_vals.row(i) << std::endl;
 		if(vec_relax[i]){
			fill_param(data, input_name, prior_name, vec_inputs[i],  prior_vals.row(i), i+pos, i0_prior_if);
		} else{
			fill_param(data, input_name,  "Fix", vec_inputs[i],  prior_vals.row(i), i+pos, i0_prior_else);		
		}
	}
	//show_param(*data);
	return 0;
}

short int IO_models::add_param(Input_Data *data, Input_Data *data_param, const int pos){//, const int Nparams){
/* 
 *
 * This function adds at the position pos0 of data_all a new block of parameters
 * Note that data must have been initialised before using initialise_param()
 * In particular data.plength must be defined
 * Beware that the use should be cautious: pos0 can take any value such that an improper use may overwrite existing useful data
 *
*/
	const int Nmax_prior_params=(*data).priors.rows();
	const int Nparams=(*data_param).inputs.size(); 

	for(int i=0; i<Nparams; i++){
		(*data).inputs_names[pos + i]=(*data_param).inputs_names[i];
		(*data).priors_names[pos +  i]=(*data_param).priors_names[i];
	}
	
	(*data).inputs.segment(pos , Nparams)=(*data_param).inputs;
	(*data).relax.segment(pos , Nparams)=(*data_param).relax;
	(*data).priors.block(0, pos, Nmax_prior_params, Nparams)=(*data_param).priors;

	return 0; // Returns an error code
}

	
short int IO_models::add_extra_priors(Input_Data *data, const VectorXd& extra, const int pos0){
	
	std::cout << "Yet to be written" << std::endl;
	exit(EXIT_SUCCESS);

	return 0;
}

short int IO_models::append_param(Input_Data *data, Input_Data *data_param, const int Nparams){
/* 
 *
 * This function append an existing data set 'data' with 'data_param', a new block of parameters
 * contrary to add_param, not initialization is required, and the data block is always at the end of the data set
 *
*/

	std::cout << "Yet to be written" << std::endl;
	exit(EXIT_SUCCESS);
	return 0; // Returns an error code
}


short int IO_models::initialise_param(Input_Data *data, const int size_vec, const int Nrows, const VectorXi& plength, const int Nextra_priors){
/* 
 * Define sizes of vectors/matrixes and set the provided plength vector. Other parameters are set to safe default values 
*/	
	initialise_param(data, size_vec, Nrows, plength.size(), Nextra_priors);
	(*data).plength=plength;	
	return 0;
}

short int IO_models::initialise_param(Input_Data *data, const int size_vec, const int Nrows, const int Nplength, const VectorXd& extra_priors){
/* 
 * Define sizes of vectors/matrixes and set the provided extra_priors vector. Other parameters are set to safe default values 
*/	
	initialise_param(data, size_vec, Nrows, Nplength, extra_priors.size());
	(*data).extra_priors=extra_priors;	
	return 0;
}


short int IO_models::initialise_param(Input_Data *data, const int size_vec, const int Nrows, const VectorXi& plength, const VectorXd& extra_priors){
/* 
 * Define sizes of vectors/matrixes and set the provided plength and extra_priors vector. Other parameters are set to safe default values 
*/	
	initialise_param(data, size_vec, Nrows, plength.size(), extra_priors.size());
	(*data).plength=plength;
	(*data).extra_priors=extra_priors;	
	return 0;
}

short int IO_models::show_param(const Input_Data& data, const bool show_metadata){
/* 
 *
 * Function that shows the content of a Input_Data structure in a structured maner
 * The function requires the data block
 *
*/

//	std::cout << "-------------------------------------" << std::endl;	
	std::cout << std::left << "       Index    " <<  std::setw(45) << "Parameter name" << std::setw(15) << "Input value" << std::setw(8) << "Relax" << std::setw(20) << "Prior name";
	std::cout << std::setw(15) << "Prior parameters" << std::endl;   	
	for(int row=0; row<data.inputs.size(); row++){
			if(row < 10){std::cout << std::left << "       [00" << row << "]    ";}
			if((row >=10) & (row < 100)){std::cout << std::left << "       [0" << row << "]    ";}
			if(row >= 100){std::cout << std::left << "       [" << row << "]    ";}
			
			std::cout << std::setw(45) << data.inputs_names[row] << std::setw(15) <<  data.inputs[row] << std::setw(8)  << data.relax[row] << std::setw(20) << data.priors_names[row];
			for(int col=0; col<data.priors.rows(); col++){
				std::cout << std::setw(15) << data.priors(col, row);
			}
			std::cout << std::endl;
		}

	if(show_metadata){
		std::cout << "-------------------------------------" << std::endl;	
		std::cout << "          Extra parameters           " << std::endl;
		for (int i=0; i<data.extra_priors.size(); i++){
			if(i < 10){std::cout << std::left << "       [00" << i << "]    ";}
			if((i >=10) & (i < 100)){std::cout << std::left << "       [0" << i << "]    ";}
			if(i >= 100){std::cout << std::left << "       [" << i << "]    ";}
			std::cout << std::setw(15) << data.extra_priors(i) << std::endl;
		}
		std::cout << "-------------------------------------" << std::endl;	
		std::cout << "    Parameters structure pointer     " << std::endl;
		std::cout << "       " << data.plength.transpose() << std::endl;
	}
	return 0;
}

short int IO_models::show_param(const Input_Data& data){
/*
 * An alias of show_param that do not show the metadata
 *
*/
	return show_param(data, 0);
}

short int IO_models::initialise_param(Input_Data * data, const int size_vec, const int Nrows, const int Nplength, const int Nextra_priors){
/* 
 * Defines the size of the vectors/matrixes within the data structure and initialise it
 * to safe values
 * This function does not modify the 'model_fullname' variable
 *
*/

	if(size_vec > 0 && Nrows > 0){
		(*data).inputs_names.resize(size_vec);
    	(*data).priors_names.resize(size_vec);
    	(*data).priors.resize(Nrows,size_vec);
    	(*data).inputs.resize(size_vec); 
    	(*data).relax.resize(size_vec);
		(*data).priors_names_switch.resize(size_vec);
		
		(*data).priors.setConstant(-9999); // Put a recognizable values to indicate empty slots
		(*data).inputs.setZero(); // At the end, nothing should have this dummy value
		(*data).relax.setZero(); // By default, nothing is relaxed.

	} else{
		if ((size_vec > 0 && Nrows <0) || (size_vec < 0 && Nrows >0)){
			msg_handler("Fatal Error in function IO_models::initialise_param. You need both size_vec and Nrows to be positive or negative", -1);
		} // If both values are negative, no need to initialize
	}
	
	for(int i=0; i<size_vec;i++){
		(*data).inputs_names[i]="Empty";
		(*data).priors_names[i]="Fix"; // Default is a Fixed parameter
	}

 
	if(Nplength > 0){
		(*data).plength.resize(Nplength);
		(*data).plength.setConstant(0);
	}
	if(Nextra_priors > 0){
		(*data).extra_priors.resize(Nextra_priors);
		(*data).extra_priors.setConstant(0); // By default, any extra prior is set to 0
	}
	return 0;
} 
 
void IO_models::msg_handler(const std::string msg, const short int severity){


	std::cout << msg << std::endl;
	if(severity < 0){
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} 
	
}
