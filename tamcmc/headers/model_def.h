/*
 * model_def.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "config.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

//using namespace std;

class Model_def{
/*
 * This class contains all the elements (function and parameters) to generate models.
 * The elements containing the parameters can contain as many models as wished. 
 * The number of models in the class is defined by Nmodels.
 * We assume here that all models are generated using the same function. Thus
 * while the parameters of that function may change, the number of parameters 
 * Nparams, Nvars, Ncons must be invariant quantity in function of Nmodels.
*/
		VectorXd cons; // 1D because invariant over Nmodels
		long Ncons, Nvars, Nparams, Nmodels; // Nmodels is the number of models that we have in the class
		std::vector<std::string> vars_names;
		std::vector<std::string> cons_names;
		std::vector<int> index_to_relax; // all positions that correspond to variables into the params VectorXd

		std::string model_fct_name; // Defines which function we for calculating the model
		std::string likelihood_fct_name; // Defines which function we call for calculating the likelihood
		std::string prior_fct_name; // Defines which function we call for calculating the prior on the model

		int model_fct_name_switch; // the names are encoded into integers so that we can use the switch statement. See Config::convert_model_fct_name_to_switch()
		int likelihood_fct_name_switch; // See Config::convert_likelihood_fct_name_to_switch() for further details
		int prior_fct_name_switch; // See Config::convert_prior_fct_name_to_switch()

		std::vector<std::string> likelihood_params_name; // Name of the parameters in likelihood_params
		std::vector<std::string> priors_params_names; // This may have to be included into a structure / class
		std::vector<std::string> params_names;
		VectorXi priors_params_names_switch; // names of the priors were converted into a table of integer (readable by the CASE statement)
		MatrixXd priors_params; // All required parameters for handling priors. This may have to be included into a structure / class of its own.
		//VectorXd likelihood_params; // Any extra parameters required to define the likelihood statistics (e.g. p for a chi(2,2p))
		double likelihood_params; // Any extra parameters required to define the likelihood statistics (e.g. p for a chi(2,2p))
		VectorXi relax; // 1D because invariant over Nmodels
		VectorXi plength; // 1D because invariant over Nmodels
		VectorXd extra_priors;
	public:
		
		MatrixXd params; // could be 2D because could varies over Nmodels
		MatrixXd vars; // could be 2D because could varies over Nmodels
		MatrixXd model; // The model calculated using call_model(). 2D as it varies with Nmodels
		MatrixXd init_model; // INITIAL model calculated using call_model(). 2D as it varies with Nmodels
		VectorXd logLikelihood; // calculated using call_likelihood()
		VectorXd init_logLikelihood; // Added on 7 Dec 2021: Use of the initial loglikelihood in replacement of loglikelihood when logPrior is Infinity ==> Improve performance
		VectorXd logPrior; // calculated using call_prior()
		VectorXd logPosterior; // calculated using generate_model()... This also calls call_model(), call_likelihood() and call_prior()
		VectorXd Pmove; // probability of move... usefull for diagnostics
		std::vector<bool> moved; // 0 if we did not moved, 1 if we moved. Used to define the average acceptance rate
		bool swaped;
		long double Pswap;
		VectorXd comparator_MH; // Uniform number generated for the Metropolis Hasting (FOR DEBUG ONLY)
		long double comparator_PT; // Uniform number generated for the Parallel tempering (FOR DEBUG ONLY)
		//VectorXd Pswap;
		//vector<bool> swaped; // Tells us if we have swaped two chains. The two swaped chains (of index ind_A and ind_A+1) will be at true if swaping happened
		
		Model_def(Config *config, const VectorXd& Tcoefs, const bool verbose); // The constructor
		Model_def(); // Empty constructor (to use only to call functions that do not use internal public/private variables
		~Model_def(); // The destructor		

		VectorXd call_model(Data *data_struc, const int m, bool outparams=false); // call a model using its name. 
		VectorXd call_model_explicit(Data *data_struc, const VectorXi& plength0, const VectorXd& params0, const int model_case, bool outparams=true); // same as call_model(Data *data_struc, int m) but can be called explicitly
		void update_params_with_vars(const long m); // simple loop that refresh values of params(m, *) with vars(m,*)
		inline long double call_likelihood(Data *data_struc, const int m, const VectorXd& Tcoefs); // call the likelihood using its name. Before returning the log(Likelihood), it saves it into 'logLikelihood[m]'
		inline long double call_prior(Data *data_struc, const int m); // call the prior using its name. Before returning the log(prior), we save it into 'logPrior[m]'
		long double generate_model(Data *data_struc, const long m, const VectorXd& Tcoefs); // call successively call_model, call_likelihood and call_prior and then calculates 'model[m,*]', 'logLikelihood[m]', 'logPrior[m]' logPosterior[m]. This is also returned.
};
