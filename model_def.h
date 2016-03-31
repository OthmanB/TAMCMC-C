/*
 * model_def.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <string>
#include <vector>
#include "config.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using namespace std;

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
		vector<string> vars_names;
		vector<string> cons_names;
		vector<int> index_to_relax; // all positions that correspond to variables into the params VectorXd
		string model_fct_name; // Defines which function we for calculating the model
		string likelihood_fct_name; // Defines which function we call for calculating the likelihood
		string prior_fct_name; // Defines which function we call for calculating the prior on the model
		vector<string> likelihood_params_name; // Name of the parameters in likelihood_params
		vector<string> priors_params_names; // This may have to be included into a structure / class
		vector<string> params_names;
		VectorXd priors_params; // All required parameters for handling priors. This may have to be included into a structure / class of its own.
		VectorXd likelihood_params; // Any extra parameters required to define the likelihood statistics (e.g. p for a chi(2,2p))
		VectorXi relax; // 1D because invariant over Nmodels
		VectorXi plength; // 1D because invariant over Nmodels
	public:
		
		MatrixXd params; // could be 2D because could varies over Nmodels
		MatrixXd vars; // could be 2D because could varies over Nmodels
		MatrixXd model; // The model calculated using call_model(). 2D as it varies with Nmodels
		VectorXd logLikelihood; // calculated using call_likelihood()
		VectorXd logPrior; // calculated using call_prior()
		VectorXd logPosterior; // calculated using generate_model()... This also calls call_model(), call_likelihood() and call_prior()
		VectorXd Pmove; // probability of move... usefull for diagnostics
		vector<bool> moved; // 0 if we did not moved, 1 if we moved. Used to define the average acceptance rate
		bool swaped;
		long double Pswap;
		//VectorXd Pswap;
		//vector<bool> swaped; // Tells us if we have swaped two chains. The two swaped chains (of index ind_A and ind_A+1) will be at true if swaping happened

		Model_def(string m_name, string l_name, string p_name, VectorXi rlx, VectorXi plgth, VectorXd params_in, 
		     vector<string> params_in_names, VectorXd priors_in, vector<string> priors_nme, long Nmdls, 
		     Config *config, VectorXd Tcoefs); // the config is required to initialize model... when all finished, this should be the only input of the function
		inline VectorXd call_model(Data *data_struc, int m); // call a model using its name. 
		void update_params_with_vars(long m); // simple loop that refresh values of params(m, *) with vars(m,*)
		inline long double call_likelihood(Data *data_struc, int m, VectorXd Tcoefs); // call the likelihood using its name. Before returning the log(Likelihood), it saves it into 'logLikelihood[m]'
		inline long double call_prior(Data *data_struc, int m); // call the prior using its name. Before returning the log(prior), we save it into 'logPrior[m]'
		long double generate_model(Data *data_struc, long m, VectorXd Tcoefs); // call successively call_model, call_likelihood and call_prior and then calculates 'model[m,*]', 'logLikelihood[m]', 'logPrior[m]' logPosterior[m]. This is also returned.
};
