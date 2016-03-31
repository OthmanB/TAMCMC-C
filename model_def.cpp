/*
 * model_def.cpp
 *
 * Contains the methods associated to the class 'Model'
 * The header file of that class is in encapsulators.h
 * Note that all models should be put in models.cpp
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include "model_def.h"
#include "models.h"   // This contains the models that can be used
#include "likelihoods.h"  // This contains the likelihoods that can be used
#include "priors_calc.h"  // This contains the priors that can be used
//#include "config.h"
//#include "data.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using namespace std;


Model_def::Model_def(string m_name, string l_name, string p_name, VectorXi rlx, VectorXi plgth, VectorXd params_in, 
		     vector<string> params_in_names, VectorXd priors_in, vector<string> priors_nme, long Nmdls, Config *config, VectorXd Tcoefs){


	Nmodels=Nmdls;
	model_fct_name=m_name;
	likelihood_fct_name=l_name;
	prior_fct_name=p_name;
	params_names=params_in_names;
	priors_params_names=priors_nme;
	relax=rlx;
	plength=plgth;
	
	Nparams=plength.sum();
	params.resize(Nmodels, Nparams);

	for (int m=0; m<Nmodels; m++){
		params.row(m)=params_in;
	}
	priors_params=priors_in;
	
	Pmove.setZero(Nmodels); // initialize the value of the move probabilty for all models
	moved.assign(Nmodels, 0); // initialize to False the move action indicator for all models
	
	//swaped.assign(Nmodels, 0); // initialize to False the swap action indicator for all models
	//Pswap.setZero(Nmodels);
	swaped=0;
	Pswap=0;

	// --- initialise vectors of variables/constant ------
	Nvars=0; // initialize the number of variables to 0
	Ncons=0; // initialize the number of constant to 0
	vars.resize(Nmodels, Nparams); // set to maximum possible size. This must be 2D
	cons.resize(Nparams);
	index_to_relax.resize(Nparams);
	vars_names.resize(Nparams); // set to maximum possible size this is 1D
	cons_names.resize(Nparams);

	if( relax.sum() <= 1){
		std::cout << "Less than one parameter free! No minimization possible!" << std::endl;
		std::cout << "You need to set relax to 1 for at least two parameters" << std::endl;
	}
	for (int i=0; i<Nparams;i++){
		if(relax[i] == 1){
			vars(0,Nvars)=params(0,i);
			vars_names[Nvars]=params_names[i];
			index_to_relax[Nvars]=i; // used to quickly refresh params when calculating all the models
			Nvars=Nvars+1;
		} else{
			cons[Ncons]=params(0,i);
			cons_names[Ncons]=params_names[i];
			Ncons=Ncons+1;
		}
	}
	vars.conservativeResize(Nmodels, Nvars); // Resize by keeping the elements inside
	vars_names.resize(Nvars); // This is a vector... we can do a normal resize
	cons.resize(Ncons); // This is also a vector, not a VectorXd or a MatrixXd
	cons_names.resize(Ncons); // again a vector...
	index_to_relax.resize(Nvars); // vector with position of params that vary

//	Debug only
	std::cout << " --- Debug Only --- " << std::endl;
	std::cout << "Nmodels =" << Nmodels << std::endl;
	std::cout << "Number of variables: " << Nvars << std::endl;
	std::cout << "Number of constants: " << Ncons << std::endl;
	std::cout << "Number of parameters " << plength.sum() << std::endl;
	std::cout << "List of the parameters : " << std::endl;
	cout << "        ------------- " << endl;
	std::cout << "   - Variables " << std::endl;
	for(int i=0; i<Nvars; i++){
		cout << vars_names[i]  << " | "  << vars(0,i) << endl;
	}
	std::cout << "   - Constant " << std::endl;
	for(int i=0; i<Ncons; i++){
		cout << cons_names[i] << " | " << cons[i] << endl;
	}
	cout << " ------------- " << endl;

	// Initialize the model, logLikelihood, logPrior, logPosterior variables
	model.resize(Nmodels, (*config).data.data.Nx);
	logLikelihood.resize(Nmodels);
	logPrior.resize(Nmodels);
	logPosterior.resize(Nmodels);
	
	Data data_in=config->data.data;
	bool empty_container=0; // Used to know whether we fill MatrixXd/VectorXd for all Nmodels
	if(empty_container == 0){
		for(int m=0; m<Nmodels; m++){
			generate_model(&data_in, m, Tcoefs); // No need of the returned value
		} 
	}
//

}

VectorXd Model_def::call_model(Data *data_struc, int m){

	bool passed=0;

	if(model_fct_name == "model_MS_Global"){
		passed=1;
		return model_MS_Global(params.row(m), plength, (*data_struc).x);
	} 
	if(model_fct_name == "model_Test_Gaussian"){
		passed=1;
		/*std::cout << "  In call model " << std::endl;
		std::cout << "  params..." << std::endl;
		std::cout << params.row(m) << std::endl;
		std::cout << "-------" << std::endl;
		std::cout << "  plength..." << std::endl;
		std::cout << plength << std::endl;
		*/
		return model_Test_Gaussian(params.row(m), plength, (*data_struc).x);
	}
	
	// after the full list of possibility, we check that passed is not yet at 0
	if(passed == 0){
		cout << "model_name = " << model_fct_name << " was not recognized" << endl;
		cout << "Please only use know reference keywords for the model_name" << endl;
		cout << "Accepted keywords so far:" << endl;
		cout << "    - MS_Global" << endl;
		cout << "    - Test_Gaussian (For Debug only)" << endl;
		std::cout << "The program will exit now" << std::endl;
		exit ( EXIT_FAILURE );
	}
}

long double Model_def::call_likelihood(Data *data_struc, int m, VectorXd Tcoefs){
/* 
 * call the likelihood using its name. Before returning the log(Likelihood), it saves it into 'logLikelihood'
*/
	bool passed;
	double p;
	long double logL;

	passed=0;
	if(likelihood_fct_name == "chi(2,2p)"){
		passed=1;
		p=likelihood_params[0]; // to construct the chi(2,2p) statistics, we need p
		logL=likelihood_chi22p((*data_struc).y, model.row(m), p); // 'Model_def::model' must be calculated prior to the execution of this function
		return logL/Tcoefs[m];
	}
	if(likelihood_fct_name == "chi_square"){
		passed=1;
		logL=likelihood_chi_square((*data_struc).y, model.row(m), (*data_struc).sigma_y); // Here no need of likelihood_params
		return logL/Tcoefs[m];
	}	
	// after the full list of possibility, we check that passed is not yet at 0
	if(passed == 0){
		cout << "likelihood_name = " << likelihood_fct_name << " was not recognized" << endl;
		cout << "Please only use know reference keywords for the model_name" << endl;
		cout << "Accepted keywords so far:" << endl;
		cout << "    - chi(2,2p)" << endl;
		std::cout << "The program will exit now" << std::endl;
		exit ( EXIT_FAILURE );
	}
	
}

long double Model_def::call_prior(Data *data_struc, int m){
/* 
 * call the prior using its name. Before returning the log(prior), we save it into 'logPrior'
*/
	bool passed=0;
	double p;

	if(prior_fct_name == "prior_MS_Global"){
		passed=1;
		return priors_MS_Global(params, plength, priors_params);
	}
	if(prior_fct_name == "prior_Test_Gaussian"){
		passed=1;
		//std::cout << "In prior Test case " << std::endl;
		return priors_Test_Gaussian(params.row(m), plength, priors_params);
	}
	// after the full list of possibility, we check that passed is not yet at 0
	if(passed == 0){
		cout << "prior_name = " << prior_fct_name << " was not recognized" << endl;
		cout << "Please only use know reference keywords for the model_name" << endl;
		cout << "Accepted keywords so far:" << endl;
		cout << "    - prior_MS_Global" << endl;
		cout << "    - prior_Test_Gaussian (For Debug only)" << endl;
		std::cout << "The program will exit now" << std::endl;
		exit ( EXIT_FAILURE );
	}
}

long double Model_def::generate_model(Data *data_struc, long m, VectorXd Tcoefs){
/*
 * call successively call_model, call_likelihood and call_prior and then calculates the logPosterior. This is also returned.
*/
	//std::cout << "In generate model " << std::endl;
	model.row(m)=call_model(data_struc, m);
	//std::cout << "      - After call model " << std::endl;

	logLikelihood[m]=call_likelihood(data_struc, m, Tcoefs); // logLikelihood saved at element m of the vector
	//std::cout << "      - logLikelihood[" << m << "]=" << logLikelihood[m]  << std::endl;
	logPrior[m]=call_prior(data_struc, m); // logprior saved at element m of the vector
	//std::cout << "      - logPrior[" << m << "]=" << logPrior[m]  << std::endl;
	logPosterior[m]= logLikelihood[m] + logPrior[m]; // logPosterior saved at element m of the vector
	//std::cout << "      - logPosterior[" << m << "]=" << logPosterior[m]  << std::endl;

return logPosterior[m];
}


void Model_def::update_params_with_vars(long m){
	//std::cout << "    In update params" << std::endl;
	//std::cout << "    chain = " << m << std::endl;
	
	/*for(int i=0; i<index_to_relax.size(); i++){
		std::cout << index_to_relax[i] << std::endl;
	}
	std::cout << "---------- " << std::endl;

	for(int i=0; i<params.row(m).size(); i++){
		std::cout << params(m,i) << std::endl;
	}
	std::cout << "---------- " << std::endl;

	for(int i=0; i<vars.row(m).size(); i++){
		std::cout << vars(m,i) << std::endl;
	}
	std::cout << "---------- " << std::endl;
	*/
	for( int i=0; i<index_to_relax.size(); i++){
		//std::cout << "i=" << i << std::endl;
		params(m,index_to_relax[i])=vars(m,i);
	}
	

}


