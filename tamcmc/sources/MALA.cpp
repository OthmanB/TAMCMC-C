/*
 * MALA.cpp
 *
 *  This file contains the core program for the MCMC
 *  as it is defined in Atchade Y.F., 2006, Meth. Comp. In Applied Probab. 8, 235
 *  This MCMC is an adaptive truncated-Langevin Metropolis Hasting algorithm.
 *  The strategy involves a Robins-Monroe stochastic algorithm during the learning
 *  process and optionaly can use the local hyperspace gradient to enhance the sampling
 *  (MALA option). 
 *  The adpative part of the algorithm tries to minimize differences between successive 
 *  samples. The minimized quantities are:
 *      - The first moment (mean). It should be stable and converge towards mu when
 *        Nsamples ---> infinity
 *      - The second moment (covarmat) should be stable and converge towards 
 *        The covariance matrix near the global maxima
 *      - A scaled variance (mean variance). It is here to ease the search
 *        for the best covariance matrix. It should be stable when N ---> infinity
 *  WARNING: The current implementation does not use the Gradient.
 *  Remark: The algorithm proposed by Atchade does not implement parallel tempering
 *          Here, we added the parallel tempering which greatly improve the performance.
 *
 *  Created on: 02 Mar 2016
 */
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "unistd.h"
#include "MALA.h"
//#include <boost/random.hpp>
//#include <boost/random/normal_distribution.hpp>

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;


//MatrixXd** initialize_3dMatrix(int Nchains, int Nrows, int Ncols);
//void set_3dMatrix(MatrixXd** matrix3d, MatrixXd m_in, int position);
MatrixXd** copy_3dMatrix(MatrixXd** m3d_in, int depth);
MatrixXd multivariate_eigen(VectorXd, MatrixXd, int);
double *r8vec_normal_01 (const int n, int *seed );
double *uniform_01 (const int n, int *seed );
extern long ben_clock();

MALA::MALA(Config *cfg){ // Constructor
/* 
 * This procedure loads the initial configuration.
 * The configuration is read from a configuration file.

*/
	Nsamples=cfg->outputs.Nsamples;
	Nchains=cfg->MALA.Nchains;
	
	seed=(unsigned)time(NULL);
	srand(seed);

	// Identify the variables which in turns will enable us to initialize the proposal parameters (e.g. covarmat)
	Nvars=0; // This is a variable of the class... see MALA.h
	VectorXd vars(cfg->modeling.inputs.inputs.size());
	std::vector<std::string> vars_names(cfg->modeling.inputs.inputs.size());
	for (int i=0; i<cfg->modeling.inputs.inputs.size();i++){
		if(cfg->modeling.inputs.relax[i] == 1){
			vars(Nvars)=cfg->modeling.inputs.inputs(i);
			vars_names[Nvars]=cfg->modeling.inputs.inputs_names[i];
			Nvars=Nvars+1;
		}
	}
	vars.conservativeResize(Nvars);
	vars_names.resize(Nvars);
  
	VectorXd tmp(Nvars);

	epsilon1=cfg->MALA.epsilon1;
	tmp.setConstant(cfg->MALA.epsi2);
	epsilon2=tmp.asDiagonal();
	A1=cfg->MALA.A1;
	delta=cfg->MALA.delta;
	delta_x=cfg->MALA.delta_x;
	c0=cfg->MALA.c0;
	lambda_temp=cfg->MALA.lambda_temp;
	use_drift=cfg->MALA.use_drift;
	Nt_learn=cfg->MALA.Nt_learn;
	periods_learn=cfg->MALA.periods_learn;
	dN_mixing=cfg->MALA.dN_mixing;
	target_acceptance=cfg->MALA.target_acceptance;
	gamma=c0; // Forcing initialization of gamma
      
    initial_i=0;
	if(cfg->outputs.do_restore_last_index == 1) {
		initial_i=cfg->restored_vals.iteration;
	} else {
		initial_i=0;
	}
	Tcoefs.resize(Nchains);	
	for (int m=0; m<Nchains; m++){Tcoefs[m]=std::pow(lambda_temp, m);} // initialize the list of temperatures

    if (cfg->outputs.do_restore_proposal == 0) {
        // Import the list of error values as defined in the defaut file used to setup the covariance matrix AND initialise the covariance matrix
        init_proposal(vars, vars_names, cfg->MALA.var_names_errors, cfg->MALA.fraction_errors, cfg->MALA.offset_errors); // Initialize the covariance matrix and all parameters relative to the proposal OR
    } else {
        restore_proposal(vars, cfg); // Restore the proposal from a file
    }

}

MALA::~MALA(){ // Destructor

	Nt_learn.resize(0);
	periods_learn.resize(0);
	epsilon2.resize(0,0);
	sigma.resize(0);
	Tcoefs.resize(0);
	mu.resize(0,0);

	destroy_3dMatrix(covarmat, Nchains);
	covarmat=NULL; // set point to 0
}

void MALA::destroy_3dMatrix(MatrixXd** m3d, const int depth){

	for(int d1=0; d1 < depth; d1++){
		delete m3d[d1];
	}
	delete m3d;
}

long double MALA::p1_fct(long double x){
	
	long double scalar_p;
	if(x >= epsilon1 && x<=A1){
		scalar_p=x;
	}
	if(x < epsilon1){
		scalar_p=epsilon1;
		std::cout << "Warning: Hit epsilon1 in p1, value of x="<< x << std::endl;	
		std::cout << "                               epsilon1=" << epsilon1 << std::endl;
	}
	if(x > A1){
		scalar_p=A1;
		std::cout << "Warning: Hit A1 in p1" << std::endl;
	}
return scalar_p;
}

MatrixXd MALA::p2_fct(const MatrixXd& x){

	MatrixXd matrix_p;

	if(x.norm() <= A1){
		matrix_p=x;
	} else {
		matrix_p=x * A1/x.norm();
		std::cout << "Warning: Hit A1 in p2" << std::endl;
	}
	return matrix_p;
}

VectorXd MALA::p3_fct(const VectorXd& x){
	VectorXd vector_p;
	
	if(x.norm() <=A1){
		vector_p=x;
	} else{
		vector_p=x * A1/x.norm();
		std::cout << "Warning: Hit A1 in p3" << std::endl;
	}
	return vector_p;
}


long MALA::random_int_vals(int N0, int Nmax){
/* 
 * Generate a uniform random number between 1 and Nmax
 * WARNING: srand(time(NULL)) must be put in the main
 *          program in order to generate proper
 *          pseudo-random number!
 */

return N0 + rand()%Nmax; // srand(seed) is executed by the constructor of the MALA class
}


void MALA::restore_proposal(const VectorXd& vars, Config *cfg){
/* 
 * Used to setup the initial values for the proposal law when it is requested
 * to restore the variables (do_restore=1 in the configuration file)
*/
	
	MatrixXd cov0;
	MatrixXd* covtmp;

	std::cout << "     ---------------------------------------- " << std::endl;
	std::cout << "     Nsamples that remains to process = " << Nsamples - initial_i << std::endl;
	std::cout << "     Nvars = " << Nvars << std::endl;
	std::cout << "     Nchains = " << Nchains << std::endl;
	std::cout << "     ---------------------------------------- " << std::endl;
 
	if(cfg->outputs.do_restore_proposal_mean == 0){
		std::cout << "             Use the last values of the user-specified previous analysis..." << std::endl;
		std::cout << " " << std::endl;
		sigma=cfg->restored_vals.sigma;
		mu=cfg->restored_vals.mu;
	} else{
		std::cout << "             Use the averaged values over Nbuffer of the user-specified previous analysis..." << std::endl;
		std::cout << " " << std::endl;
		sigma=cfg->restored_vals.sigma_mean;
		mu=cfg->restored_vals.mu_mean;
	}
	covarmat=initialize_3dMatrix(Nchains, Nvars, Nvars);

	for(int m=0; m<Nchains; m++){
		if(cfg->outputs.do_restore_proposal_mean == 0){
			covtmp=cfg->restored_vals.covarmat[m];
		} else{
			covtmp=cfg->restored_vals.covarmat_mean[m];
		}
		set_3dMatrix(covarmat, *covtmp, m);

		std::cout << "scaling factor sigma[" << m << "]=" << sigma[m] << std::endl;
		std::cout << "Initial Covariance Matrix covarmat["<< m << "]=" << std::endl;
		if(Nvars < 10){
			std::cout << *covarmat[m] << std::endl;
		} else {
			cov0=*covarmat[m];
			std::cout << "Diagonal only is shown due to the large size of the matrix" << std::endl;
			std::cout << cov0.diagonal() << std::endl;
		}
		std::cout << " -------------------------------- "  << std::endl;
		
		std::cout << "Initial Vectors mu[" << m << ",*]=" << mu.row(m) << std::endl;
	}
	
	std::cout << " Vector of initial variables: " << std::endl;
	std::cout << vars.transpose() << std::endl;

}

void MALA::init_proposal(const VectorXd vars, const std::vector<std::string> var_names, 
			 const std::vector<std::string> s_inerror, const VectorXd fracerr, const VectorXd offseterr){

	VectorXd error(Nvars);
	MatrixXd cov0;

	error.setConstant(1); // Defaut non-zero value;
	for(long i=0; i<vars.size(); i++){
		for(int j=0; j<s_inerror.size();j++){
			if(var_names[i] == s_inerror[j]){
				error[i]=vars[i]*fracerr[j] + offseterr[j];
			}
		}
	}
    for(long i=0; i<vars.size(); i++){
        std::cout << " var[" << i << "]=" << vars[i] <<"   ===> err[" << i << "]=" << error[i]  <<  std::endl;
    }

	std::cout << " ---------------------------------------------------------------" << std::endl;
	std::cout << " Initializing the matrix of covariance and the scaling factor..." << std::endl;
        error=error.array().square();
	cov0=error.asDiagonal();
	covarmat=initialize_3dMatrix(Nchains, Nvars, Nvars);
	sigma.setZero(Nchains);
	mu.setZero(Nchains, Nvars);
	for(int m=0; m<Nchains; m++){
		set_3dMatrix(covarmat, cov0, m); // the error variance is the initial matrix
		sigma[m]=std::pow(2.38,2)*std::pow(Tcoefs[m],0.2)/Nvars; // If gaussian pdfs, this should ensure proper scaling (no security ==> better initial acceptance but more risk of hitting the bundary p1)
		std::cout << "scaling factor sigma[" << m << "]=" << sigma[m] << std::endl;
		std::cout << "Initial Covariance Matrix covarmat["<< m << "]=" << std::endl;
		//std::cout << "        " << *covarmat[m] << std::endl;
		if(Nvars < 10){
			std::cout << *covarmat[m] << std::endl;
		} else {
			std::cout << "Diagonal only is shown due to the large size of the matrix" << std::endl;
			std::cout << cov0.diagonal() << std::endl;
		}
		std::cout << " -------------------------------- "  << std::endl;
		mu.row(m)=vars;
	}
	std::cout << " Initial Vectors mu(m,*) initialized to the initial values of the variables " << std::endl;
	
	std::cout << "     ----------------" << std::endl;
}


void MALA::update_proposal(const VectorXd& vars, long double acceptance, int m){
/*
 * This function update the proposal law parameters
*/
	std::ofstream debug_file_stream;

	VectorXd var_p3; //, mu_m;
	MatrixXd var_p2, mat; //, covarmat_m;
	long double var_p1; //, sigma_m;

	// Update mu
	var_p3=mu.row(m) + gamma*( vars.transpose() - mu.row(m));
	mu.row(m)=p3_fct(var_p3);

	// Update the covariance matrix
	mat=(vars.transpose() - mu.row(m)).transpose() * (vars.transpose() - mu.row(m)); // This must corresponds to |x1><x2|
	var_p2=*covarmat[m] + gamma*(mat - *covarmat[m]);
	*covarmat[m]=p2_fct(var_p2);

	// Update sigma
	var_p1=sigma[m] + gamma*(acceptance - target_acceptance);
	sigma[m]=p1_fct(var_p1);
	
}

VectorXd MALA::D_MALA(){
/* The Langevin term of the Algorithm in MALA mode
 * For the moment, this is not implemented
*/ 
	VectorXd drift(Nvars);
	drift.setZero();
return drift;
}

long double MALA::multinormal_logpdf(const VectorXd& deltavars, const VectorXd& drift1, const MatrixXd& covmat1){
/* The function that ensure proper balance.
 * For the moment, this is not implemented so that only the Metropolis algorithm would work
 * For its implementation, check the function of the same name in stat.py
*/ 

return 0.;
}

VectorXd MALA::new_prop_values(const VectorXd& vars, int m){
/*
 * This function generate new samples according to the current proposal law
*/
	bool is_finite=true;
	int c;
	VectorXd ran, shift(Nvars);
	MatrixXd tmpmat;

	tmpmat=(*covarmat[m]+ epsilon2)*sigma[m];
	// Adds a small quantity in order to ensure the positivity of the matrix	
	Eigen::MatrixXd Lchol( tmpmat.llt().matrixL() ); // compute the Cholesky decomposition of the covariance matrix
	
	double *y = r8vec_normal_01 ( epsilon2.rows(), &seed ); // generate epsilon2.rows() gaussian random values
	
	// ***** Update of the proposal laws *****
	ran=vars + Lchol*Eigen::Map<VectorXd>(y, epsilon2.rows()); // Compute the multivariate distribution and project it into the parameter space
	
	c=0;
	while (is_finite == true && c<ran.size()){
		is_finite=std::isfinite(ran[c]);
		if (is_finite == false){
			std::cout << "Warning : Generation of a new vector lead to non-finite value: Recomputing..." << std::endl;
			ran=vars + Lchol*Eigen::Map<VectorXd>(y, epsilon2.rows()); // Compute the multivariate distribution and project it into the parameter space	
		}
		c=c+1;
	}
	delete y;
return ran;
}

void MALA::read_proposal(std::string readwhat, int m){
	std::cout << "                 --------------------" << std::endl;
	std::cout << "                  Read in chain m=" << m << std::endl;
	std::cout << "                 --------------------" << std::endl;
	if(readwhat == "current"){
		std::cout << "  + Current values of " << std::endl;
		std::cout << "     - mu " << std::endl;
		std::cout << mu.row(m) << std::endl;
		std::cout << "     - sigma " << std::endl;
		std::cout << sigma[m] << std::endl;
		std::cout << "     - covarmat " << std::endl;
		std::cout << *covarmat[m] << std::endl;
	}
	/*if(readwhat == "previous"){
		std::cout << "  + Previous values of " << std::endl;
		std::cout << "     - mu_prev " << std::endl;
		std::cout << mu_prev.row(m) << std::endl;
		std::cout << "     - sigma_prev " << std::endl;
		std::cout << sigma_prev[m] << std::endl;
		std::cout << "     - covarmat_prev " << std::endl;
		std::cout << *covarmat_prev[m] << std::endl;
	}
	*/
}


int MALA::parallel_tempering(Model_def *model){

	int ind_A, ind_B;	
	double* comparator=uniform_01(1, &seed );
	
	VectorXd ri_T(2), p_A, v_A, m_A;
	long double r_T, logL_A_TB, logL_B_TA, logPost_A, logPost_B, logPost_A_cross, logPost_B_cross;
	long double logPr_A, move_A, Pmove_A;

	ind_A=random_int_vals(0, Nchains-1); // integer generated from 0 to Nchains-1
	ind_B=ind_A+1; //this strategy swaps two adjacent chain: weak tunneling

	logL_A_TB=(*model).logLikelihood[ind_A]*Tcoefs[ind_A]/Tcoefs[ind_B]; // Apply the temperature Tcoefs[ind_B] on the chain of temperature Tcoef[ind_A]
	logL_B_TA=(*model).logLikelihood[ind_B]*Tcoefs[ind_B]/Tcoefs[ind_A]; // Apply the temperature Tcoefs[ind_A] on the chain of temperature Tcoef[ind_B]

	ri_T << 1., exp( logL_A_TB + logL_B_TA - (*model).logLikelihood[ind_A] - (*model).logLikelihood[ind_B]); // PRIORS DO NOT INTERVENE IF WE DEVELOP THE MATHS
	r_T=ri_T.minCoeff(); // take the minimum of the two possible values

	(*model).swaped=0;
	(*model).Pswap=0; 
		
	if (*comparator <= r_T){ // need to swap the content of the chains
		// Save A
		p_A=(*model).params.row(ind_A);
		v_A=(*model).vars.row(ind_A);
		m_A=(*model).model.row(ind_A);
		//logL_A=(*models).logLikelihood[ind_A]; // This was in my older algorithm. But here we need logL_A=logL_A_T_cross because we save the tempered likelihood
		logPr_A=(*model).logPrior[ind_A];
		move_A=(*model).moved[ind_A]; // This may have to be bypassed for the sake of the Learning algorithm (TEST ONLY AT THE MOMENT!)
		Pmove_A=(*model).Pmove[ind_A]; // This	 may have to be bypassed for the sake of the Learning algorithm (TEST ONLY AT THE MOMENT!)

		// Do A <-- B
		(*model).params.row(ind_A)=(*model).params.row(ind_B);
		(*model).vars.row(ind_A)=(*model).vars.row(ind_B);
		(*model).model.row(ind_A)=(*model).model.row(ind_B);
		(*model).logLikelihood[ind_A]=logL_B_TA; //we put the likelihood that was in B, with the temperature Tcoef[ind_A]
		(*model).logPrior[ind_A]=(*model).logPrior[ind_B];
		(*model).logPosterior[ind_A]=logL_B_TA + (*model).logPrior[ind_B];
		(*model).moved[ind_A]=(*model).moved[ind_B]; 
		(*model).Pmove[ind_A]=(*model).Pmove[ind_B]; 

		// Do B <-- A
		(*model).params.row(ind_B)=p_A;
		(*model).vars.row(ind_B)=v_A;
		(*model).model.row(ind_B)=m_A;
		(*model).logLikelihood[ind_B]=logL_A_TB;
		(*model).logPrior[ind_B]=logPr_A;
		(*model).logPosterior[ind_B]=logL_A_TB + (*model).logPrior[ind_A]; 
		(*model).moved[ind_B]=move_A; 
		(*model).Pmove[ind_B]=Pmove_A; 

		(*model).swaped=1;
		(*model).Pswap=r_T;

	} else{
		//if(ind_A == 0) {	
		//	std::cout << "No mixing" << std::endl;
		//}
	}

	model->comparator_PT=*comparator;

	delete comparator;
return ind_A; // returns the switched chain
}

void MALA::update_position_MH(Model_def *model_current, Model_def *model_propose, Data *data_struc, Config *cfg_class, int m){

	VectorXd vars_new;
	VectorXd ri(2);
	double* comparator=uniform_01(1, &seed );
	long double logproba_cur;
	long double logproba_prop;
	long double r;

	// Makes use of a grid in order to sample on the nodes of the grid
	if((*cfg_class).MALA.proposal_type == "Exact_On_Grid"){
		// Need a function that get random values on a provided grid
		// A solution can be to get random values and then to discretise
		// them in order to stick to a grid. This solution is the one
		// implemented in my Python code and might be adapted here later
		msg_handler("", "ongrid_option", "MALA::update_position_MH", (*cfg_class).MALA.proposal_type, 1);
	}
	if((*cfg_class).MALA.proposal_type == "Random"){
		vars_new=new_prop_values((*model_current).vars.row(m), m);
	}
	if((*cfg_class).MALA.proposal_type != "Random" && (*cfg_class).MALA.proposal_type != "Exact_On_Grid"){
		msg_handler("", "ongrid_option", "MALA::update_position_MH", (*cfg_class).MALA.proposal_type, 1);
	}
	(*model_propose).vars.row(m)=vars_new.transpose();
	(*model_propose).update_params_with_vars(m); // required as the model depends on params, not on vars...
	(*model_propose).generate_model(data_struc, m, Tcoefs); // This will calculate (*model_class).logLikelihood (tempered), (*model_class).logPriors and (*model_class).logPosterior
	
	if (std::isnan((*model_propose).logLikelihood[m]) == 0){
		if ((*model_propose).logPosterior[m] == -INFINITY){ // Some priors should generate infinities (e.g. Uniform priors), in that case the solution is rejected 100% of the time

 			r=0.;

		} else{
			if (use_drift == 1){
				///*
				// * If we use the drift, we need to ensure the proper balance. NOTE THAT THIS IS NOT YET IMPLEMENTED!
				//*/
				msg_handler("", "drift_option", "MALA::update_position_MH", "", 1);
				
				// ---- These lines are for handling the drift ----
				VectorXd drift1=sigma[m]*D_MALA();
				MatrixXd covmat1=sigma[m]* (*covarmat[m]);
				logproba_cur=multinormal_logpdf((*model_propose).vars.row(m) -(*model_current).vars.row(m), drift1, covmat1);
				logproba_prop=multinormal_logpdf((*model_current).vars.row(m) -(*model_propose).vars.row(m), drift1, covmat1);
				// ------------------------------------------------
			} else { 
				// nothing to do if use_drift == 0
				logproba_cur=0;
				logproba_prop=0;
			}

			
			ri << 1., exp( (*model_propose).logPosterior[m] - (*model_current).logPosterior[m] + logproba_cur - logproba_prop );
			r=ri.minCoeff(); // take the minimum of the two possible values
		}

		if (std::isnan(r) != 0){
			msg_handler("", "rejectionrate_isnan", "MALA::update_position_MH", int_to_str(m), 1);
		}
	} else { // This is if isnan((*model_propose).logLikelihood[m]) == 1
		msg_handler("", "model_isnan", "MALA::update_position_MH", int_to_str(m), 0);
		r=0.; // If we had a NaN the value is always rejected... 
			std::cout << " Information for debug " << std::endl;
			std::cout << " 		        [1] params.row(" << m << ")=" << (*model_current).params.row(m) << std::endl;
			std::cout << " 			[4] logLikelihood[" << m << "]=" << (*model_current).logLikelihood[m] << std::endl;
			std::cout << " 			[5] logPrior[" << m << "]=" << (*model_current).logPrior[m] << std::endl;
			std::cout << " 			[6] logPosterior[" << m << "]=" << (*model_current).logPosterior[m] << std::endl;
			std::cout << "     - Model propose: " << std::endl;
			std::cout << " 			[1] params.row(" << m << ")=" << (*model_propose).params.row(m) << std::endl;
			std::cout << " 			[4] logLikelihood[" << m << "]=" << (*model_propose).logLikelihood[m] << std::endl;
			std::cout << " 			[5] logPrior[" << m << "]=" << (*model_propose).logPrior[m] << std::endl;
			std::cout << " 			[6] logPosterior[" << m << "]=" << (*model_propose).logPosterior[m] << std::endl;
	}
	if (*comparator <= r){ // Case of a move in the parameter space (new model is accepted)
		(*model_current).params.row(m)=(*model_propose).params.row(m);
		(*model_current).vars.row(m)=(*model_propose).vars.row(m);
		(*model_current).model.row(m)=(*model_propose).model.row(m);
		(*model_current).logLikelihood[m]=(*model_propose).logLikelihood[m];
		(*model_current).logPrior[m]=(*model_propose).logPrior[m];
		(*model_current).logPosterior[m]=(*model_propose).logPosterior[m];
		(*model_current).moved[m]=1;

	} else{
		(*model_current).moved[m]=0;

	}
	(*model_current).Pmove[m]=r; // In all case we save the move probability
	(*model_current).comparator_MH[m]=*comparator; // In all case we save the move probability

	delete comparator;
}

void MALA::execute(Model_def *model_current, Model_def *model_propose, Data *data_struc, Config *cfg_class, Outputs *outputs_class, Diagnostics *diags){
/*
 * The class MALA must have been fully initialized before proceeding.
 * The model_current must also have been fully initialized.
 * Content of the passed classes:
 *     - model_current: Fully initialized model. This includes priors
 *     - model_propose: Fully initialized model. The configuration must be the same as model_current
 *     - data_struc: A structure with the observational data
 *     - cfg_class: a class which contains the full configuration. Typically, this
 *                  would have been initialized by reading a configuration file.
 *     - outputs_class: a class that contains all function handling output data files (not the graphic diagnostics)
 *     - diags: a class that contains all function handling the diagnostics (graphics and non summary of the results)
*/

	long time0, time1, timedelta;
	bool logic, logic_old; // booleans that tell when we learn the var-covariance and when we don't. logic_old is used to output the Learning status...
	bool tempted_mixing; // 1 if we fulfill the condition for parallel-chains mixing
	int which, which_old;
	long Nswap=0, swap_period=1000;
	double Nestimate=0, swaping_rate=0;
	//VectorXd vec_tmp; // required to let the compiler happy when passing (*model_current).vars.row(chain) into MALA::update_proposal()
	VectorXd mean_params_Nbuf(model_current->params.cols());
	VectorXd mean_params_Nsamp(model_current->params.cols());
	VectorXd mean_model; // This contains the averaged model either over Nbuffer or Nsamples
	int chain_A;
        const int MaxChain = 24;
	int chain;
	bool logicA[MaxChain];
	int  whichA[MaxChain];

        if (Nchains > MaxChain) {
		msg_handler("", "Nchains_error", "MALA::execute", "Number of chains="+ int_to_str(Nchains) + " exceeds MaxChain=" + int_to_str(MaxChain), 1);
	}

	if (Nsamples < outputs_class->get_Nbuffer()){
		 swap_period=Nsamples;
	} else{
		swap_period=outputs_class->get_Nbuffer();
	}

	mean_params_Nbuf.setZero(); // This is used to generate plot over Nbuffer (intermediate plots) 
	mean_params_Nsamp.setZero(); // This is used to generate plot over Nsamples (final plot)

	// Show the initial plot
	diags->gnuplt_model_diags(data_struc, model_current->model.row(0), "init"); // We show only the fit for the model with the lowest temperature
	
	long i=initial_i;
	long count=0;
	logic_old=0;
	which_old=-1;
	if (i>=Nsamples-1){
		msg_handler("", "iteration_error", "MALA::execute", "initial i="+ int_to_str(i) + " >= Nsamples-1=" + int_to_str(Nsamples-1), 1);
	}
	std::cout << " ------------------------------------------------" << std::endl;
	std::cout << "             BEGINING THE MCMC PROCESS           " << std::endl;
	std::cout << " ------------------------------------------------" << std::endl;
	std::cout << "                 Learning is OFF                 " << std::endl;
	if (dN_mixing >= Nsamples){
		std::cout << "        No Mixing of the parallel chains        " << std::endl;
	} else{
		std::cout << "      Mixing rate of the parallel chains : " << dN_mixing << std::endl;
	}
	
	time0=time1=ben_clock(); // Time at the begining of the Sampling
	Nestimate=1000.;
	if(Nsamples < Nestimate){
		Nestimate=Nsamples/2.;
	};
	while (i<Nsamples){
		if (i%cfg_class->outputs.Nbuffer == 0 || (i==Nsamples-1)){ 
			timedelta=ben_clock() - time1; // Time interval since last writting operation
			time1=ben_clock();
			std::cout << "[" << i << "]  " << std::endl; 
			if( i!=0 ){
				std::cout << "   - Total time so far: " << float(time1-time0)/ 60. << " min = ";
				std::cout << float(time1-time0)/ 3600. << " hours = ";
				std::cout << float(time1-time0)/ 86400. << " days" << std::endl;
				if(i != Nsamples-1){
					std::cout << "   - Processing rate: " << cfg_class->outputs.Nbuffer/(float(timedelta)/ 60.) << " samples/min (";
					std::cout << cfg_class->outputs.Nbuffer/(float(timedelta)/ 3600.) << " samples/hours)" << std::endl;				
				}
			} else{
				std::cout << "An estimate of the time required before writting next buffer will be given once " << Nestimate << " samples have been processed" << std::endl;
			}
		}
		if (i == Nestimate){
			std::cout << "   - Estimated processing rate: " << Nestimate /(float(ben_clock() - time0)/ 3600.);
			std::cout << " samples/hour" << std::endl;
		}

		gamma=c0/(1. + i);

#pragma omp parallel for num_threads(Nchains) default(shared) private(chain)
		for (chain=0; chain<Nchains; chain++){
			VectorXd vec_tmp;

			// [1] Propose a new position and test it so that the model_current is updated (if new position accepted) using model_propose
			//     Note that model_propose must have been initialized in order to work
			update_position_MH(model_current, model_propose, data_struc, cfg_class, chain);
			// [2] Test whether we update the variance-covariance parameters
			logicA[chain]=0; // reset
			whichA[chain]=0; 
			for(int l=0; l<periods_learn.size(); l++){ // note that periods.size() must be Nt_learn.size() - 1
				logicA[chain]= logicA[chain] || ((i>=Nt_learn[l]) && (i<Nt_learn[l+1])); //Test whether we have to learn
				if ((i>=Nt_learn[l]) && (i<Nt_learn[l+1])){ // Tells which phase we are actually doing
					whichA[chain]=l;
				}
			}
			if (logicA[chain] == 1 && ( i%periods_learn[whichA[chain]]) == 0){
				vec_tmp=(*model_current).vars.row(chain);
				update_proposal(vec_tmp, (*model_current).Pmove[chain], chain);	
			} //else { std::cout << " NO Learning for chain " << chain << std::endl; }
		}
		for (int chain=0; chain<Nchains; chain++){
			if(logicA[chain] > logic_old){ // case of logic=1 and logic_old=0
				std::cout << "        ["<< i << "] Learning turned ON...Learning periodicity: " << periods_learn[whichA[chain]] << std::endl;
				logic_old=logicA[chain];
				which_old=whichA[chain];
			}
			if(logicA[chain] < logic_old){ // case of logic=0 and logic_old=1
				std::cout << "        ["<< i << "] Learning turned OFF" << std::endl;
				logic_old=logicA[chain];
				which_old=whichA[chain];
			}
			if(logicA[chain] == 1 && whichA[chain] != which_old){ // case of logic=0 and logic_old=1
				std::cout << "        ["<< i << "] Learning maintained ON... Learning periodicity changed to: "<< periods_learn[which] << std::endl;
				logic_old=logicA[chain];
				which_old=whichA[chain];
			}
		}
		
		// [3] Parallel tempering
		if (i%dN_mixing == 0 && i!=0){
			chain_A=parallel_tempering(model_current);
			if (model_current->swaped == 1){
				swaping_rate=swaping_rate + 1.0; // Used to calculate the average swaping rate
			}
			Nswap=Nswap+1;
			tempted_mixing=1;
			if( (Nswap == swap_period) && (i!=0)){
				std::cout << "        ["<< i << "] Average swaping rate for the last " << swap_period << " swaps (in %): " << 100. * swaping_rate/Nswap << std::endl;
				Nswap=0; // re-initialize the average counter
				swaping_rate=0; // re-initialize the average swaping rate
			}			
		} else {
			tempted_mixing=0;
			chain_A=-1;
		}
		// [4] Show / save results
		// Update the restoration point when counts == Nbuffer (conditon handled within the method)
		outputs_class->update_buffer_restore(sigma, mu, covarmat, model_current->vars);
		// Generate the outputs... If requested
		outputs_class->update_buffer_stat_criteria(model_current->logLikelihood, model_current->logPrior, model_current->logPosterior); 	
		//std::cout << "Before writting params" << std::endl;
		outputs_class->update_buffer_params(model_current->vars);
		//std::cout << "Before writting ptempering" << std::endl;
		outputs_class->update_buffer_ptempering(tempted_mixing, chain_A, model_current->Pswap, model_current->swaped);
		//std::cout << "Before writting proposals" << std::endl;
		outputs_class->update_buffer_proposals(sigma, mu, model_current->Pmove, 
			    			model_current->moved, covarmat);
		outputs_class->update_buffer_models(model_current->model); //IF ACTIVATED IT SLOWS DOWN THE PROCESS AND CREATE HUGE FILES!!! OVERRID THIS FUNCTION WHEN Ndata_Obs > Ndata_Obs_limit
		// Generate diagnostics... If requested... REQUIRES THAT THE APPROPRIATE OUTPUTS EXIST
		if(diags->get_model_buffer_diags() == 1){
			mean_params_Nbuf=(mean_params_Nbuf.transpose() + model_current->params.row(0)); // average made only over the coldest chain
			if(count == (outputs_class->get_Nbuffer())){
				model_propose->params.row(0)=mean_params_Nbuf/outputs_class->get_Nbuffer();
				mean_model=model_propose->call_model(data_struc, 0); // model for the coldest chain m=0
				diags->gnuplt_model_diags(data_struc, mean_model, "buffer"); // plot on the eps file for the buffer

				mean_params_Nbuf.setZero(); // initialize the buffer for the mean parameters
				count=0;			
			}
			count=count+1;
		}
		if (diags->get_model_final_diags() == 1){
			mean_params_Nsamp=(mean_params_Nsamp.transpose() + model_current->params.row(0)); // average made only over the coldest chain
		}
		if((i == Nsamples-1) && (diags->get_model_final_diags() == 1)){
			model_propose->params.row(0)=mean_params_Nsamp/Nsamples;
			mean_model=model_propose->call_model(data_struc, 0); // model for the coldest chain m=0
			diags->gnuplt_model_diags(data_struc, mean_model, "final"); // plot on the eps file for the final plot

			mean_params_Nsamp.setZero(); // initialize the buffer for the mean parameters
		}
		diags->gnuplt_chains_diags(i, Tcoefs); // diagnostics for the acceptance rate and the likelihoods
	  	diags->gnuplt_pdfs_diags_main(i); // diagnostics for the pdfs. Note that if requested, this slows down the computation

		i=i+1;

	}
	
}

/*
 * Function that show an error message depending on error_type
 * It receives the filename, fonction name specific message values (arguments) and returns an error code
 * A negative error code correspond to a fatal error. A positive
 * code correspond to a warning.
 * Exit from the program is controlled by the boolean fatal
*/
int MALA::msg_handler(const std::string file, const std::string error_type, const std::string fct_name, const std::string arguments, const bool fatal){

	bool err_msg=0;

	if (fatal == 1){
		std::cout << "Warning Fatal error :" << std::endl;
	}
	if (file != ""){ std::cout << "          Accessing file: " << file << std::endl;}
	if (fct_name !=""){ std::cout << "          Fonction name: " << fct_name << std::endl;}

	if(error_type == "rejectionrate_isnan"){
			std::cout << " ******* FATAL ERROR ******" << std::endl;
			std::cout << " got r=min([1, exp(logP_propose - logP_current)] = NaN  for chain m=" << arguments << std::endl;
			std::cout << " Need a serious debuging here!!! " << std::endl;
			std::cout << " ------------------------" << std::endl;
	}
	if(error_type == "model_isnan"){
		std::cout << " ******* Warning ********" << std::endl;
		std::cout << " logL = NaN detected... for chain m=" << arguments << "We highly recommend to avoid such case to happen as its slows down the computation " << std::endl;
		std::cout << " To avoid this you need to modify your model so that NaN are outside the values allowed by your priors (e.g. take abs(param[i])) " << std::endl;
		std::cout << " Solution rejected " << std::endl;
		std::cout << " ------------------------" << std::endl;

	}
	if(error_type == "ongrid_option"){
		std::cout << "proposal_type = " << arguments << std::endl;
		std::cout << "This functionality is not yet available!" << std::endl;
		std::cout << "Use proposal_type = Random instead" << std::endl;
	}
	if(error_type == "drift_option"){
				std::cout << " ------------------------" << std::endl;
				std::cout << " use_drift is set to 1 : " << std::endl;
				std::cout << " The MALA algorithm has not been implementing the drift yet!" << std::endl;
				std::cout << " This will not work! " << std::endl;
				std::cout << " Please set use_drift to 0 so that the programs works in MH mode" << std::endl;
				std::cout << " In MH mode, the following (unfinished) functions are not used: " << std::endl;
				std::cout << "         - MALA::D_MALA() " << std::endl;
				std::cout << "         - MALA::multinomal_logpdf() " << std::endl;
				std::cout << " ------------------------" << std::endl;
	}
	if(error_type == "iteration_error"){
		std::cout << " This is due to an improper use of restored proposal law" << std::endl;
		std::cout << "   - Set Nsamples to a greater value than the initial i in order to properly complete the previously unfinished run" << std::endl;
		std::cout << "   - Beware of the setups for do_restore_proposal and do_restore_variables" << std::endl;
		err_msg=1;
	}
	if(error_type == "Nchains_error"){
		std::cout << "You need to set a lower value for Nchains or change the limit in the code" << std::endl;
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

    	if(error_type == "" || err_msg == 0){
    	    std::cout << "Unknown Fatal Error in " << fct_name << std::endl;
    	    std::cout << "Could not open the file: " << file << std::endl;
    	    std::cout << "Neeed debuging..." << std::endl;
    	}
	if(fatal == 1){
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
		return -1;
	} else{
		return 1;
	}
	
}

