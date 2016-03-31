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
#include "MALA.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using namespace std;

MatrixXd** initialize_3dMatrix(int Nchains, int Nrows, int Ncols);
void set_3dMatrix(MatrixXd** matrix3d, MatrixXd m_in, int position);
MatrixXd** copy_3dMatrix(MatrixXd** m3d_in, int depth);
MatrixXd multivariate_eigen(VectorXd, MatrixXd, int);
double *r8vec_normal_01 ( int n, int *seed );
double *uniform_01 ( int n, int *seed );

MALA::MALA(long Nspl, long Nchs, VectorXd params, vector<string> params_names, VectorXi relax){
/* 
 * This procedure loads the initial configuration.
 * For the moment all is hard coded. But ultimately,
 * the configuration might be read from a configuration file
*/
	Nsamples=Nspl;
	Nchains=Nchs;
	seed=(unsigned)time(NULL);
	srand(seed);

	// Identify the variables which in turns will enable us to initialize the proposal parameters (e.g. covarmat)
	Nvars=0; // This is a variable of the class... see MALA.h
	VectorXd vars(params.size());
	vector<string> vars_names(params.size());
	for (int i=0; i<params.size();i++){
		if(relax[i] == 1){
			vars(Nvars)=params(i);
			vars_names[Nvars]=params_names[i];
			Nvars=Nvars+1;
		}
	}
	vars.conservativeResize(Nvars);
	vars_names.resize(Nvars);
	init_proposal(vars, vars_names); // Initialize the covariance matrix and all parameters relative to the proposal

	VectorXd tmp(Nvars);
	tmp.setConstant(1e-8);

	epsilon1=1e-8; // Bundary for the Robins-Monroe stochastic algorithm. Must be as small as the numerical precision of the machine. With long double, can THEORITICALLY be 1e-19 (1e-16 for double)
	epsilon2=tmp.asDiagonal(); // Bundary for the Robins-Monroe stochastic algorithm. Must be as small as the numerical precision of the machine.
	A1=1e10; // Bundary for the Robins-Monroe stochastic algorithm. Must be as large as the numerical precision of the machine allows.
	delta=0; // This parameter is only for the MALA algorithm (ie, when the gradient is available). For a RWM, set to 0...
	delta_x=1e-10;
	c0=5;
	lambda_temp=1.8; //value should always be adjusted so that Tmax>120
	use_drift=0; // When 0, use the MH scheme instead of the MALA scheme

	// Just for test purpose, define three zone where Learning is made in different ways
	Nt_learn.resize(3); // Controls when we do the learning
	periods_learn.resize(2);  // For each values of Nt_learn, tells at which rate the learning is done
	Nt_learn[0]=10;  // DO NOT LEARN BELLOW THAT VALUE
	Nt_learn[1]=100; 
	Nt_learn[2]=Nsamples; // DO NOT LEARN AFTER THAT VALUE
	periods_learn[0]=5; // if Nt_learn[0]< i <Nt_learn[1], learn every 10 times
	periods_learn[1]=1; // for all Nt_learn[1]< i <Nt_learn[2], learn every times
        // Note that if you set Nt_learn[2] <  Nsamples, then no learning is done after Nt_learn[2]

	dN_mixing=5.; // Perform a chain mixing every dN_mixing

	Tcoefs.resize(Nchains);	
	for (int m=0; m<Nchains; m++){Tcoefs[m]=std::pow(lambda_temp, m-1);} // initialize the list of temperatures

}

double MALA::p1_fct(double x, double epsilon, double A){
	
	double scalar_p;
	if(x >= epsilon && x<=A){
		scalar_p=x;
	}
	if(x < epsilon){
		scalar_p=A;
	}
	if(x > A){
		scalar_p=A;
	}
return scalar_p;
}

MatrixXd MALA::p2_fct(MatrixXd x, double A){

	MatrixXd matrix_p;
	double norm_x=x.norm(); // Frobenius norm

	if(norm_x <= A){
		matrix_p=x;
	}
	if(norm_x > A){
		matrix_p=x * A/norm_x;
	}
	return matrix_p;
}

VectorXd MALA::p3_fct(VectorXd x, double A){
	VectorXd vector_p;
	double norm_x=x.norm(); // Frobenius norm
	
	if(norm_x <=A){
		vector_p=x;
	}
	if(norm_x > A){
		vector_p=x * A/norm_x;
	}
	return vector_p;
}

/*
double MALA::set_temperature(const long double likelihood, const long double lambda, const long j){

	long double T=std::pow(lambda, j-1);  // j is the chain number

return likelihood/T;
}
*/

long MALA::random_int_vals(int N0, int Nmax){
/* 
 * Generate a uniform random number between 1 and Nmax
 * WARNING: srand(time(NULL)) must be put in the main
 *          program in order to generate proper
 *          pseudo-random number!
 */

return N0 + rand()%Nmax; // srand(seed) is exectuded by the constructor of the MALA class
}

void MALA::init_proposal(VectorXd vars, vector<string> var_names){

	VectorXd error(Nvars);
	MatrixXd cov0;

	error.setConstant(1); // Defaut non-zero value;
	for(long i=0; i<vars.size(); i++){
		if(var_names[i] == "Height_l0" || var_names[i] == "Height_l1" || var_names[i] == "Height_l2" || var_names[i] == "Height_l3"){
			error[i]=std::abs(vars[i])*0.2 + 0.1;
		}
		if(var_names[i] == "Width_l0" || var_names[i] == "Width_l1" || var_names[i] == "Width_l2" || var_names[i] == "Width_l3"){
			error[i]=std::abs(vars[i])*0.15 + 0.05;
		}
		if(var_names[i] == "Visibility_l0" || var_names[i] == "Visibility_l1" || var_names[i] == "Visibility_l2" || var_names[i] == "Visibility_l3"){
			error[i]=0.2;
		}
		if(var_names[i] == "Frequency_l0" || var_names[i] == "Frequency_l1"){
			error[i]=1.5;
		}	
		if(var_names[i] == "Frequency_l2"){
			error[i]=2.0;
		}
		if(var_names[i] == "Frequency_l3"){
			error[i]=2.5;
		}
		if(var_names[i] == "Splitting_a1"){
			error[i]=std::abs(vars[i])*0.2 + 0.2;
		}
		if(var_names[i] == "Splitting_eta"){
			error[i]=std::abs(vars[i])*2. + 1e-5;
		}
		if(var_names[i] == "Splitting_a3"){
			error[i]=std::abs(vars[i])*1.5 + 0.02;
		}
		if(var_names[i] == "Splitting_mag_b"){
			error[i]=std::abs(vars[i])*0.5 + 0.5;
		}
		if(var_names[i] == "Splitting_mag_alfa"){
			error[i]=1.0; //std::abs(vars[i]);
		}
		if(var_names[i] == "Lorentzian_asymetry"){
			error[i]=std::abs(vars[i])*0.1 + 5.;
		}
		if(var_names[i] == "Inclination"){
			error[i]=20;
		}
		if(var_names[i] == "Harvey-Noise_H"){
			error[i]=std::abs(vars[i])*0.3;
		}
		if(var_names[i] == "Harvey_Noise_tc"){
			error[i]=std::abs(vars[i])*0.15;
		}
		if(var_names[i] == "Harvey_Noise_p"){
			error[i]=std::abs(vars[i])*0.15;
		}
		if(var_names[i] == "White_Noise_N0"){
			error[i]=std::abs(vars[i])*0.05;
		}
	}
	
	std::cout << "Nsamples = " << Nsamples << std::endl;
	std::cout << "Nvars = " << Nvars << std::endl;
	std::cout << "Nchains = " << Nchains << std::endl;
	std::cout << " ------- " << Nchains << std::endl;
        error=error.array().square();
	cov0=error.asDiagonal();
	covarmat=initialize_3dMatrix(Nchains, Nvars, Nvars);
	covarmat_prev=initialize_3dMatrix(Nchains, Nvars, Nvars);
	sigma.setZero(Nchains);
	sigma_prev.setZero(Nchains);
	mu.setZero(Nchains, Nvars);
	mu_prev.setZero(Nchains, Nvars);
	for(int m=0; m<Nchains; m++){
		set_3dMatrix(covarmat, cov0, m); // the error variance is the initial matrix
		set_3dMatrix(covarmat_prev, cov0, m); // the error variance is the initial matrix
		sigma[m]=std::pow(2.38,2)*(0.1*m+1.)/Nvars;
		sigma_prev[m]=std::pow(2.38,2)*(0.1*m+1.)/Nvars;
		std::cout << "scaling factor sigma[" << m << "]=" << sigma[m] << std::endl;
		std::cout << "Initial Covariance Matrix covarmat["<< m << "]=" << std::endl;
		std::cout << *covarmat[m] << std::endl;
		std::cout << " -------------------------------- "  << std::endl;
		mu.row(m)=vars;
	}
	mu_prev=mu;
	std::cout << " Initial Vectors mu(m,*) initialized to the initial values of the variables " << std::endl;
	std::cout << " Vector of initial variables: " << std::endl;
	std::cout << vars.transpose() << std::endl;

}

void MALA::update_proposal(VectorXd vars, double acceptance, int m){
	
/*
 * This function update the proposal law parameters
*/
	VectorXd var_p3, mu_m;
	MatrixXd var_p2, mat, covarmat_m;
	double var_p1, sigma_m;
	
	// Store old values for the mu vector, sigma and the covariance matrix
	sigma_prev[m]=sigma[m];
	mu_prev.row(m)=mu.row(m);
	set_3dMatrix(covarmat_prev, *covarmat[m], m);

	// Update mu
	var_p3=mu.row(m) + gamma*( vars.transpose() - mu.row(m));
	mu.row(m)=p3_fct(var_p3, A1);

	// Update the covariance matrix
	mat=(vars.transpose() - mu.row(m)).transpose() * (vars.transpose() - mu.row(m)); // This must corresponds to |x1><x2|
	var_p2=*covarmat[m] + gamma*(mat - *covarmat[m]);
	*covarmat[m]=p2_fct(var_p2, A1);

	// Update sigma
	var_p1=sigma[m] + gamma*(acceptance - target_acceptance);
	sigma[m]=p1_fct(var_p1, epsilon1, A1);

	//read_proposal("previous", m);
	//read_proposal("current", m);

}

VectorXd MALA::D_MALA(){
/* The Langevin term of the Algorithm in MALA mode
 * For the moment, this is not implemented
*/ 
	VectorXd drift(Nvars);
	drift.setZero();
return drift;
}

long double MALA::multinormal_logpdf(VectorXd deltavars, VectorXd drift1, MatrixXd covmat1){
/* The function that ensure proper balance.
 * For the moment, this is not implemented so that only the Metropolis algorithm would work
 * For its implementation, check the function of the same name in stat.py
*/ 

return 0.;
}

VectorXd MALA::new_prop_values(VectorXd vars, int m){
/*
 * This function generate new samples according to the current proposal law
*/
	VectorXd ran, shift(Nvars);

	//for(int ii=0; ii<Nchains; ii++){
	//	std::cout << "covarmat[" << ii << "] =" << *covarmat[m] << std::endl;
	//}
	covarmat_prev=copy_3dMatrix(covarmat, Nchains); // do a hard copy
	//covarmat_prev=covarmat; // copy by reference cannot be done here
	sigma_prev=sigma; // This works as it seems to do a hard copy

	//std::cout << "seed value =" << seed << std::endl;
	Eigen::MatrixXd Lchol( (*covarmat[m]+ epsilon2).llt().matrixL() );	
	
	double *y = r8vec_normal_01 ( (*covarmat[m]+ epsilon2).rows(), &seed );
	
	// ***** Update of the proposal laws *****
	// Adds a small quantity in order to ensure the positivity of the matrix
	// And generate a vector of random numbers
	ran=mu.row(m).transpose() + Lchol*Eigen::Map<VectorXd>(y, (*covarmat[m]+ epsilon2).rows());

return ran;
}

void MALA::read_proposal(string readwhat, int m){
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
	if(readwhat == "previous"){
		std::cout << "  + Previous values of " << std::endl;
		std::cout << "     - mu_prev " << std::endl;
		std::cout << mu_prev.row(m) << std::endl;
		std::cout << "     - sigma_prev " << std::endl;
		std::cout << sigma_prev[m] << std::endl;
		std::cout << "     - covarmat_prev " << std::endl;
		std::cout << *covarmat_prev[m] << std::endl;
	}

}

int MALA::parallel_tempering(Model_def *model){

	int ind_A, ind_B;	
	double* comparator=uniform_01(1, &seed );
	VectorXd ri_T(2), p_A, v_A, m_A;
	long double r_T, logL_A_T_cross, logL_B_T_cross, logPost_A, logPost_B, logPost_A_cross, logPost_B_cross;
	long double logPr_A, move_A, Pmove_A;


	ind_A=random_int_vals(0, Nchains-1); // integer generated from 0 to Nchains-1
	ind_B=ind_A+1; //this strategy swaps two adjacent chain: weak tunneling

	//logL_A_T=(*model).logLikelihood[ind_A]; // We do not need this because we have (*model).logPosterior[ind_A]
	//logL_B_T=(*model).logLikelihood[ind_B]; // We do not need this because we have (*model).logPosterior[ind_B]
	logL_A_T_cross=(*model).logLikelihood[ind_A]*Tcoefs[ind_A]/Tcoefs[ind_B]; // Apply the temperature Tcoefs[ind_B] on the chain of temperature Tcoef[ind_A]
	logL_B_T_cross=(*model).logLikelihood[ind_B]*Tcoefs[ind_B]/Tcoefs[ind_A]; // Apply the temperature Tcoefs[ind_A] on the chain of temperature Tcoef[ind_B]
	
	logPost_A=(*model).logPosterior[ind_A]; // Retrieve the logProba for chain ind_A
	logPost_B=(*model).logPosterior[ind_B]; // Retrieve the logProba for chain ind_B
	logPost_A_cross=logL_A_T_cross + (*model).logPrior[ind_A]; // Posterior of chain ind_A, heated at Tcoefs[ind_B]
	logPost_B_cross=logL_B_T_cross + (*model).logPrior[ind_B]; // Posterior of chain ind_B, heated at Tcoefs[ind_A]

	ri_T << 1., exp( logPost_A_cross + logPost_B_cross - logPost_A - logPost_B);
	r_T=ri_T.minCoeff(); // take the minimum of the two possible values

	(*model).swaped=0;
	//(*model).swaped.assign(Nchains, 0); // Set all values to 0. Required because swaped chains are randomly chosen
	(*model).Pswap=0; // Set all values to 0.

	if (*comparator <= r_T){ // need to swap the content of the chains
		// Save A
		p_A=(*model).params.row(ind_A);
		v_A=(*model).vars.row(ind_A);
		m_A=(*model).model.row(ind_A);
		//logL_A=(*models).logLikelihood[ind_A]; // This was in my older algorithm. But here we need logL_A=logL_A_T_cross because we save the tempered likelihood
		logPr_A=(*model).logPrior[ind_A];
		move_A=(*model).moved[ind_A]; // This may have to be bypassed for the sake of the Learning algorithm (TEST ONLY AT THE MOMENT!)
		Pmove_A=(*model).Pmove[ind_A]; // This may have to be bypassed for the sake of the Learning algorithm (TEST ONLY AT THE MOMENT!)

		// Do A <-- B
		(*model).params.row(ind_A)=(*model).params.row(ind_B);
		(*model).vars.row(ind_A)=(*model).vars.row(ind_B);
		(*model).model.row(ind_A)=(*model).model.row(ind_B);
		(*model).logLikelihood[ind_A]=logL_B_T_cross; //we put the likelihood that was in B, with the temperature Tcoef[ind_A]
		(*model).logPrior[ind_A]=(*model).logPrior[ind_B];
		(*model).logPosterior[ind_A]=logPost_B_cross;
		(*model).moved[ind_A]=1; // EXPERIMENTAL... We put override the result from the MH by the Parallel tempering! IF IT DOES NOT WORK, PUT: (*model).moved[ind_B]
		(*model).Pmove[ind_A]=r_T; // EXPERIMENTAL... We put override the result from the MH by the Parallel tempering! IF IT DOES NOT WORK, PUT: (*model).Pmove[ind_B]
	
		// Do B <-- A
		(*model).params.row(ind_B)=p_A;
		(*model).vars.row(ind_B)=v_A;
		(*model).model.row(ind_B)=m_A;
		(*model).logLikelihood[ind_B]=logL_A_T_cross; 
		(*model).logPrior[ind_B]=logPr_A;
		(*model).logPosterior[ind_B]=logPost_A_cross; 
		(*model).moved[ind_B]=1; // EXPERIMENTAL... We put override the result from the MH by the Parallel tempering! IF IT DOES NOT WORK, PUT: move_A
		(*model).Pmove[ind_B]=1.-r_T; // EXPERIMENTAL... We put override the result from the MH by the Parallel tempering! IF IT DOES NOT WORK, PUT: Pmove_a

		(*model).swaped=1;
		//(*model).swaped[ind_A]=1;
		//(*model).swaped[ind_B]=1;
		(*model).Pswap=r_T;
		//std::cout << "Mixing done with a probability : " << r_T << std::endl;
	}

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
		std::cout << "proposal_type = Exact_On_Grid" << std::endl;
		std::cout << "This functionality is not yet available!" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit ( EXIT_FAILURE );
	}
	if((*cfg_class).MALA.proposal_type == "Random"){
		//std::cout << " Before new prop values " << std::endl;
		//std::cout << " (*model_current).vars.row(m) = " << (*model_current).vars.row(m) << std::endl;
		vars_new=new_prop_values((*model_current).vars.row(m), m);
		//std::cout << " vars_new = " << vars_new << std::endl;
		
	}
	if((*cfg_class).MALA.proposal_type != "Random" && (*cfg_class).MALA.proposal_type != "Exact_On_Grid"){
		std::cout << "proposal_type =" << (*cfg_class).MALA.proposal_type << std::endl;
		std::cout << "This proposal_type is not recognized!" << std::endl;
		std::cout << "Allowed keywords: Random or Exact_On_Grid" << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit ( EXIT_FAILURE );
	}

	(*model_propose).vars.row(m)=vars_new.transpose();
	// ---- Some debug lines ----
	//std::cout << "     params(" << m << ") before update..." << std::endl;
	//std::cout << (*model_propose).params.row(m) << std::endl;
	(*model_propose).update_params_with_vars(m); // required as the model depends on params, not on vars...
	// ---- Some debug lines ----	
	//std::cout << "     params(" << m << ") after update..." << std::endl;
	//std::cout << (*model_propose).params.row(m) << std::endl;
	long double logP=(*model_propose).generate_model(data_struc, m, Tcoefs); // This will calculate (*model_class).logLikelihood (tempered), (*model_class).logPriors and (*model_class).logPosterior

	//std::cout << "logP=" << logP << std::endl;
	
	if (isnan((*model_propose).logLikelihood[m]) == 0){
		if ((*model_propose).logPosterior[m] == -INFINITY){ // Some priors should generate infinities (e.g. Uniform priors), in that case the solution is rejected 100% of the time

 			r=0.;

		} else{
			if (use_drift == 1){
				///*
				// * If we use the drift, we need to ensure the proper balance. NOTE THAT THIS IS NOT YET IMPLEMENTED!
				//*/
				std::cout << " ******* Warning ******" << std::endl;
				std::cout << " use_drift is set to 1 : " << std::endl;
				std::cout << " The MALA algorithm has not been finished/tested yet!" << std::endl;
				std::cout << " This will not work! " << endl;
				std::cout << " Please set use_drift to 0 so that the programs works in MH mode" << std::endl;
				std::cout << " In MH mode, the following (unterminated) functions are not used: " << std::endl;
				std::cout << "         - MALA::D_MALA() " << std::endl;
				std::cout << "         - MALA::multinomal_logpdf() " << std::endl;
				std::cout << " The program will exit now " << std::endl;
				std::cout << " ------------------------" << std::endl;
				exit ( EXIT_FAILURE );
				
				VectorXd drift1=sigma[m]*D_MALA();
				MatrixXd covmat1=sigma[m]* (*covarmat[m]);
				logproba_cur=multinormal_logpdf((*model_propose).vars.row(m) -(*model_current).vars.row(m), drift1, covmat1);
				logproba_prop=multinormal_logpdf((*model_current).vars.row(m) -(*model_propose).vars.row(m), drift1, covmat1);

			} else { 
				// nothing to do if use_drift == 0
				logproba_cur=0;
				logproba_prop=0;
			}

			//std::cout << "(*model_propose).logPosterior[" << m << "]=" << (*model_propose).logPosterior[m] << std::endl;
			//std::cout << "(*model_current).logPosterior[" << m << "]=" << (*model_current).logPosterior[m] << std::endl;
			
			ri << 1., exp( (*model_propose).logPosterior[m] - (*model_current).logPosterior[m] + logproba_cur - logproba_prop );
			//std::cout << "ri=" << ri << std::endl;
			r=ri.minCoeff(); // take the minimum of the two possible values
		}
		if (isnan(r) != 0){
			std::cout << " ******* FATAL ERROR ******" << std::endl;
			std::cout << " got r=min([1, exp(logP_propose - logP_current)] = NaN  " << std::endl;
			std::cout << " Need a serious debuging here!!! " << std::endl;
			std::cout << " The program will exit now " << std::endl;
			std::cout << " ------------------------" << std::endl;
			std::cout << " Information for debug " << std::endl;
			std::cout << "     - Model current: " << std::endl;
			std::cout << " 			[1] params.row(" << m << ")=" << (*model_current).params.row(m) << std::endl;
			std::cout << " 			[2] vars.row(" << m << ")=" << (*model_current).vars.row(m) << std::endl;
			//std::cout << " 			[3] model.row(" << m << ")=" << (*model_current).model.row(m) << std::endl;
			std::cout << " 			[4] logLikelihood[" << m << "]=" << (*model_current).logLikelihood[m] << std::endl;
			std::cout << " 			[5] logPrior[" << m << "]=" << (*model_current).logPrior[m] << std::endl;
			std::cout << " 			[6] logPosterior[" << m << "]=" << (*model_current).logPosterior[m] << std::endl;
			std::cout << " 			[7] moved[" << m << "]=" << (*model_current).moved[m] << std::endl;
			std::cout << "     - Model propose: " << std::endl;
			std::cout << " 			[1] params.row(" << m << ")=" << (*model_propose).params.row(m) << std::endl;
			std::cout << " 			[2] vars.row(" << m << ")=" << (*model_propose).vars.row(m) << std::endl;
			//std::cout << " 			[3] model.row(" << m << ")=" << (*model_propose).model.row(m) << std::endl;
			std::cout << " 			[4] logLikelihood[" << m << "]=" << (*model_propose).logLikelihood[m] << std::endl;
			std::cout << " 			[5] logPrior[" << m << "]=" << (*model_propose).logPrior[m] << std::endl;
			std::cout << " 			[6] logPosterior[" << m << "]=" << (*model_propose).logPosterior[m] << std::endl;
			std::cout << " 			[7] moved[" << m << "]=" << (*model_propose).moved[m] << std::endl;
			exit ( EXIT_FAILURE );
		}
	} else { // This is if isnan((*model_propose).logLikelihood[m]) == 1
		r=0.; // If we had a NaN the value is always rejected... 
		std::cout << " ******* Warning ********" << std::endl;
		std::cout << " logL = NaN detected... We highly recommend to avoid such case to happen as its slows down the computation " << std::endl;
		std::cout << " To avoid this you need to modify your model so that NaN are outside the values allowed by your priors (e.g. take abs(param[i])) " << std::endl;
		std::cout << " Solution rejected " << std::endl;
		std::cout << " ------------------------" << std::endl;
	}
	if (*comparator <= r){ // Case of a move in the parameter space (new model is accepted)
		//std::cout << "                     chain [" << m << "] proposed model ACCEPTED! comparator=" << *comparator << "  r=" << r << std::endl;
		(*model_current).params.row(m)=(*model_propose).params.row(m);
		(*model_current).vars.row(m)=(*model_propose).vars.row(m);
		(*model_current).model.row(m)=(*model_propose).model.row(m);
		(*model_current).logLikelihood[m]=(*model_propose).logLikelihood[m];
		(*model_current).logPrior[m]=(*model_propose).logPrior[m];
		(*model_current).logPosterior[m]=(*model_propose).logPosterior[m];
		(*model_current).moved[m]=1;
	} else{
		//std::cout << "                     chain [" << m << "] proposed model REJECTED! comparator=" << *comparator << "  r=" << r << std::endl;
		(*model_current).moved[m]=0;
	}
	(*model_current).Pmove[m]=r; // In all case we save the move probability

}

void MALA::execute(Model_def *model_current, Model_def *model_propose, Data *data_struc, Config *cfg_class, Outputs *outputs_class){
/*
 * Main function that execute the process according to a given configuration.
 * The class MALA must have been fully initialized before proceeding.
 * The model_current must also have been fully initialized.
 * Content of the passed classes:
 *     - model_current: Fully initialized model. This includes priors
 *     - model_propose: Fully initialized model. The configuration must be the same as model_current
 *     - data_struc: A structure with the observational data
 *     - cfg_class: a class which contains the full configuration. Typically, this
 *                  would have been initialized by reading a configuration file.
 *                  Note however that for the moment, we only use the string 
 *                 cfg_class.proposal_type in MALA::update_position_MH().
*/

	bool logic; // a boolean that tells when we learn the var-covariance and when we don't
	bool tempted_mixing; // 1 if we fulfill the condition for parallel-chains mixing
	int which;
	VectorXd vec_tmp; // required to let the compiler happy when passing (*model_current).vars.row(chain) into MALA::update_proposal()
	int chain_A;


	long i=0;
	while (i<Nsamples){
		if (i%1000 == 0){ std::cout << "[" << i << "]" << std::endl;}
		gamma=c0/(1. + i);
		for (int chain=0; chain<Nchains; chain++){

			// [1] Propose a new position and test it so that the model_current is updated (if new position accepted) using model_propose
			//     Note that model_propose must have been initialized in order to work
			update_position_MH(model_current, model_propose, data_struc, cfg_class, chain);


			// [2] Test whether we update the variance-covariance parameters
			logic=0; // reset
			//which=0; 
			for(int l=0; l<periods_learn.size(); l++){ // note that periods.size() must be Nt_learn.size() - 1
				logic= logic || ((i>=Nt_learn[l]) && (i<Nt_learn[l+1])); //Test whether we have to learn
				if ((i>=Nt_learn[l]) && (i<Nt_learn[l+1])){ // Tells which phase we are actually doing
					which=l;
				}
			}
			if (logic == 1 && ( i%periods_learn[which]) == 0){
				vec_tmp=(*model_current).vars.row(chain);
				update_proposal(vec_tmp, (*model_current).Pmove[chain], chain);
				// ---------- For Debug only -----------
				/*
				std::cout << " Learning chain " << chain << "... " << std::endl;
				std::cout << "(*model_current).vars.row(" << chain << ") = " << (*model_current).vars.row(chain) << endl;
				std::cout << "sigma[" << chain << "] = " << sigma[chain] << endl;
				std::cout << "covarmat[" << chain << "] = " << endl;
				std::cout << *covarmat[chain] << endl;
				*/				
				// -------------------------------------				
				
			} //else { std::cout << " NO Learning for chain " << chain << std::endl; }
		}
		// [3] Parallel tempering
		tempted_mixing=0;
		chain_A=-1;
		if (i%dN_mixing == 0){
			//std::cout << "Attempting Mixing..." << std::endl;
			chain_A=parallel_tempering(model_current);
			tempted_mixing=1;
		}
		// [4] Show / save results
		//          A lot of implementation still required here!
		//std::cout << "[" << i << "]" << std::endl;
		//std::cout << "    - current parameters of chain 0" << std::endl;
		//std::cout << (*model_current).params.row(0) << std::endl;
		//std::cout << "    - proposed parameters chain 0" << std::endl;
		//std::cout << (*model_propose).params.row(0) << std::endl;
		//std::cout << "------------------------" << std::endl;

		outputs_class->update_buffer_stat_criteria(model_current->logLikelihood, model_current->logPrior, model_current->logPosterior); 	
		outputs_class->update_buffer_params(model_current->vars);
		outputs_class->update_buffer_ptempering(tempted_mixing, chain_A, model_current->Pswap, model_current->swaped);
		outputs_class->update_buffer_proposals(sigma, mu, model_current->Pmove, 
			    			model_current->moved, covarmat);
		outputs_class->update_buffer_models(model_current->model); //IF ACTIVATED IT SLOWS DOWN THE PROCESS AND CREATE HUGE FILES!!!
		
		i=i+1;
	}


}


int main(){
/*
	long Nch=3;
	VectorXd vars(4), vars_new(4), relax(4);
	vector<string> vars_names(4);
	vars_names.assign("");
	MALA params_algo(1000, 3, vars, vars_names, relax);
	
	vars(0)=0.25;
	vars(1)=0.5;
	vars(2)=0.75;
	vars(3)=1.;

	std::cout << "Here" << std::endl;
	var_names.assign(4, "");
	var_names.at(0)="Height l=0";
	var_names.at(1)="Width l=0";
	var_names.at(2)="Frequency l=0";
	var_names.at(3)="White Noise N0";
	std::cout << "-----" << std::endl;

	std::cout << " --------- Individual tests ----------" << std::endl;
	std::cout << " [1] MALA::init_proposal " << std::endl;
	//params_algo(1000, 3, vars.size());
	params_algo.init_proposal(vars, var_names);
	//for(int m=0; m<Nch; m++){
	//	std::cout << *params_algo.covarmat[m] << std::endl;
	//}
	
	std::cout << " [2] MALA::new_prop_values and MALA::update_proposal " << std::endl;
	double acceptance=0.23;
	int m=0;
	vars_new=params_algo.new_prop_values(vars, m);
	params_algo.update_proposal(vars_new, acceptance, m);
*/
// -------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------
	long Nchains=6;

	std::cout << " --------- FULL MALA TEST ----------" << std::endl;
	string model_fct_name, likelihood_fct_name, prior_fct_name;
	vector<string> params_names, priors_params_names;
	VectorXd params, priors_params;
	VectorXi relax, plength;

	params.resize(4);
	params(0)=1.;
	params(1)=0.5;
	params(2)=0.;
	params(3)=0.3;

	params_names.assign(4, "");
	params_names.at(0)="Height_l0";
	params_names.at(1)="Width_l0";
	params_names.at(2)="Frequency_l0";
	params_names.at(3)="White_Noise_N0";

	priors_params.setZero(4); // 4 priors of params are set to 0
	
	priors_params_names.assign(4, "");
	priors_params_names[0]="Height_Prior";
	priors_params_names[1]="Width_Prior";
	priors_params_names[2]="Frequency_Prior";
	priors_params_names[3]="White_Noise_Prior";

	plength.resize(2);
	plength[0]=3; // Parameters of the Gaussian
	plength[1]=1; // Noise parameters

	relax.resize(4);
	relax << 1, 1, 1, 1;

	model_fct_name="model_Test_Gaussian";
	likelihood_fct_name="chi_square";
	prior_fct_name="prior_Test_Gaussian";

	// Load the configuration
	//string file_in_data="/home/obenomar/Dropbox/Temporary/Cpp-Playground/data_tests/test-gauss.txt";
	string file_in_data="/Users/obenomar/Dropbox/Temporary/Cpp-Playground/data_tests/test-gauss.txt";
    
    Config config(file_in_data, "","");

	// Initialise the class of the MALA algorithm
	int Nsamples=5000;
	MALA TAMCMC(Nsamples, Nchains, params, params_names, relax); // Nsamples, Nchains, Nvars

	// Initialise the classes model_propose and model_current
	// This is necessary because we want to avoid to have to reset all the 
	// class in the loop
	Model_def model_current(model_fct_name, likelihood_fct_name, prior_fct_name,
		relax, plength, params, params_names, 
		priors_params, priors_params_names, Nchains, &config, TAMCMC.Tcoefs);  // last we put the number of point in the model to initialize model_current.model
	Model_def model_propose(model_fct_name, likelihood_fct_name, prior_fct_name,
		relax, plength, params, params_names, 
		priors_params, priors_params_names, Nchains, &config, TAMCMC.Tcoefs); 

/*	std::cout << " Does models are correctly initialized? " << std::endl;

	for (int m=0; m<Nchains;m++){
		std::cout << " Information for debug " << std::endl;
		std::cout << "     - Model current: " << std::endl;
		std::cout << " 			[1] params.row(" << m << ")=" << model_current.params.row(m) << std::endl;
		std::cout << " 			[2] vars.row(" << m << ")=" << model_current.vars.row(m) << std::endl;
		//std::cout << " 			[3] model.row(" << m << ")=" << model_current.model.row(m) << std::endl;
		std::cout << " 			[4] logLikelihood[" << m << "]=" << model_current.logLikelihood[m] << std::endl;
		std::cout << " 			[5] logPrior[" << m << "]=" << model_current.logPrior[m] << std::endl;
		std::cout << " 			[6] logPosterior[" << m << "]=" << model_current.logPosterior[m] << std::endl;
		std::cout << " 			[7] moved[" << m << "]=" << model_current.moved[m] << std::endl;
		std::cout << "     - Model propose: " << std::endl;
		std::cout << " 			[1] params.row(" << m << ")=" << model_propose.params.row(m) << std::endl;
		std::cout << " 			[2] vars.row(" << m << ")=" << model_propose.vars.row(m) << std::endl;
		//std::cout << " 			[3] model.row(" << m << ")=" << model_propose.model.row(m) << std::endl;
		std::cout << " 			[4] logLikelihood[" << m << "]=" << model_propose.logLikelihood[m] << std::endl;
		std::cout << " 			[5] logPrior[" << m << "]=" << model_propose.logPrior[m] << std::endl;
		std::cout << " 			[6] logPosterior[" << m << "]=" << model_propose.logPosterior[m] << std::endl;
		std::cout << " 			[7] moved[" << m << "]=" << model_propose.moved[m] << std::endl;

		std::cout << " ------------------------" << std::endl;
	}			
*/
//	// Initialise the classes 'out' of outputs in files
	string diag_cfg_file="";
	Outputs out(diag_cfg_file, Nchains, config.data.data.Nx, Nsamples, params, params_names, relax, plength, TAMCMC.Tcoefs);

	//------------------------------
	clock_t begin_time = clock();
	TAMCMC.execute(&model_current, &model_propose, &config.data.data, &config, &out);
	std::cout << "    Calculation finished in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;


}



