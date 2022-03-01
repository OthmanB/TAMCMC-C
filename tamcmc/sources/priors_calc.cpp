/*
 * priors_calc.cpp
 *
 *  Contains the function that calculate the penalization
 *  according to presets. So far, the only priors preset is for
 *  a main sequence star: priors_MS_Global()
 * 
 *  Created on: 02 Mar 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "math.h"
#include "priors_calc.h"
#include "stats_dictionary.h"
#include "derivatives_handler.h"
#include "linfit.h"
#include "linspace.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

long double priors_MS_Global(const VectorXd& params, const VectorXi& params_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch, const VectorXd& extra_priors){

	long double f=0;

	const int smooth_switch=extra_priors[0];
	const double scoef=extra_priors[1];
	//const double a3ova1_limit=extra_priors[2];
	//const double a2ova1_limit=extra_priors[4]; // limit for a2 when this is in use
	//const int impose_normHnlm=extra_priors[3];
	//const int model_index=extra_priors[5]; // An index specifying specific rules that might be applied for specific models
	//const int numax=extra_priors[3]; 'Prior on numax, only applied if non-zero Not used.
	VectorXd ajova1_limit(6);
	ajova1_limit << extra_priors[2], extra_priors[3],extra_priors[4],extra_priors[5], extra_priors[6], extra_priors[7];
	const int impose_normHnlm=extra_priors[8];
	const int model_index=extra_priors[9]; // An index specifying specific rules that might be applied for specific models
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 3 for a global MS model (a1,a2,a3)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3; // Total number of modes
	int i0, j;
	double fl, fl1, fl2,fl3;
	double Dnu, d02, aj, a1, a2, a3, alfa, b, fmax, Q11, max_b, el, em;
	VectorXd tmp, fit, a2_terms, epsilon_terms, thetas, Nfl(4), a1_n;
	Deriv_out frstder, scdder;

	Nfl << Nfl0, Nfl1, Nfl2, Nfl3;

//	std::cout << "[0] f=" << f << std::endl;

	// ----- Add a positivity condition on visibilities ------
	for(int i=Nmax; i<=Nmax+lmax; i++){
		if(params[i] < 0){
			f=-INFINITY;
			goto end;
		}
	}
//	std::cout << "[1] f=" << f << std::endl;
	// Specific cases
	switch (model_index){
		case 0: // *.a1etaa3 family of models
			a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);
			a3=params[Nmax+lmax+Nf+2];
			break;
		case 1: // *a1a2a3 family 
			a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);
			a2_terms=params.segment(Nmax + lmax + Nf + 6,3);
			a3=params[Nmax+lmax+Nf+2];
			//std::cout << "[In priors] a2_terms = " << a2_terms.transpose() << std::endl;
			for (int n=0; n<Nfl1; n++){
				fl1=params[Nmax+lmax+Nfl0+n];
				a2=a2_terms[0] + a2_terms[1]*(fl1*1e-3) + a2_terms[2]*(fl1*fl1*1e-6);
				//std::cout << "      abs(a2/a1) = " << std::abs(a2/a1)  << "   " << (std::abs(a2/a1) >= a2ova1_limit) << std::endl;
				if (std::abs(a2/a1) >= ajova1_limit[1]){
					f=-INFINITY;
					goto end;
				}
			}
			//std::cout << " f after l=1:  " << f << std::endl;
			for (int n=0; n<Nfl2; n++){
				fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
				a2=a2_terms[0] + a2_terms[1]*(fl2*1e-3) + a2_terms[2]*(fl2*fl2*1e-6);
				//std::cout << "      abs(a2/a1) = " << std::abs(a2/a1)  << "   " << (std::abs(a2/a1) >= a2ova1_limit) << std::endl;
				if (std::abs(a2/a1) >= ajova1_limit[1]){
					f=-INFINITY;
					goto end;
				}
			}
			//std::cout << " f after l=2:  " << f << std::endl;
			for (int n=0; n<Nfl3; n++){
				fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
				a2=a2_terms[0] + a2_terms[1]*(fl3*1e-3) + a2_terms[2]*(fl3*fl3*1e-6);
				//std::cout << "      abs(a2/a1) = " << std::abs(a2/a1)  << "   " << (std::abs(a2/a1) >= a2ova1_limit) << std::endl;
				if (std::abs(a2/a1) >= ajova1_limit[1]){
					f=-INFINITY;
					goto end;
				}
			}
 		   // Prior on a3/a1 ratio. a3 << a1 is enforced here by putting a3ova1_limit
			if(std::abs(a3/a1) >= ajova1_limit[2]){
				f=-INFINITY;
				goto end;
			}
			//std::cout << " f after l=3:  " << f << std::endl;
			break;
		case 2: //model_MS_Global_a1l_etaa3_HarveyLike
			a1=(std::abs(params[Nmax + lmax + Nf]) + std::abs(params[Nmax + lmax + Nf+6]))/2; // a1_avg = (a11 + a12)/2
			a3=params[Nmax+lmax+Nf+2];
			std::cout << "in priors_calc.cpp: case 2 needs checks. Exiting" << std::endl;
			exit(EXIT_SUCCESS);
			break;
		case 3: // model_MS_Global_a1n_etaa3_HarveyLike
			a1=params.segment(Nmax + lmax + Nf + 6, Nmax).sum()/Nmax; // a1_avg = avg(a12) = avg(a11) 
			a3=params[Nmax+lmax+Nf+2];
			std::cout << "in priors_calc.cpp: case 3 needs checks. Exiting" << std::endl;
			exit(EXIT_SUCCESS);
			break;
		case 4: // model_MS_Global_a1nl_etaa3_HarveyLike
			a1=params.segment(Nmax + lmax + Nf + 6, 2*Nmax).sum()/(2.*Nmax); // a1_avg = avg(a11, a12) with size(a11)=size(a12) = Nmax 
			a3=params[Nmax+lmax+Nf+2];
			std::cout << "in priors_calc.cpp: case 4 needs checks. Exiting" << std::endl;
			exit(EXIT_SUCCESS);
			break;
		case 5:  // model_MS_Global_a1n_a2a3_HarveyLike
			a1=params.segment(Nmax + lmax + Nf + 6, Nmax).sum()/Nmax; // a1_avg = avg(a12) = avg(a11) 
			a2_terms=params.segment(Nmax + lmax + Nf + 6 + Nmax, Nmax);
			a3=params[Nmax+lmax+Nf+2];
			std::cout << "in priors_calc.cpp: case 5 needs checks. Exiting" << std::endl;
			exit(EXIT_SUCCESS);
			break;
		case 6: //model_MS_Global_a1l_a2a3_HarveyLike
			 a1=(std::abs(params[Nmax + lmax + Nf]) + std::abs(params[Nmax + lmax + Nf+6]))/2; // a1_avg = (a11 + a12)/2
			 a2_terms=params.segment(Nmax + lmax + Nf + 7, lmax*3 );
			 a3=params[Nmax+lmax+Nf+2];
			std::cout << "in priors_calc.cpp: case 6 needs checks. Exiting" << std::endl;
			 exit(EXIT_SUCCESS);
			 break;
		case 7: // model_MS_Global_a1nl_a2a3_HarveyLike
			a1=params.segment(Nmax + lmax + Nf + 6, 2*Nmax).sum()/(2.*Nmax); // a1_avg = avg(a11, a12) with size(a11)=size(a12) = Nmax 
			a2_terms=params.segment(Nmax + lmax + Nf + 6 + Nmax + Nmax , Nmax); // after a11 and a12
			a3=params[Nmax+lmax+Nf+2];
			std::cout << "in priors_calc.cpp: case 7 needs checks. Exiting" << std::endl;
			exit(EXIT_SUCCESS);
			break;
		case 8: // model_MS_Global_a1etaGlma3_HarveyLike
			a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);
			a3=params[Nmax+lmax+Nf+2];
 		   // Prior on a3/a1 ratio. a3 << a1 is enforced here by putting a3ova1_limit
			if(std::abs(a3/a1) >= ajova1_limit[2]){
				f=-INFINITY;
				goto end;
			}
			//epsilon_terms=params.segment(Nmax + lmax + Nf + 6,3);
			//thetas=params.segment(Nmax+lmax+Nf+6+3,2);
			//std::cout << "[In priors] epsilon_terms = " << epsilon_terms.transpose() << std::endl;
			//std::cout << "[In priors] thetas = " << thetas.transpose() << std::endl;
			/*if (thetas[0]>thetas[1]){ // Impose theta_min <= theta_max
					std::cout << " theta condition issue: delta=" << thetas[1]-thetas[0] << std::endl;
					f=-INFINITY;
					goto end;				
			}
			*/
			break;
		case 9:
			// Ensure that none of the aj coefficient exceeds the set limits at any given frequency
			i0=Nfl[0]; // We start at l=1... need to get the pointer to Nmax+lmax+Nfl0 then
			//std::cout << " Nfl = " << Nfl.transpose() << std::endl;
			for (int el=1; el<lmax+1; el++){
				for (int j=1; j<6;j++){ // compare a2, a3, a4, a5, a6,to a1 >> set the limit
					for (int n=0; n<Nfl[el]; n++){
						fl=params[Nmax+lmax+i0+n];
						a1=params[Nmax + lmax + Nf] + params[Nmax+lmax+Nf+1]*(fl*1e-3);
						aj=params[Nmax + lmax + Nf+ 2*j] + params[Nmax+lmax+Nf+ 2*j+1]*(fl*1e-3);
						//std::cout << "  f(l=" << el << ", n=" << n << ")  =" << fl << std::endl;
						//std::cout << "  a1(l=" << el << ", n=" << n << ") ="<< a1 << std::endl;
						//std::cout << "  a" << j+1 << "(l=" << el << ", n=" << n << ") = " << aj << std::endl;
						//std::cout << "         ajova1_limit[" << j << "] = " << ajova1_limit[j] << std::endl;
						//std::cout << "     ----- " << std::endl;
						if (std::abs(aj/a1) >= ajova1_limit[j]){
							f=-INFINITY;
							goto end;
						}
					}
				}
				i0=i0+Nfl[el];
				//std::cout << " el = " << el << "       Nmax+lmax+i0 = " << Nmax+lmax+i0 << std::endl;
				//std::cout << "     ----- " << std::endl;
			}
		default:
			a1=params[Nmax+lmax+Nf]; // By default it is the slot reserved for a1 (Classic models)
			switch(impose_normHnlm){ 
				case 1: // Case specific to model_MS_Global_a1etaa3_HarveyLike_Classic_v2
					//l=1: m=0, m=+/-1
					f=f+logP_uniform(0, 1.+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise] + 2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+1]); // The sum must be positive
					//l=2: m=0, m=+/-1, m=+/-2
 					// The sum must be positive
					f=f+logP_uniform(0, 1+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+2]+ 
										   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+3]+ 
										   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+4]
									); // The sum must be positive
					//l=3: m=0, m=+/-1, m=+/-2, m=+/-3
 					// The sum must be positive
					f=f+logP_uniform(0, 1+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+5]+ 
										   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+6]+ 
										   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+7]+
										   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+8]
										   ); // The sum must be positive
		
				break;
				case 2:
					std::cout << "priors_MS_Global: impose_normHnlm=2 YET TO BE IMPLEMENTED!" << std::endl;
					exit(EXIT_SUCCESS);
				break;
			}
			break;
	}
//	std::cout << "[2] f=" << f << std::endl;

	// Implement securities to avoid unphysical quantities that might lead to NaNs
	if(params[Nmax+lmax+Nf+4] < 0){ // Impose that the power coeficient of magnetic effect is positive
		f=-INFINITY;
		goto end;
	}
	if(priors_names_switch[Nmax+lmax+Nf+Nsplit+Nwidth+3] != 0){
		if( (params[Nmax+lmax+Nf+Nsplit+Nwidth+3] < 0) || // Harvey profile height
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+4] < 0) || // Harvey profile tc
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+5] < 0) ){ // Harvey profile p
				f=-INFINITY;
				goto end;
		}
	}
	if(priors_names_switch[Nmax+lmax+Nf+Nsplit+Nwidth+6] != 0){
		if( (params[Nmax+lmax+Nf+Nsplit+Nwidth+6] < 0) || // Harvey profile height
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+7] < 0) || // Harvey profile tc
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+8] < 0) ){// Harvey profile p 
				f=-INFINITY;
				goto end;
		}
	}
	if((priors_names_switch[Nmax+lmax+Nf+9] != 0) && (params[Nmax+lmax+Nf+Nsplit+Nwidth+9] < 0)){
		f=-INFINITY;
		goto end;
	}
//	std::cout << "[3]  f=" << f << std::endl;

	// Apply the priors as defined in the configuration defined by the user and read by 'io_MS_global.cpp'
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);
//	std::cout << "[4] f=" << f << std::endl;

	// Determine the large separation
	//frstder=Frstder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // First derivative of fl0 gives Dnu
	//Dnu=frstder.deriv.sum();
	tmp=linspace(0, params.segment(Nmax+lmax, Nfl[0]).size()-1, params.segment(Nmax+lmax, Nfl[0]).size());
	fit=linfit(tmp, params.segment(Nmax+lmax, Nfl[0])); // fit[0] is the slope ==> Dnu and fit[1] is the ordinate at origin ==> fit[1]/fit[0] = epsilon
	Dnu=fit[0];

//	std::cout << "[5] f=" << f << std::endl;

	// Apply a prior on the d02
	//std::cout << "--- d02 --" << std::endl;
	if(Nfl[0] == Nfl[2]){
		for(int i=0; i<Nfl[0]; i++){
			d02=params[Nmax+lmax+i] - params[Nmax+lmax+Nfl[0]+Nfl[1]+i];
			f=f+logP_gaussian_uniform( 0, Dnu/3., 0.015*Dnu, d02); // This is mainly for F stars
		}
	}
//	std::cout << "[6] f=" << f << std::endl;

	// Set the smootheness condition handled by derivatives_handler.cpp
	switch(smooth_switch){
			case 1: // Case with smoothness
				i0=0;
				for (int el=0; el<lmax+1; el++){
					if(Nfl[el] != 0){
						scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+i0, Nfl[el])); // The l=0 frequencies
						for(int i=0; i<Nfl[el]; i++){
							f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
						}
					}
					i0=i0+Nfl[el];			
				}
			  	break;
	}
	end:	
//	std::cout << " Final f:  "<< f << std::endl;
//	exit(EXIT_SUCCESS);
return f;
} 

long double priors_asymptotic(const VectorXd& params, const VectorXi& params_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch, const VectorXd& extra_priors){
	// The priors_asymptotic() function is basically the same as the prior_ms_global() but:
	//    (1) Need to exclude l=1 from the smoothing
	//    (2) Replace a3/a1 ratio by a3/rot_env (A slight change for the a3 index because we have Snlm= [rot_env, rot_core, eta, a3, asym] here instead of [a1, eta, a3, asym]
	
	long double f=0;

	const int smooth_switch=extra_priors[0];
	const double scoef=extra_priors[1];
	const double a3ova1_limit=extra_priors[2];
	const int impose_normHnlm=extra_priors[3];
	const int model_switch=extra_priors[4];  // Switch that specify the model that is used... its definition is made inside io_asymptotic.cpp


	//const int numax=extra_priors[3]; 'Prior on numax, only applied if non-zero Not used.
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of parameters for l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 3 for a global MS model (rot_env, rot_core, eta, a3, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3; // Total number of modes
	const int Nmixedmodes_g_params=7; // THIS IS SOMEWHAT PART OF NFL1. BUT WE COUNT HERE ONLY PARAMETERS FOR G MODES: DP, q, alpha, delta01, sigma_Hl1, sigma_fl1g, sigma_fl1m

	int i, i0=Nmax+lmax+Nf+Nsplit;

	double Dnu, d02, rot_env, a3, alfa, b, fmax, Q11, max_b, el, em;
	double Dnu_l1;
	VectorXd tmp, fit;
	Deriv_out frstder, scdder;
	
	a3=params[Nmax+lmax+Nf+3];
	rot_env=std::abs(params[Nmax + lmax + Nf]);

	// ----- Add a positivity condition on visibilities ------
	for(int i=Nmax; i<=Nmax+lmax; i++){
		if(params[i] < 0){
			f=-INFINITY; 
			goto end;
		}
	}

	// Prior on a3/a1 ratio. a3 << a1 is enforced here by putting a3ova1_limit
	if(std::abs(a3/rot_env) >= a3ova1_limit){
		f=-INFINITY;
		goto end;
	}
	// Implement securities to avoid unphysical quantities that might lead to NaNs
	if(params[Nmax+lmax+Nf+4] < 0){ // Impose that the power coeficient of magnetic effect is positive
		f=-INFINITY;
		goto end;
	}
	if(priors_names_switch[Nmax+lmax+Nf+Nsplit+Nwidth+3] != 0){
		if( (params[Nmax+lmax+Nf+Nsplit+Nwidth+3] < 0) || // Harvey profile height
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+4] < 0) || // Harvey profile tc
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+5] < 0) ){ // Harvey profile p
				f=-INFINITY;
				goto end;
		}
	}
	if(priors_names_switch[Nmax+lmax+Nf+Nsplit+Nwidth+6] != 0){
		if( (params[Nmax+lmax+Nf+Nsplit+Nwidth+6] < 0) || // Harvey profile height
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+7] < 0) || // Harvey profile tc
		    (params[Nmax+lmax+Nf+Nsplit+Nwidth+8] < 0) ){// Harvey profile p 
				f=-INFINITY;
				goto end;
		}
	}
	if((priors_names_switch[Nmax+lmax+Nf+9] != 0) && (params[Nmax+lmax+Nf+Nsplit+Nwidth+9] < 0)){
		f=-INFINITY;
		goto end;
	}
	// ------ Add a positiviy condition on Width parameters -----
	for (i=i0; i<i0 + Nwidth; i++){
		switch(priors_names_switch[i]){
			case 2:
				if (params[i] < 0){
					f=-INFINITY;
					goto end;
				}
			break;
		}
	}

	// Apply the priors as defined in the configuration defined by the user and read by 'io_MS_global.cpp'
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);

	// ----- Add a positivity condition on inclination -------
	// The prior could return values -90<i<90. We want it to give only 0<i<90
	//f=f+logP_uniform(0., 90., params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise]);
	/*switch(impose_normHnlm){ 
		case 1: // Case specific to model_MS_Global_a1etaa3_HarveyLike_Classic_v2
			//l=1: m=0, m=+/-1
			f=f+logP_uniform(0, 1.+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise] + 2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+1]); // The sum must be positive
			//l=2: m=0, m=+/-1, m=+/-2
 			// The sum must be positive
			f=f+logP_uniform(0, 1+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+2]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+3]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+4]
							); // The sum must be positive
			//l=3: m=0, m=+/-1, m=+/-2, m=+/-3
 			// The sum must be positive
			f=f+logP_uniform(0, 1+1e-10, params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+5]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+6]+ 
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+7]+
								   2*params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+8]
								   ); // The sum must be positive

		break;
		case 2:
			std::cout << "priors_MS_Global: impose_normHnlm=2 YET TO BE IMPLEMENTED!" << std::endl;
			exit(EXIT_SUCCESS);
		break;
	}
	*/
	// Determine the large separation
	//frstder=Frstder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // First derivative of fl0 gives Dnu
	//Dnu=frstder.deriv.sum();

	tmp=linspace(0, params.segment(Nmax+lmax, Nfl0).size()-1, params.segment(Nmax+lmax, Nfl0).size());
	fit=linfit(tmp, params.segment(Nmax+lmax, Nfl0)); // fit[0] is the slope ==> Dnu and fit[1] is the ordinate at origin ==> fit[1]/fit[0] = epsilon
	Dnu=fit[0];
	//epsilon=fit[1]/fit[0];
	//epsilon=epsilon - floor(epsilon);
	
	switch(model_switch){ // SPECIFIC HANDLING FOR SOME MODELS
		case 1: // Case of model model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3 or model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3
			//std::cout << "fl1p:" << params.segment(Nmax+lmax+Nfl0+Nmixedmodes_g_params, Nfl1-Nmixedmodes_g_params).transpose() << std::endl;
			tmp=linspace(0, params.segment(Nmax+lmax+Nfl0+Nmixedmodes_g_params, Nfl1-Nmixedmodes_g_params).size()-1, params.segment(Nmax+lmax+Nfl0+Nmixedmodes_g_params, Nfl1-Nmixedmodes_g_params).size());
			fit=linfit(tmp, params.segment(Nmax+lmax+Nfl0+Nmixedmodes_g_params, Nfl1-Nmixedmodes_g_params)); // fit[0] is the slope ==> Dnu and fit[1] is the ordinate at origin ==> fit[1]/fit[0] = epsilon
			Dnu_l1=fit[0];
			f=f+ logP_gaussian(Dnu, 0.003*Dnu,Dnu_l1);  // IMPOSES Dnu(l=0) = Dnu(l=1) + N(0, 0.01*Dnu(l=0))
			if(Nfl1-Nmixedmodes_g_params >= 3){ // We need suficient number of g modes to apply a smoothness condition
				scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nmixedmodes_g_params, Nfl1-Nmixedmodes_g_params)); // The l=0 frequencies
				for(int i=0; i<scdder.deriv.size(); i++){
					f=f+ logP_gaussian(0, 0.05*Dnu,scdder.deriv[i]); // Penalize the value to enforce the smoothness
				}				
			}	
 		break;
	}

	// Apply a prior on the d02 : NOTE CANNOT DUE TO POTENTIAL l=2 MIXED MODES
	//std::cout << "--- d02 --" << std::endl;
	/*if(Nfl0 == Nfl2){
		for(int i=0; i<Nfl0; i++){
			d02=params[Nmax+lmax+i] - params[Nmax+lmax+Nfl0+Nfl1+i];
			f=f+logP_gaussian_uniform(0, Dnu/3., 0.015*Dnu, d02); // This is mainly for F stars
		}
	}
	*/
	// Set the smootheness condition handled by derivatives_handler.cpp
	switch(priors_names_switch[i]){
			case 1: // Case with smoothness
				//std::cout << " ------- Frequency derivatives ------" << std::endl;	
				if(Nfl0 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax, Nfl0)); // The l=0 frequencies
				}
				//std::cout << "--- Fl0 --" << std::endl;
				for(int i=0; i<Nfl0; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
				}

				// -- NOT APPLICABLE DUE TO MIXED MODES --
				//if(Nfl1 != 0){
				//	scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0, Nfl1)); // The l=1 frequencies
				//}
				//std::cout << "--- Fl1 --" << std::endl;
				//for(int i=0; i<Nfl1; i++){
				//	f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
				//}
	
				/* NOTE: CANNOT DUE POTENTIAL l=2 MIXED MODES... MAY NEED A SMARTER SWITCH
 				if(Nfl2 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1, Nfl2)); // The l=2 frequencies
				}
				
				//std::cout << "--- Fl2 --" << std::endl;
				for(int i=0; i<Nfl2; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
	
				}
				*/
				if(Nfl3 != 0){
					scdder=Scndder_adaptive_reggrid(params.segment(Nmax+lmax+Nfl0+Nfl1+Nfl2, Nfl3)); // The l=3 frequencies
				}
				//std::cout << "--- Fl3 --" << std::endl;
				for(int i=0; i<Nfl3; i++){
					f=f+ logP_gaussian(0, scoef,scdder.deriv[i]); // Penalize the value
				}
			  	break;
	}
	/*
	std::cout << "a3 =" << a3 << std::endl;
	std::cout << "rot_env ="	<< params[Nmax+lmax+Nf] << std::endl;
	std::cout << "ratio =" << a3/rot_env<< std::endl;
	std::cout << "a3ova1_limit = " << a3ova1_limit << std::endl; 
	std::cout << "after a3ova1_limit " << f << std::endl;
	std::cout << "extra_priors = " << extra_priors << std::endl;
	std::cout << " MODIFICATION MADE ON 21 JAN 2021 (rot_env ~ a1 to get a3/a1... Need checks. Exiting" << std::endl;
	exit(EXIT_SUCCESS);
	*/
	end:

	return f;
}

long double priors_local(const VectorXd& params, const VectorXi& params_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch, const VectorXd& extra_priors){

	long double f=0;

	//const int smooth_switch=extra_priors[0];
	const double a3ova1_limit=extra_priors[2];
	const int impose_normHnlm=extra_priors[3];
	//const int numax=extra_priors[3]; 'Prior on numax, only applied if non-zero Not used.
	
	const int Nmax=params_length[0]; // Number of Heights
	const int Nvis=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 3 for a global MS model (a1,a2,a3)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
	const int lmax=3; // lmax is hardcoded here

	int pos0;
	VectorXi Nf_el(4);  
	Nf_el[0]=Nfl0; Nf_el[1]=Nfl1; Nf_el[2]=Nfl2; Nf_el[3]=Nfl3; 

	// Prior on a3/a1 ratio. a3 << a1 is enforced here by putting a3ova1_limit
	if (params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3] !=0){ // If there is a a1 to be fitted
		if(std::abs(params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+2]/params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3]) >= a3ova1_limit){
			f=-INFINITY;
			goto end;
		}
	} else{
		if ((params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+3] !=0) &&  (params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+4] !=0)){ // If the sqrt(a1).cosi and sqrt(a1).sin(i) are not 0
			//std::cout << "Case sqrt(a1).cos and sin" << std::endl;
			if(std::abs(params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+2]/(pow(params[Nmax + Nvis + Nfl0+Nfl1+Nfl2+Nfl3+3],2)+ pow(params[Nmax + Nvis + Nfl0+Nfl1+Nfl2+Nfl3+4],2))) >= a3ova1_limit){
				f=-INFINITY;
				goto end;
			}
		}
	}
	//std::cout << f << std::endl;
	
	// Prior in inclination positivity (if relevant, ie if it is directly fitted)
	//std::cout << " priors_names_switch[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] = " << priors_names_switch[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] << std::endl;
	//std::cout << " params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] =" << params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] << std::endl;
	if((priors_names_switch[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] != 0) && (params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise] < 0)){
		f=-INFINITY;
		goto end;
	}

	//double Dnu, d02, scoef, a1, alfa, b, fmax, Q11, max_b, el, em;
	
 	// Apply the priors as defined in the configuration defined by the user and read by 'io_MS_global.cpp'
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);
	// ----- Add a positivity condition on inclination -------
	// The prior could return values -90<i<90. We want it to give only 0<i<90
	//f=f+logP_uniform(0., 90., params[Nmax+Nvis+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+Nnoise]);


// ALL THIS IS NOT NECESSARY BECAUSE WE DO NOT DEAL WITH VISIBILITIES... NO NORMALISATION TO 1 REQUIRED
// THEREFORE NO NEED TO EXCLUDE VALUES BEYOND 1.
/*	switch(impose_normHnlm){ 
		case 1: // Case specific to model_MS_Global_a1etaa3_HarveyLike_Classic_v2
			//l=1: m=0, m=+/-1
			std::cout << "priors_MS_local: impose_normHnlm=1 Not required!" << std::endl;
			std::cout << "                 Likely bug in io_local.cpp ... Exiting now" << std::endl; 
			exit(EXIT_SUCCESS);	
			break;
		case 2:
			std::cout << "In case 2" << std::endl;
			std::cout << "lmax=" << lmax << std::endl;
			for(int el=1; el<=lmax; el++){
					for(int en=0; en<Nf_el[el]; en++){
						for(int em=0; em<=el; em++){
							std::cout << "el=" << el << "   ==> " <<  Nf_el.segment(0,el+1)  << std::endl;
							pos0=Nf_el.segment(0,el+1).sum() + (el+1)*en;
							std::cout << pos0 << std::endl;
							switch(el){
								case 1: //l=1: m=0, m=+/-1
									std::cout << "In l=1: m=0, m=+/-1" << std::endl;
									f=f+logP_uniform(0, 1, params[pos0] + 2*params[pos0+1]); // The sum must be positive
									break;
								case 2: //l=2: m=0, m=+/-1, m=+/-2
									std::cout << "In l=2: m=0, m=+/-1, m=+/-2" << std::endl;
									f=f+logP_uniform(0, 1, params[pos0] + 2*params[pos0+1] + 2*params[pos0+2]); // The sum must be positive
									break;
								case 3:
									std::cout << "In l=3: m=0, m=+/-1, m=+/-2, m=+/-3" << std::endl;
									f=f+logP_uniform(0, 1, params[pos0]+ 2*params[pos0+1]+ 2*params[pos0+2] + 2*params[pos0+3]); // The sum must be positive
									break;
							}
						}
					}
				}

			break;
	}
*/
	//std::cout << f << std::endl;

	// BELOW: COMMENTED ON 13 AUG 2020 AS THESE LINES MAKE NO SENSE TO ME AND INDUCE A CRASH IN SOME OCCASION. I SUSPECT THAT THIS WAS POSITIVITY CONDITION ON INC
	// WAS REPLACE BY THE LINE ABOVE BUT KEPT HERE IN CASE I FINALY UNDERSTAND THOSE LINES
	//if((priors_names_switch[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+9] != 0) && (params[Nmax+lmax+Nfl0+Nfl1+Nfl2+Nfl3+Nsplit+Nwidth+9] < 0)){
	//	f=-INFINITY;
	//}

	//exit(EXIT_SUCCESS);
	//std::cout << "Priors done" << std::endl;
	end:
return f;
} 







//long double priors_Test_Gaussian(const VectorXd params, const VectorXi param_length, const MatrixXd priors_params, const std::vector<std::string> priors_names){
long double priors_Test_Gaussian(const VectorXd& params, const VectorXi& params_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch){

	long double f=0;

	f=f + apply_generic_priors(params, priors_params, priors_names_switch);

	//// More rigid strategy... How much time is spent on the cases ??
	//for(i=0; i<params.size(); i++){
	//	f=f + logP_uniform(priors_params(0, i), priors_params(1, i), params[i]);
	//}
	//exit(EXIT_SUCCESS);
return f;
} 

//long double priors_Test_Gaussian(const VectorXd params, const VectorXi param_length, const MatrixXd priors_params, const std::vector<std::string> priors_names){
long double priors_Harvey_Gaussian(const VectorXd& params, const VectorXi& params_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch){

	long double f=0;
	// Stello2009 scaling relation between Dnu and numax to remove spurious spiky fits
	const long double beta0=0.263;
	const long double beta1=0.77; 
	const long double Dnu_expected=beta0*std::pow(params[8],beta1);  // param2[8] must be numax here
	//std::cout << f << std::endl;
	//std::cout << "Dnu_expected =" << Dnu_expected << std::endl;
	//std::cout << "params[9] =" << params[9] << std::endl;
	if (params[9] < Dnu_expected/2){
		f=f -INFINITY;  // If the width of gaussian is smaller than 1.5 times the numax, reject the solution
		//std::cout << "f = " << f << std::endl;
		goto end;
	}
	f=f + apply_generic_priors(params, priors_params, priors_names_switch);
	end:

return f;
} 


long double apply_generic_priors(const VectorXd& params, const MatrixXd& priors_params, const VectorXi& priors_names_switch){
/*
 * This function apply the generic priors, that are defined in the configuration file and that have been translated into
 * The constant priors_params (the parameters that define the priors) and into priors_names_switch (translation of 
 * priors_names into integers in order to be handled by the switch / case statement.
 *
*/
	int i;
	long double pena=0;
	for(i=0; i<params.size(); i++){
		switch(priors_names_switch[i]){
			case 0: // No Prior Applied --> The parameter is Fixed or for some reason the prior is set to NONE
			  //std::cout << "CASE 0" << std::endl;
			  //std::cout << "[" << i << "] FIX OR NONE, penalizationa=" << pena << std::endl;
			  break;
			case 1: // Uniform Prior
			  pena=pena + logP_uniform(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "CASE 1" << std::endl;
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;
			  break;
			case 2: // Gaussian Prior
			  pena=pena + logP_gaussian(priors_params(0, i), priors_params(1, i),params[i]);
			  //std::cout << "CASE 2" << std::endl;
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;
			  break;
			case 3: // Multivariate Gaussian Prior
			  std::cout << "The Multivariate Gaussian Prior IS NOT HANDLED" << std::endl;
		      std::cout << "The program will exit now" << std::endl;
			  exit(EXIT_FAILURE);
			  //pena=pena + logP_multivariate_gaussian( priors_params(0, i), Matrix, params[i]);
			  break;
			case 4: // Jeffreys Prior
			  //std::cout << "CASE 4" << std::endl;
			  pena=pena + logP_jeffrey(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;	  
			  break;
			case 5: // Uniform-Gaussian Prior
			  //std::cout << "CASE 5" << std::endl;
			  pena=pena + logP_uniform_gaussian( priors_params(0, i), priors_params(1, i), priors_params(2, i), params[i]); // bmin, bmax, sigma, x
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "  " << "priors_params(2, " << i << ")=" << priors_params(2, i) << std::endl;
			  //std::cout << "    pena=" << pena << std::endl;			  
			  break;
			case 6: // Gaussian-Uniform Prior
			  //std::cout << "CASE 6" << std::endl;
			  pena=pena + logP_gaussian_uniform( priors_params(0, i), priors_params(1, i), priors_params(2, i), params[i]); // bmin, bmax, sigma, x
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "  " << "priors_params(2, " << i << ")=" << priors_params(2, i) << std::endl;
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 7: // Gaussian-Uniform-Gaussian Prior
			  //std::cout << "CASE 7" << std::endl;
			  pena=pena + logP_gaussian_uniform_gaussian( priors_params(0, i), priors_params(1, i), priors_params(2, i),  priors_params(3, i), params[i]); // bmin, bmax, sigma1, sigma2, x
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "  " << "priors_params(2, " << i << ")=" << priors_params(2, i) << std::endl;
			  //std::cout << "  " << "priors_params(3, " << i << ")=" << priors_params(3, i) << std::endl;
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 8: // Uniform Prior on abs(params[i])
			  //std::cout << "CASE 8" << std::endl;
			  pena=pena + logP_uniform_abs(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 9: // Uniform Prior on cos(params[i]*pi/180)
			  //std::cout << "CASE 9" << std::endl;
			  pena=pena + logP_uniform_cos(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i);
			  //std::cout << "    pena=" << pena << std::endl;
			  break;
			case 10: // Jeffreys Prior on abs(param[i])
			  //std::cout << "CASE 4" << std::endl;
			  pena=pena + logP_jeffrey_abs(priors_params(0, i), priors_params(1, i), params[i]);
			  //std::cout << "params[" << i << "]=" <<  params[i] << "  " << "priors_params(0, " << i << ")=" << priors_params(0, i);
			  //std::cout << "  " << "priors_params(1, " << i << ")=" << priors_params(1, i) << std::endl;
			  //std::cout << "    pena=" << pena  << std::endl;	  
			  break;
			case 11: // Case Auto: No Prior Applied here --> The user must set his own prior
			  //std::cout << "CASE 10" << std::endl;
			  //std::cout << "[" << i << "] Auto, pena=" << pena << std::endl;
			  break;
			default:
			  std::cout << " Problem in priors_calc.cpp! " << std::endl;
			  std::cout << " priors_names_switch[" << i << "]=" << priors_names_switch[i] << std::endl;
		          std::cout << " This value is not associated to any known case statement " << std::endl;
			  std::cout << " The program will exit now" << std::endl;
			  exit(EXIT_FAILURE);
		}
	}

return pena;
}


