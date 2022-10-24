/*
 * io_asymptotic.cpp
 *
 * Contains all kind of functions
 * used to arrange/handle the strings for the MS_Global models
 * 
 *  Created on: 27 Nov 2020
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "io_asymptotic.h"
#include "io_models.h"
#include "function_rot.h"
#include "io_ms_global.h" // We include this only because io_asymptotic will use several of its functions such as read_MCMC_file_MS_Global() and set_width_App2016_params_v2()
#include "../../external/ARMM/solver_mm.h" // To be able to use some function that generate fl1 modes for examples

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

// Asymptotic models will use the same MCMC file structure than MS_Global models
MCMC_files read_MCMC_file_asymptotic(const std::string cfg_model_file, const bool verbose){
	return  read_MCMC_file_MS_Global(cfg_model_file, verbose);
}

Input_Data build_init_asymptotic(const MCMC_files inputs_MS_global, const bool verbose, const double resol){

	const long double pi = 3.141592653589793238L;
	const double G=6.667e-8;
	const double Teff_sun= 5777; 
	const double Dnu_sun=135.1;
	const double numax_sun=3150.;
	const double R_sun=6.96342e5; //in km
	const double M_sun=1.98855e30; //in kg
	const double rho_sun=M_sun*1e3/(4*pi*std::pow(R_sun*1e5,3)/3); //in g.cm-3
	const int Nmax_prior_params=4; // The maximum number of parameters for the priors. Should be 4 in all my code

	const double Hmin=1, Hmax=10000; // Define the default lower and upper boundary for the Jeffreys priors applied to heights
	const int Nmixedmodes_g_params=7; // Number of global parameters used to describe the l=1 mixed modes
	const double sigma_limit=inputs_MS_global.Dnu/10.;

	double rho=pow(inputs_MS_global.Dnu/Dnu_sun,2.) * rho_sun;
	double Dnl=0.75, trunc_c=-1;
	double numax=inputs_MS_global.numax;
	double err_numax=inputs_MS_global.err_numax;

	// All Default booleans
	bool do_a11_eq_a12=1, do_avg_a1n=1, do_amp=0;
	bool bool_a1sini=0, bool_a1cosi=0;
	bool status_model = false;
	int lmax, Nmixedmodes_params, en, ind, Ntot, p0, p0_effective, cpt;
	//uint8_t do_width_Appourchaux=0; // We need more than a boolean here, but no need to use a 64 bit signed int
	int do_width_Appourchaux=0, model_type=-1, bias_type=-1;
	double tol=1e-2, tmp;
	VectorXi pos_el, pos_relax0, els_eigen, Nf_el(4), plength;
	VectorXd ratios_l, tmpXd, extra_priors, fl1p_all, fref, ferr;
	std::vector<int> pos_relax;
	std::vector<double> f_inputs, h_inputs, w_inputs, f_priors_min, f_priors_max, f_el;
	std::vector<bool> f_relax, h_relax, w_relax; 
	std::vector<int> rf_el, rw_el, rh_el;
	std::vector<std::string> tmpstr_vec;

	std::string tmpstr_h, tmpstr;
	
	Input_Data Snlm_in, Vis_in, Inc_in, Noise_in, freq_in, height_in, width_in; //, width_App2016_params; // This is by block, each category of parameters		
	Input_Data all_in; // The final structure of parameters, using the standards of my code
	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure

	// Flatening and ordering of all the inputs/relax variables
	lmax=inputs_MS_global.els.maxCoeff();

	// -- Initialisation of structures --
    // --- Look for common instruction That must be run before the setup ---------
	all_in.model_fullname=" "; // Default is an empty string
	//all_in.prior_fullname="prior_MS_Global"; // Default set of prior
    for(int i=0; i<inputs_MS_global.common_names.size(); i++){
        if(inputs_MS_global.common_names[i] == "model_fullname" ){ // This defines if we assume S11=S22 or not (the executed model will be different)
        	all_in.model_fullname=inputs_MS_global.common_names_priors[i];
            if(all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v2" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3" ||
            		"model_RGB_asympt_a1etaa3_AppWidth_HarveyLike" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            	do_width_Appourchaux=2;
//            	if(numax != -9999){ 
            	if (numax <= 0){ // Check that if there is a value, this one is a valid input for numax
            		std::cout << "    The model: " << all_in.model_fullname << " is supposed to have numax for argument" << std::endl;
            		std::cout << "    However, you provided a negative input. Please:" << std::endl;
//            		std::cout << "           [1] Not specify any numax or set !n -9999 in the .model file. In that case, the code will calculate a numax using weighted average of frequencies using Heights as weights" << std::endl;
            		std::cout << "           Specify a valid (positive) numax." << std::endl;
            		std::cout << "    The program will exit now. " << std::endl;
            		exit(EXIT_SUCCESS);
            	}
//            	}	
            }
            if(all_in.model_fullname == "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v3"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            	do_width_Appourchaux=1;
            	if(numax != -9999){ 
            		if (numax <= 0){ // Check that if there is a value, this one is a valid input for numax
            			std::cout << "    The model: " << all_in.model_fullname << " is supposed to have numax for (optional) argument" << std::endl;
            			std::cout << "    However, you provided a negative input. Please either:" << std::endl;
            			std::cout << "           [1] Not specify any numax or set !n -9999 in the .model file. In that case, the code will calculate a numax using weighted average of frequencies using Heights as weights" << std::endl;
            			std::cout << "           [2] Specify a valid (positive) numax." << std::endl;
            			std::cout << "    The program will exit now. " << std::endl;
            			exit(EXIT_SUCCESS);
            		}
            	}	
            }
            if(all_in.model_fullname == "model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            	do_width_Appourchaux=0;
            	if(numax != -9999){ 
            		if (numax <= 0){ // Check that if there is a value, this one is a valid input for numax
            			std::cout << "    The model: " << all_in.model_fullname << " is supposed to have numax for (optional) argument" << std::endl;
            			std::cout << "    However, you provided a negative input. Please either:" << std::endl;
            			std::cout << "           [1] Not specify any numax or set !n -9999 in the .model file. In that case, the code will calculate a numax using weighted average of frequencies using Heights as weights" << std::endl;
            			std::cout << "           [2] Specify a valid (positive) numax." << std::endl;
            			std::cout << "    The program will exit now. " << std::endl;
            			exit(EXIT_SUCCESS);
            		}
            	}	
            }
        }
        if(all_in.model_fullname == "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4"){
           	do_a11_eq_a12=1;
           	do_avg_a1n=1;
           	do_width_Appourchaux=1;
           	if (numax <= 0){ // Check that if there is a value, this one is a valid input for numax
           		std::cout << "    The model: " << all_in.model_fullname << " is supposed to have numax for argument" << std::endl;
           		std::cout << "    However, you provided a negative input. Please:" << std::endl;
           		std::cout << "           Specify a valid (positive) numax." << std::endl;
           		std::cout << "    The program will exit now. " << std::endl;
           		exit(EXIT_SUCCESS);
          	}
        }
        if(inputs_MS_global.common_names[i] == "fit_squareAmplitude_instead_Height" ){ 
        		do_amp=1;
        	if(inputs_MS_global.common_names_priors[i] != "bool"){
				fatalerror_msg_io_MS_Global("fit_squareAmplitude_instead_Height", "bool", "[0/1]  ", "1 " );
			} else{
				do_amp=inputs_MS_global.modes_common(i,0);
				std::cout << "Using do_amp = " << do_amp << std::endl;
			}

        }
    }
    if(all_in.model_fullname == " "){
    	std::cout << "Model name empty. Cannot proceed. Check that the .model file contains the model_fullname variable." << std::endl;
    	exit(EXIT_FAILURE);
    }
 	
 	io_calls.initialise_param(&Vis_in, lmax, Nmax_prior_params, -1, -1);
	io_calls.initialise_param(&Inc_in, 1, Nmax_prior_params, -1, -1);
		
	// -----------------------------------------------------------------
	// ------------ Handling Frequencies/Widths/Heights ----------------
	// -----------------------------------------------------------------
	Nf_el.setZero();
	int excluded_el=1; // The asymptotic model does not take into account any l=1 modes that may be provided into the configuration file. l=1 modes are handled by the global parameters
	for(int el=0; el<lmax+1; el++){
		//if (el != excluded_el){		
			els_eigen.resize(inputs_MS_global.eigen_params.col(0).size());
			for(int i=0; i<els_eigen.size();i++){
					els_eigen[i]=inputs_MS_global.eigen_params(i,0);
			}
			// --- sublist of relax taken from the first frequency list at the given el----
			pos_relax0=where_int(inputs_MS_global.els, el); // Get all the positions for the modes of degree el in the relax tab
			for(int i=0; i<pos_relax0.size(); i++){ // fill temporary variables for the given el 
				f_el.push_back(inputs_MS_global.freqs_ref[pos_relax0[i]]);
				rf_el.push_back(inputs_MS_global.relax_freq[pos_relax0[i]]);
				rw_el.push_back(inputs_MS_global.relax_gamma[pos_relax0[i]]);
				rh_el.push_back(inputs_MS_global.relax_H[pos_relax0[i]]);
			}
			
			pos_el=where_int(els_eigen, el); // find where we have modes of degree el in the eigen vectors
			Nf_el[el]=pos_el.size(); // Get the Number of frequencies of a given el... this will be used in plength
			for(int en=0; en<pos_el.size(); en++){
				f_inputs.push_back(inputs_MS_global.eigen_params(pos_el[en],1));
				f_priors_min.push_back(inputs_MS_global.eigen_params(pos_el[en],2));
				f_priors_max.push_back(inputs_MS_global.eigen_params(pos_el[en],3));
				if(el ==0){
					w_inputs.push_back(inputs_MS_global.eigen_params(pos_el[en],4)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE WIDTH ARE INTERPOLATED USING l=0
					h_inputs.push_back(inputs_MS_global.eigen_params(pos_el[en],5)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE HEIGHT ARE SCALED USING VISBILITIES				
				}
	
				pos_relax=where_dbl(f_el, inputs_MS_global.eigen_params(pos_el[en],1), tol); // Get the (unique) position in the relax tab for the seeked frequency
	
				if(pos_relax[0] != -1 && pos_relax.size() == 1){ // We impose uniqueness of the solution, within the tolerance
					f_relax.push_back(rf_el[pos_relax[0]]);
					if(el ==0){
						w_relax.push_back(rw_el[pos_relax[0]]);
						h_relax.push_back(rh_el[pos_relax[0]]);
					}
				} else{
					std::cout << "Error when preparing the vector of parameters using the MCMC file" << std::endl;
					std::cout << "The uniqueness of the frequency " << inputs_MS_global.eigen_params(pos_el[en],1) << " is not respected" << std::endl;
					std::cout << "     - pos_relax.size() =" << pos_relax.size() << std::endl;
					std::cout << "     - pos_relax[0] =" << pos_relax[0] << std::endl;
					if(pos_relax.size() != 1){
						std::cout << "Check that you do not have duplicated values in your MCMC file " << std::endl;
	
					}
					if(pos_relax[0] == -1){
						std::cout << "Check that you do not the same values in the relax list and in the eigen list of your MCMC file " << std::endl;
					}
					std::cout << "The program will exit now" << std::endl;
					exit(EXIT_FAILURE);
				}
			}
			// Empty the temporary variables for the relax and the freqs_ref at a given el
			f_el.resize(0);
			rf_el.resize(0);
			rw_el.resize(0);
			rh_el.resize(0);
		//}
	}	
	// ------------------------------------------------------------------------------------------
	// ------------------------------- Handling the Common parameters ---------------------------
	// ------------------------------------------------------------------------------------------
	if( do_amp){
		std::cout << "   ===> Requested to fit squared amplitudes instead of Height... Converting height inputs into A_squared = pi*Height*Width..." << std::endl;
		tmpstr_h="Amplitude_l0_rgb";
		for(int i=0; i<h_inputs.size(); i++){
			h_inputs[i]=pi*w_inputs[i]*h_inputs[i]; 
		}
	} else{
		tmpstr_h="Height_l0_rgb";
	}
    // Set default value of priors for Height Width and frequency
	io_calls.initialise_param(&height_in, h_relax.size(), Nmax_prior_params, -1, -1);
	if (do_width_Appourchaux == 0){
		io_calls.initialise_param(&width_in, w_relax.size(), Nmax_prior_params, -1, -1); 
	} 
	if (do_width_Appourchaux == 1){ // CONSTANT WIDTH CASE
		io_calls.initialise_param(&width_in, 1, Nmax_prior_params, -1, -1);
	}
	if (do_width_Appourchaux == 2){
		io_calls.initialise_param(&width_in, 6, Nmax_prior_params, -1, -1);
	}
	if (do_width_Appourchaux > 2 ){
		std::cout << "Fatal error on io_asymptotic.cpp: Allowed Appourchaux Widths models are only model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3 or model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3. This because other models" << std::endl;
		std::cout << "are not yet implemented. However the logical switches do_width_Appourchaux did not match and is not consistent with non-Appourchaux width"<< std::endl;
		std::cout << "Serious debug required here. The program will exit now" << std::endl;
		exit(EXIT_FAILURE);

	}
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3"
			|| "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v3"){
		// Generate initial set of l=1		
		// The initial values at fl0 + Dnu/2 + delta0l = inputs_MS_global.Dnu/100
		VectorXd fl0(Nf_el[0]); 
		for (int i=0; i<Nf_el[0] ; i++){
			fl0[i]=f_inputs[i];
		}
		fl1p_all=asympt_nu_p_from_l0_Xd(fl0, inputs_MS_global.Dnu, 1, inputs_MS_global.Dnu/100, fl0.minCoeff() - inputs_MS_global.Dnu, fl0.maxCoeff() + inputs_MS_global.Dnu);
		Nmixedmodes_params=Nmixedmodes_g_params + fl1p_all.size(); // Adding the fl1p modes in the parameter list 
		status_model=true;
	}
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4" || all_in.model_fullname =="model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4"){
		ferr.resize(inputs_MS_global.hyper_priors.rows()); 
		fref.resize(inputs_MS_global.hyper_priors.rows());		
		// Generate initial set of ferr and fref. fref is basically a grid of nodes for the x-axis of the bias splines		
		// We use the hyper-prior section of the .model file in order to get the list of fref. It is up to the user to provide that after examination of the spectrum
//		if (inputs_MS_global.hyper_priors.cols() > 3){
//			for (int i=0; i<inputs_MS_global.hyper_priors.rows()-1;i++){ // Check that we have values that are monotonically increasing
//				if (inputs_MS_global.hyper_priors(i+1,0) < inputs_MS_global.hyper_priors(i,0)){
//					std::cout << "Error in io_asymptotic, model " << all_in.model_fullname << std::endl;
//					std::cout << " You must ensure that values in the hyper prior section are non-zero and are monotonically increasing" << std::endl;
//				} 
//			}
//			for(int i=0; i<inputs_MS_global.hyper_priors.rows(); i++){
//				fref[i]=inputs_MS_global.hyper_priors(i,0);
//			}
//		} else{
		// MOD ON 24 OCT 2022: Check that it works for both scenario of extra params
			for (int i=0; i<inputs_MS_global.hyper_priors.rows()-1;i++){ // Check that we have values that are monotonically increasing
				if (inputs_MS_global.hyper_priors(i+1,0) < inputs_MS_global.hyper_priors(i,0)){
					std::cout << "Error in io_asymptotic, model " << all_in.model_fullname << std::endl;
					std::cout << " You must ensure that values in the hyper prior section are non-zero and are monotonically increasing" << std::endl;
				} 
			}
			for(int i=0; i<inputs_MS_global.hyper_priors.rows(); i++){
				fref[i]=inputs_MS_global.hyper_priors(i,0);
			}

			if (inputs_MS_global.hyper_priors.cols() !=1 ){ // If there is not all of the prior information, we expect a single column
				std::cout << "Error in io_asymptotic, model " << all_in.model_fullname << std::endl;
				std::cout << " You must have at least 4 non-zero data points in the hyper_priors section in order to describe the bias frequencies used as hanchors" << std::endl;
				std::cout << " Or you just specify the anchor points and default prior will be set" << std::endl;
				exit(EXIT_FAILURE);
			}
//		}
		if (inputs_MS_global.hyper_priors.cols() ==1){
			ferr.setZero(); // The initial bias vector is set to 0
		} else{
			for(int i=0; i<inputs_MS_global.hyper_priors.rows(); i++){
				ferr[i]=inputs_MS_global.hyper_priors(i,1);
			}
		}
		Nmixedmodes_params=Nmixedmodes_g_params + 2*fref.size() + 1; // Adding the fref and ferr modes in the parameter list + Hfactor : BEWARE: WHEN WRITTING THIS, I REALIZE THAT the +1 MAY BE MISSING IN v3 models... 
		status_model=true;
	} 
	if (status_model == false) {// status_model checks whether we already got a model defined. If not, it set the parameters to the default value
		Nmixedmodes_params=Nmixedmodes_g_params;
	}	
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike"){
		// delta0l=params[Nmax + lmax + Nfl0];
		//DPl=std::abs(params[Nmax + lmax + Nfl0 + 1])
		//alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2])
		// q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
		//sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]); // Evaluation of the inaccuracy in terms of Heights
		//sigma_m_0_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]); // Evaluation of the inaccuracy in terms of frequencies: Ordinate at origin ==> sigma_p
		//sigma_m_1_l1=std::abs(params[Nmax + lmax + Nfl0 + 6]); // Evaluation of the inaccuracy in terms of frequencies: slope ==> zeta=1 ==> sigma_g
    	//sigma_a1_0_l1=std::abs(params[Nmax + lmax + Nfl0 + 7]); // Evaluation of the inaccuracy in terms of splitting a1: Ordinate at origin ==> sigma_p
   		//sigma_a1_1_l1=std::abs(params[Nmax + lmax + Nfl0 + 8]); // Evaluation of the inaccuracy in terms of splitting a1: slope
    	//nmax=std::abs(params[Nmax + lmax + Nfl0 + 9]); // l=1 p modes: nmax ~ numax/Dnu + epsilon
    	//alpha_p=std::abs(params[Nmax + lmax + Nfl0 + 10])
		Nmixedmodes_params=11; 
	}

	int f_size=Nf_el[0] + Nmixedmodes_params + Nf_el[2] + Nf_el[3];

	//io_calls.initialise_param(&freq_in, f_relax.size(), Nmax_prior_params, -1, -1); 
	io_calls.initialise_param(&freq_in, f_size, Nmax_prior_params, -1, -1); 
	// DEFAULT HEIGHTS
	tmpXd.resize(4);
	tmpXd << Hmin, Hmax, -9999., -9999.; // default hmin and hmax for the Jeffreys prior
	for(int i=0; i<h_inputs.size(); i++){
		if(h_relax[i]){
			io_calls.fill_param(&height_in, tmpstr_h, "Jeffreys", h_inputs[i], tmpXd, i, 0);	
		} else{
			io_calls.fill_param(&height_in, tmpstr_h, "Fix", h_inputs[i], tmpXd, i, 0);			
		}
	}
		
	// --- Default setup for frequencies ---
	cpt=0;
	for(int i=0; i<f_inputs.size(); i++){
		if ( (i<Nf_el[0]) || (i>=Nf_el[0] + Nf_el[1]) ){ // We put frequencies of the model file only if those are not l=1
			if(f_relax[i]) {
				tmpXd << f_priors_min[i], f_priors_max[i], 0.0025*inputs_MS_global.Dnu, 0.0025*inputs_MS_global.Dnu; // default parameters for a GUG prior on frequencies
				io_calls.fill_param(&freq_in, "Frequency_RGB_l", "GUG", f_inputs[i], tmpXd, cpt, 0);	
			} else{
				io_calls.fill_param(&freq_in, "Frequency_RGB_l", "Fix", f_inputs[i], tmpXd, cpt, 0);			
			}
			cpt=cpt+1;
		} else{
			if(i == Nf_el[0]){
				cpt=cpt+Nmixedmodes_params;
			}
		}
	}

	status_model=false;
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v3"){
		// Add the fl1p modes into the frequency parameters 
		cpt=Nf_el[0] + Nmixedmodes_g_params;
		for(int i=0; i<fl1p_all.size();i++){
			tmpXd << fl1p_all[i] -10.*inputs_MS_global.Dnu/100 , fl1p_all[i] + 10.*inputs_MS_global.Dnu/100, 0.005*inputs_MS_global.Dnu, 0.005*inputs_MS_global.Dnu; // default parameters for a GUG prior on frequencies
			io_calls.fill_param(&freq_in, "Frequency_RGB_l1p", "GUG",  fl1p_all[i], tmpXd, cpt, 0);
			cpt=cpt+1;
		}
		status_model=true;
	}

	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4" || all_in.model_fullname =="model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4"){
		// Add the fref into the frequency parameters 
		cpt=Nf_el[0] + Nmixedmodes_g_params + 1; // The +1 is due to Hfactor here...
		for(int i=0; i<fref.size();i++){
			tmpXd << -9999, -9999, -9999, -9999; // The reference frequencies of the spline are not free parameters (they could, but my guess is that it will get messy...)
			io_calls.fill_param(&freq_in, "fref_bias", "Fix",  fref[i], tmpXd, cpt, 0);
			cpt=cpt+1;
		}
		// Use the default prior configuration if the user did specify only the fref frequencies (1 column)
		if (inputs_MS_global.hyper_priors.cols() == 1){
			// Add the ferr into the frequency parameters
			// The FIRST frequency (lower edge) is allowed to have larger excursion to the low frequency range
			std::cout << " inputs_MS_global.Dnu =" <<  inputs_MS_global.Dnu << std::endl;
			std::cout << " cpt =" << cpt << std::endl;
			std::cout << " ferr[0] =" << ferr[0] << std::endl;
			tmpXd << -inputs_MS_global.Dnu/2 , inputs_MS_global.Dnu/20, -9999, -9999; // default parameters for a Uniform prior
			io_calls.fill_param(&freq_in, "ferr_bias", "Uniform",  ferr[0], tmpXd, cpt, 0);
			cpt=cpt+1;
			// Core
			tmpXd << -inputs_MS_global.Dnu/20 , inputs_MS_global.Dnu/20, -9999, -9999; // default parameters for a Uniform prior
			for(int i=1; i<ferr.size()-1;i++){
				//tmpXd << fref[i-1] , fref[i+1], -9999, -9999; // default parameters for a Uniform prior
				io_calls.fill_param(&freq_in, "ferr_bias", "Uniform",  ferr[i], tmpXd, cpt, 0);
				cpt=cpt+1;
			}
			std::cout << "[4.3]" << std::endl;

			// The LAST frequency (upper edge) is allowed to have larger excursion to the high frequency range
			//tmpXd << fref[ferr.size()-2] , ferr[ferr.size()-1] + inputs_MS_global.Dnu/2, -9999, -9999; // default parameters for a Uniform prior
			tmpXd << -inputs_MS_global.Dnu/20 , inputs_MS_global.Dnu/2, -9999, -9999; // default parameters for a Uniform prior
			io_calls.fill_param(&freq_in, "ferr_bias", "Uniform",  ferr[ferr.size()-1], tmpXd, cpt, 0);
			cpt=cpt+1;
			status_model=true;
		} else{
			//tmpXd << -inputs_MS_global.Dnu/20 , inputs_MS_global.Dnu/20, -9999, -9999; // default parameters for a Uniform prior
			tmpXd.resize(4);	tmpXd.setConstant(-9999);
			for(int i=0; i<ferr.size();i++){
				for(int k=0; k<inputs_MS_global.hyper_priors.cols()-2;k++){
 					tmpXd[k]=inputs_MS_global.hyper_priors(i,2+k);
				}
				io_calls.fill_param(&freq_in, "ferr_bias", inputs_MS_global.hyper_priors_names[i], ferr[i], tmpXd, cpt, 0);	
				cpt=cpt+1;
			}
			status_model=true;
		}
	}
	if (status_model == false) {// status_model checks whether we already got a model defined. If not, it set the parameters to the default value
		std::cout << " Something went wrong in io_asymptotic: The model " << all_in.model_fullname << " did not pass an expected checkpoint " << std::endl;
		std::cout << " Please check the code and/or the model name" << std::endl;
		exit(EXIT_FAILURE);
	}	
	// ----------- Calculate numax -----------
	// Flated the output vector
	//tmpXd.resize(Nf_el.sum());
	//tmpXd.resize(height_in.inputs.size()); // We assume l=0,1,2,3 as p mode to get numax... this may not be true, but this is the best we can do at that stage of the code
	cpt=0;
	if(numax <=0){
		std::cout << "numax not provided. Input numax may be required by some models... Calculating numax ASSUMING PURE l=0 P MODES ONLY..." << std::endl;
		std::cout << "         ---- DUE TO MIXED MODES WE RECOMMEND YOU TO PROVIDE NUMAX AS AN ARGUMENT ---" << std::endl;
		if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3" ||
			all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v2" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike"    || 
			all_in.model_fullname == "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4" || all_in.model_fullname == "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v3"){
			std::cout << "     Test showed high risk of wrong estimates of the widths if numax is not provided because the number of modes in evolved stars can be small" << std::endl;
			std::cout << "     It is therefore now enforced that the user provide it in the case of models with AppWidth" << std::endl;
			std::cout << "     Please add numax in your model file following this syntax:" << std::endl;
			exit(EXIT_SUCCESS);
		}
		tmpXd=height_in.inputs;  // ALLWAY Nf_el[0] here because assuming p modes (the number of total modes listed does not matter)
		numax=getnumax(freq_in.inputs.segment(0, Nf_el[0]) , height_in.inputs); // We had to flatten the Height vector and put visibilities
		std::cout << "     numax: " << numax << std::endl;
	} else {	
		std::cout << " Using provided numax: " << numax << std::endl;
		if (err_numax <= 0){
			err_numax=0.05*numax;
			std::cout << " err_numax not provided. It will be set to 5% of numax: " << std::endl;
		} else{
			std::cout << " Using provided err_numax: " << err_numax << std::endl;
		}
	}
	std::cout << " ------------------" << std::endl;

	// 
	// -------------------------------------

	// ----- Switch between the models that handle averaging over n,l or both -----
    if(do_a11_eq_a12 == 1 && do_avg_a1n == 1){
        io_calls.initialise_param(&Snlm_in, 5, Nmax_prior_params, -1, -1);
    }  
    //if(do_a11_eq_a12 == 0 && do_avg_a1n == 1){
    //    io_calls.initialise_param(&Snlm_in, 7, Nmax_prior_params, -1, -1);
    //}
    //if(do_a11_eq_a12 == 1 && do_avg_a1n == 0){
	//	if (Nf_el[1] == Nf_el[2]){
	//		io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1], Nmax_prior_params, -1, -1);
	//	} else {
	//		std::cout << "When considering a11=a22"<< std::endl;
	//		std::cout <<" You must have as many l=1 than l=2" << std::endl;
	//		std::cout <<" Here we have: Nf(l=1) = " << Nf_el[1] << " and Nf(l=2) = " << Nf_el[2] << std::endl;
	//		std::cout <<" Check the .model file" << std::endl;
	//		std::cout <<" The program will exit now" << std::endl;
	//		exit(EXIT_FAILURE);
	//	}
    //}
    //if(do_a11_eq_a12 == 0 && do_avg_a1n == 0){
    //	io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2], Nmax_prior_params, -1, -1);
    //}
 
	// -------------- Set Extra_priors ----------------	
	extra_priors.resize(5); 
	extra_priors[0]=1; // By default, we apply a smoothness condition
	extra_priors[1]=2.; // By default, the smoothness coeficient is 2 microHz
	extra_priors[2]=0.2; // By default a3/a1<=0.2
	extra_priors[3]=0; // Switch to control whether a prior imposes Sum(Hnlm)_{m=-l, m=+l}=1. Default: 0 (none). >0 values are model_dependent
	extra_priors[4]=-1; // By default the model switch is not defined
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4"){ // The spline is handling inside the model, so there is no new priors for it. Prior MUST differ from v3
		extra_priors[4]=3;
	}
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3"){
		extra_priors[4]=1;
	}
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v3"){ // Deal with Dnu priors the same way as "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3"
		extra_priors[4]=1;
	}
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3"){ // Deal with Dnu priors the same way as "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3"
		extra_priors[4]=1;
	}
	if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike"){
		extra_priors[4]=2;
	}

	// ------------------------------------------------

	for(int i=0; i<inputs_MS_global.common_names.size(); i++){
		//std::cout << "i =" << i << std::endl;
		// --- Common parameters than can be run during setup ---
		if(inputs_MS_global.common_names[i] == "freq_smoothness" || inputs_MS_global.common_names[i] == "Freq_smoothness"){
			if(inputs_MS_global.common_names_priors[i] != "bool"){
				fatalerror_msg_io_MS_Global("freq_smoothness", "bool", "[0/1]     [Smoothness coeficient in microHz]", "1     2.0" );
			} else{
				extra_priors[0]=inputs_MS_global.modes_common(i,0);
				extra_priors[1]=inputs_MS_global.modes_common(i,1);
			}
		}	
		if(inputs_MS_global.common_names[i] == "trunc_c"){
			if(inputs_MS_global.common_names_priors[i] != "Fix"){
				fatalerror_msg_io_MS_Global("trunc_c", "Fix", "[Truncation parameter]", "20" );
			} else{
				trunc_c=inputs_MS_global.modes_common(i,0); // This value is added to the input vector at the end of this function.
			}
		}	
		if(inputs_MS_global.common_names[i] == "model_type"){
			if(inputs_MS_global.common_names_priors[i] != "Fix"){
				fatalerror_msg_io_MS_Global("model_type", "Fix", "[Value]", "(0: use of fl0 to construct fl1p, 1: use a 2nd order polynomial to fit fl0 to construct O2p fl1p)" );
			} else{
				model_type=inputs_MS_global.modes_common(i,0); // This value is added to the input vector at the end of this function.
			}
		}	
		if(inputs_MS_global.common_names[i] == "bias_type"){
			if(inputs_MS_global.common_names_priors[i] != "Fix"){
				fatalerror_msg_io_MS_Global("bias_type", "Fix", "[Value]", "(0: window bias, 1: Cubic spline, 2: Hermite spline" );
			} else{
				bias_type=inputs_MS_global.modes_common(i,0); // This value is added to the input vector at the end of this function.
			}
		}	
		// --- Frequencies ---
		if(inputs_MS_global.common_names[i] == "Frequency" || inputs_MS_global.common_names[i] == "frequency"){ 
			if(inputs_MS_global.common_names_priors[i] == "GUG" || inputs_MS_global.common_names_priors[i] == "Uniform"){
				int p_effective=0;
				for(int p0=0; p0<f_inputs.size(); p0++){
					if ( (p0 < Nf_el[0]) || (p0 >= (Nf_el[0] + Nf_el[1])) ){ // Here WE SKIP THE L=1 modes TO FILL f_input for ONLY l=0, 2, 3 into freq_in. Filling l=1 uses specific keywords for global parameters
						if(f_relax[p0]){ 
							if(inputs_MS_global.common_names_priors[i] == "GUG"){
								tmpXd << f_priors_min[p0], f_priors_max[p0], inputs_MS_global.modes_common(i,3), inputs_MS_global.modes_common(i,4);
							} else{
								tmpXd << f_priors_min[p0], f_priors_max[p0], -9999, -9999; 						
							}
							io_calls.fill_param(&freq_in, "Frequency_l", inputs_MS_global.common_names_priors[i], f_inputs[p0], tmpXd, p_effective, 0);	
						} else{
							io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f_inputs[p0], tmpXd, p_effective, 0);			
						}
						p_effective=p_effective+1;
					}
					else{ // WHENEVER WE HAVE TO DEAL WITH l=1 modes, this is where the global parameters must be
						if(p0 == Nf_el[0]){
							p_effective=p_effective + Nmixedmodes_params; // Let space for the l=1 mixed modes global parameters
							// Warning here: Nf_el[1] must be updated to be equal to Nmixedmodes_params just after the loop scanning all the global parameters
						}
					}
				}
			} else{
				fatalerror_msg_io_MS_Global(inputs_MS_global.common_names[i], "GUG or Uniform", "", "" );
			}
		}
		// ---- l=1 mixed modes Global parameters ---
		if(inputs_MS_global.common_names[i] == "delta01" && all_in.model_fullname != "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3" && all_in.model_fullname != "model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3"){ // Valid parameter only for a specific model
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				p0=Nf_el[0];
//				fatalerror_msg_io_MS_Global("delta01", "Fix_Auto", "", "" );
                tmpstr="Uniform";
                tmpXd.resize(4);
                tmpXd << -1.*inputs_MS_global.Dnu/100, 1.*inputs_MS_global.Dnu/100, -9999., -9999.;
                std::cout << "Fix_Auto requested for delta01..., the prior will be with this syntax (negative delta01):" << std::endl;
                std::cout << "          " << std::left << std::setw(15) << tmpstr << " [-10*Deltanu/100]  [+10*Deltanu/100]   0.00000      -9999    -9999" << std::endl;
                std::cout << "          Initial guess: Dnu/100  = " << 0.5*inputs_MS_global.Dnu/100 << std::endl;
                io_calls.fill_param(&freq_in, "delta01", tmpstr, 0.5*inputs_MS_global.Dnu/100 ,  tmpXd, p0, 0);
			} else{
				p0=Nf_el[0];
				io_calls.fill_param(&freq_in, "delta01", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
		} 
		if(inputs_MS_global.common_names[i] == "DP1"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("DP1", "Fix_Auto", "", "" );
			}
			p0=Nf_el[0]+1;
			io_calls.fill_param(&freq_in, "DP1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "alpha_g"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("alpha_g", "Fix_Auto", "", "" );
			}
			p0=Nf_el[0]+2;
			io_calls.fill_param(&freq_in, "alpha_g", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "q"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("q", "Fix_Auto", "", "" );
			}
			p0=Nf_el[0]+3;
			io_calls.fill_param(&freq_in, "q", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "sigma_Hl1"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("sigma_Hl1", "Fix_Auto", "", "" );
			}
			p0=Nf_el[0]+4;
			io_calls.fill_param(&freq_in, "sigma_Hl1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "sigma_fl1g"){
			//std::cout << "sigma_fl1g" << std::endl;
			if (extra_priors[4] == 1){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
					fatalerror_msg_io_MS_Global("sigma_fl1g", "Fix_Auto", "", "" );
				}
				p0=Nf_el[0]+5;
				io_calls.fill_param(&freq_in, "sigma_fl1g", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
			//std::cout << "		DONE" << std::endl;
		}
		if(inputs_MS_global.common_names[i] == "sigma_fl1m" && all_in.model_fullname != "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4" && all_in.model_fullname != "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4"){
			//std::cout << "sigma_fl1m" << std::endl;
			if (extra_priors[4] == 1){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value			
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
					fatalerror_msg_io_MS_Global("sigma_fl1m", "Fix_Auto", "", "" );
				}
				p0=Nf_el[0]+6;
				io_calls.fill_param(&freq_in, "sigma_fl1m", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
		}
		if(inputs_MS_global.common_names[i] == "Wfactor"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("Wfactor", "Fix_Auto", "", "" );
			}
			p0=Nf_el[0]+6;
			io_calls.fill_param(&freq_in, "Wfactor", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "Hfactor"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("Hfactor", "Fix_Auto", "", "" );
			}
			p0=Nf_el[0]+7;
			io_calls.fill_param(&freq_in, "Hfactor", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "sigma_fl1m_0"){
			//std::cout << "sigma_fl1m_0" << std::endl;
			if (extra_priors[4] == 2){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto" || inputs_MS_global.common_names_priors[i] == "Auto"){ // ONLY the manual mode is handled at the moment
					fatalerror_msg_io_MS_Global("sigma_fl1m_0", "Fix_Auto", "", "" );
				}
				// Those models have two priors for fl1m sigma inaccuracy terms (linear function of zeta(nu), see model in models.cpp)
				p0=Nf_el[0]+5;
				io_calls.fill_param(&freq_in, "sigma_fl1m_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
			//std::cout << "		DONE" << std::endl;
		}
		if(inputs_MS_global.common_names[i] == "sigma_fl1m_1"){
			//std::cout << "sigma_fl1m_1" << std::endl;
			if (extra_priors[4] == 2){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto" || inputs_MS_global.common_names_priors[i] == "Auto"){ // ONLY the manual mode is handled at the moment
					fatalerror_msg_io_MS_Global("sigma_fl1m_1", "Fix_Auto", "", "" );
				}
				// Those models have two priors for fl1m sigma inaccuracy terms (linear function of zeta(nu), see model in models.cpp)
				p0=Nf_el[0]+6;
				io_calls.fill_param(&freq_in, "sigma_fl1m_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
			//std::cout << "		DONE" << std::endl;
		}
		if(inputs_MS_global.common_names[i] == "sigma_a1_0"){
			//std::cout << "sigma_a1_0" << std::endl;
			if (extra_priors[4] == 2){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto" || inputs_MS_global.common_names_priors[i] == "Auto"){ // ONLY the manual mode is handled at the moment
					fatalerror_msg_io_MS_Global("sigma_a1_0", "Fix_Auto", "", "" );
				}
				// Those models have two priors for fl1m sigma inaccuracy terms (linear function of zeta(nu), see model in models.cpp)
				p0=Nf_el[0]+7;
				io_calls.fill_param(&freq_in, "sigma_a1_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
			//std::cout << "		DONE" << std::endl;
		}
		if(inputs_MS_global.common_names[i] == "sigma_a1_1"){
			std::cout << "sigma_a1_1" << std::endl;
			if (extra_priors[4] == 2){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto" || inputs_MS_global.common_names_priors[i] == "Auto"){ // ONLY the manual mode is handled at the moment
					fatalerror_msg_io_MS_Global("sigma_a1_1", "Fix_Auto", "", "" );
				}
				// Those models have two priors for fl1m sigma inaccuracy terms (linear function of zeta(nu), see model in models.cpp)
				p0=Nf_el[0]+8;
				io_calls.fill_param(&freq_in, "sigma_a1_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			}
			//std::cout << "		DONE" << std::endl;
		}
		if(inputs_MS_global.common_names[i] == "nmax"){
			std::cout << "nmax" << std::endl;
			if (extra_priors[4] == 2){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value
				if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ // ONLY the AUTO mode is handled at the moment. The user only choose the initial guess
					fatalerror_msg_io_MS_Global("nmax", "Fix_Auto", "", "" );
				}
				// We impose a Gaussian, centered on numax/Dnu + epsilon with stddev of 0.5
				p0=Nf_el[0]+9;
				tmpXd.resize(4);
				tmpXd << numax/inputs_MS_global.Dnu + inputs_MS_global.C_l/inputs_MS_global.Dnu, 0.5, -9999., -9999.;
				std::cout << " Applying a Gaussian prior N(mu, sig) with :"<< std::endl;
				std::cout << "      mu = numax/inputs_MS_global.Dnu + inputs_MS_global.C_l/inputs_MS_global.Dnu = " << numax/inputs_MS_global.Dnu + inputs_MS_global.C_l/inputs_MS_global.Dnu << std::endl;
				std::cout << "      sig= " << tmpXd[1] << std::endl;
				io_calls.fill_param(&freq_in, "nmax", "Gaussian", numax/inputs_MS_global.Dnu + inputs_MS_global.C_l/inputs_MS_global.Dnu , tmpXd, p0, 0);	
			}
			//std::cout << "		DONE" << std::endl;
		}
		if(inputs_MS_global.common_names[i] == "alpha_p"){
			std::cout << "alpha_p" << std::endl;
			if (extra_priors[4] == 2){ // Deal with this parameters only for specific models as per defined by extra_priors[4] value
				if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ // ONLY the Auto mode is handled at the moment. The user only choose the initial guess
					fatalerror_msg_io_MS_Global("alpha_p", "Fix_Auto", "", "" );
				}
				// We impose a Jeffreys's prior with minimum threshold 0.01 and maxmimum threshold 0.2
				p0=Nf_el[0]+10;
				tmpXd.resize(4);
				tmpXd << 0.01, 0.2, -9999., -9999.;
				std::cout << " Applying a Jeffreys prior J(xmin, xmax) with :"<< std::endl;
				std::cout << "      xmin = " << tmpXd[0] << std::endl;
				std::cout << "      xmax= " << tmpXd[1] << std::endl;
				io_calls.fill_param(&freq_in, "alpha_p", "Jeffreys", 0.0001, tmpXd, p0, 0);	
			}
			std::cout << "		DONE" << std::endl;
		}
		// --- Height or Amplitude ---
		if(inputs_MS_global.common_names[i] == "height" || inputs_MS_global.common_names[i] == "Height" || 
		   inputs_MS_global.common_names[i] == "amplitude" || inputs_MS_global.common_names[i] == "Amplitude"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
					fatalerror_msg_io_MS_Global(inputs_MS_global.common_names[i], "Fix_Auto", "", "" );
				}
			for(p0=0; p0<h_inputs.size(); p0++){
				if(h_relax[p0]){
					io_calls.fill_param(&height_in, tmpstr_h,  inputs_MS_global.common_names_priors[i], h_inputs[p0],  inputs_MS_global.modes_common.row(i), p0, 0);	
				} else{
					io_calls.fill_param(&height_in, tmpstr_h,  "Fix", h_inputs[p0],  inputs_MS_global.modes_common.row(i), p0, 1);		
				}
			}
		}
		// -- Mode Width ---
		if((inputs_MS_global.common_names[i] == "width" || inputs_MS_global.common_names[i] == "Width") && (do_width_Appourchaux == 0)){
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
                    tmpstr="Jeffreys";
                    tmpXd.resize(4);
                    tmpXd << resol, inputs_MS_global.Dnu/3., -9999., -9999.;
                    std::cout << "Fix_Auto requested for Widths... For all free Widths, the prior will be with this syntax:" << std::endl;
                    std::cout << "          " << std::left << std::setw(15) << tmpstr << " [Spectrum Resolution]   [Deltanu / 3]   -9999    -9999" << std::endl;
                    std::cout << "          " << "Resolution: " << resol << std::endl;
                } else{
                    tmpstr=inputs_MS_global.common_names_priors[i];
                    tmpXd=inputs_MS_global.modes_common.row(i);
                }
			for(p0=0; p0<w_inputs.size(); p0++){
				if(w_relax[p0]){
					io_calls.fill_param(&width_in, "Width_l", tmpstr, w_inputs[p0],  tmpXd, p0, 0);
				} else{
					io_calls.fill_param(&width_in, "Width_l",  "Fix", w_inputs[p0],  inputs_MS_global.modes_common.row(i), p0, 1);		
				}
			}
		} 	
		if((inputs_MS_global.common_names[i] == "width" || inputs_MS_global.common_names[i] == "Width") && (do_width_Appourchaux == 1)){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
                tmpstr="Jeffreys";
                tmpXd.resize(4);
                tmpXd << resol, inputs_MS_global.Dnu/3., -9999., -9999.;
                std::cout << "Fix_Auto requested for Widths... For all free Widths, the prior will be with this syntax:" << std::endl;
                std::cout << "          " << std::left << std::setw(15) << tmpstr << " [Spectrum Resolution]   [Deltanu / 3]   -9999    -9999" << std::endl;
                std::cout << "          " << "Resolution: " << resol << std::endl;
                tmp=0;
                for(p0=0; p0<w_inputs.size(); p0++){ // Compute the mean
                	tmp=tmp + w_inputs[p0]/w_inputs.size();
                }
                std::cout << " Mean Width:" << tmp << std::endl;
                io_calls.fill_param(&width_in, "Width_l", tmpstr, tmp,  tmpXd, 0, 0);
            } else{
                std::cout << "Please use the 'Fix_Auto' option for the width when considering a constant width for the whole range of fit (do_width_Appourchaux == 1)" << std::endl;
                exit(EXIT_SUCCESS);
            }
		}	
		// --- Splittings and asymetry ---
		if(inputs_MS_global.common_names[i] == "rot_env" || inputs_MS_global.common_names[i] == "Rot_env"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("rot_env", "Fix_Auto", "", "" );
			}
			p0=0;
			io_calls.fill_param(&Snlm_in, "rot_env", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "rot_core" || inputs_MS_global.common_names[i] == "Rot_core"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("rot_core", "Fix_Auto", "", "" );
			}
			p0=1;
			io_calls.fill_param(&Snlm_in, "rot_core", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		if(inputs_MS_global.common_names[i] == "asphericity_eta"|| inputs_MS_global.common_names[i] == "Asphericity_eta"){  
			std::cout << " Warning: Update on the code (07/12/2021) does not allow to use Asphericity_eta as a keyword in io_asymptotic models" << std::endl;
			std::cout << "          The calculation of eta is now made systematically within the models and no longer requires the eta parameter" << std::endl;
			std::cout << "          Consequently, this parameter will be ignored" << std::endl;
			Snlm_in.inputs_names[1]="Empty";
			Snlm_in.priors_names[1]="Fix";
			Snlm_in.inputs_names[1]="Asphericity_eta";
			Snlm_in.relax[1]=0;
			Snlm_in.inputs[1]=0;
		}

		if(inputs_MS_global.common_names[i] == "splitting_a3" || inputs_MS_global.common_names[i] == "Splitting_a3"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("splitting_a3", "Fix_Auto", "", "" );
			}
			p0=3;
			io_calls.fill_param(&Snlm_in, "Splitting_a3", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
		}

		if(inputs_MS_global.common_names[i] == "asymetry" || inputs_MS_global.common_names[i] == "Asymetry"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("asymetry", "Fix_Auto", "", "" );
			}
			p0=4;
			io_calls.fill_param(&Snlm_in, "Lorentzian_asymetry", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
		}
		// --- Dealing with visibilities ---
		if(inputs_MS_global.common_names[i] == "visibility_l1" || inputs_MS_global.common_names[i] == "Visibility_l1"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("visibility_l1", "Fix_Auto", "", "" );
			}
			if(lmax >= 1){
				p0=0;
				io_calls.fill_param(&Vis_in, "Visibility_l1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			} else{
				std::cout << "Warning: lmax=" << lmax << " but keyword 'visibility_l1' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << std::endl;
			}
		}
		if(inputs_MS_global.common_names[i] == "visibility_l2" || inputs_MS_global.common_names[i] == "Visibility_l2"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("visibility_l2", "Fix_Auto", "", "" );
			}
			if(lmax >= 2){
				p0=1;
				io_calls.fill_param(&Vis_in, "Visibility_l2", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			} else{
				std::cout << "Warning: lmax=" << lmax << " but keyword 'visibility_l2' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << std::endl;
			}
		}
		if(inputs_MS_global.common_names[i] == "visibility_l3" || inputs_MS_global.common_names[i] == "Visibility_l3"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("visibility_l3", "Fix_Auto", "", "" );
			}
			if(lmax >= 3){
				p0=2;
				io_calls.fill_param(&Vis_in, "Visibility_l3", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			} else{
				std::cout << "Warning: lmax=" << lmax << " but keyword 'visibility_l3' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << std::endl;
			}
		}
		// --- Finally the inclination ---
		if(inputs_MS_global.common_names[i] == "inclination" || inputs_MS_global.common_names[i] == "Inclination"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("inclination", "Fix_Auto", "", "" );
			}
			if( inputs_MS_global.modes_common(i,0) >= 90){ // Avoid some nan when computing the prior
				    tmp=89.99999; 
			} else{
				tmp=inputs_MS_global.modes_common(i,0);
			}
			p0=0;
			io_calls.fill_param(&Inc_in, "Inclination", inputs_MS_global.common_names_priors[i], tmp, inputs_MS_global.modes_common.row(i), p0, 1);
		}
		
		if(inputs_MS_global.common_names[i] == "sqrt(splitting_a1).cosi"){
			std::cout << " WARNING : keyword sqrt(splitting_a1).cosi found but is incompatible with models in io_asymptotic.cpp" << std::endl;
			std::cout << "           The keyword will be ignored. Please be sure to have the keywords 'rot_core' and rot_env' to describe the rotation" << std::endl;
			//if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
			//fatalerror_msg_io_MS_Global("sqrt(splitting_a1).cosi", "Fix_Auto", "", "" );
			//}
			/*p0=3;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);

            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0  and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1cosi=1;
            */
		}
		
		if(inputs_MS_global.common_names[i] == "sqrt(splitting_a1).sini"){
			std::cout << " WARNING : keyword sqrt(splitting_a1).sini found but is incompatible with models in io_asymptotic.cpp" << std::endl;
			std::cout << "           The keyword will be ignored. Please be sure to have the keywords 'rot_core' and rot_env' to describe the rotation" << std::endl;
			//if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
			//	fatalerror_msg_io_MS_Global("sqrt(splitting_a1).sini", "Fix_Auto", "", "" );
			//}
			/*p0=4;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			
            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0 and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1sini=1;
            */
		}
	}

// ----- Required changes to ensure that sizes of l=1 arrays are fine ----
Nf_el[1]=Nmixedmodes_params;
// ---------------------------

if(bool_a1cosi != bool_a1sini){ // Case when one of the projected splitting quantities is missing ==> Problem
	std::cout << "Warning: Both 'sqrt(splitting_a1).sini' and 'sqrt(splitting_a1).cosi' keywords must appear" << std::endl;
	std::cout << "         It is forbidden to use only one of them" << std::endl;
	std::cout << "         Edit the .MODEL file accordingly" << std::endl;
	std::cout << "The program will exit now" << std::endl;
	exit(EXIT_FAILURE);
}
/*
if((bool_a1cosi == 0) && (bool_a1sini == 0)){ // Case where Inclination and Splitting_a1 are supposed to be used (CLASSIC and CLASSIC_vX models)
	if( (all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike") || (all_in.model_fullname == "model_MS_Global_a1etaa3_Harvey1985") ||
		(all_in.model_fullname == "model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1") || (all_in.model_fullname=="model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2")){
		std::cout << "Warning: Splitting_a1 and Inclination keywords detected while requested model is model_MS_Global_a1etaa3_HarveyLike... ADAPTING THE VARIABLE FOR ALLOWING THE FIT TO WORK" << std::endl;
		
		std::cout << "         Replacing variables splitting_a1 and inclination by sqrt(splitting_a1).cos(i) and sqrt(splitting_a1).sin(i)..." << std::endl;

		if( Inc_in.priors_names[0]=="Fix" && Snlm_in.priors_names[0]=="Fix"){
			std::cout << "         - Inclination and splitting_a1 were 'Fixed' ==> Fixing sqrt(splitting_a1).cos(i) and sqrt(splitting_a1).sin(i)" << std::endl; 
			
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", "Fix", sqrt(Snlm_in.inputs[0])*cos(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 3, 0);
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", "Fix", sqrt(Snlm_in.inputs[0])*sin(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 4, 0);
		} else{
			std::cout << "	       - Inclination and/or Splitting_a1 requested as a free parameter ==> Freeing both sqrt(splitting_a1).cos(i) and sqrt(splitting_a1).sin(i)" << std::endl;
			std::cout << "         - Setting the prior to the one requested for the Splitting_a1 assuming non-informative prior on inclination: " << Snlm_in.priors_names[0] << std::endl;
			std::cout << "         - Value of the prior are set for cosi=sini=1, ie the max is set by sqrt(splitting_a1)" << std::endl;
			std::cout << "           Beware that this may not be what you want!" << std::endl;

			Snlm_in.priors(1,0)=std::sqrt(Snlm_in.priors(1,0)); // Convert the max boundary of splitting_a1 into the max boundary of sqrt(splitting_a1)
			
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", Snlm_in.priors_names[0], sqrt(Snlm_in.inputs[0])*cos(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 3, 0);
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", Snlm_in.priors_names[0], sqrt(Snlm_in.inputs[0])*sin(Inc_in.inputs[0]*pi/180.), Snlm_in.priors.col(0), 4, 0);
		}

		if(Snlm_in.inputs[3] < 1e-2){ // Avoid issues with the edge of uniform priors
			Snlm_in.inputs[3]=1e-2;
		}
		if(Snlm_in.inputs[4] < 1e-2){ // Avoid issues with the edge of uniform priors
			Snlm_in.inputs[4]=1e-2;
		}

        std::cout << "            Slots for the variables Inclination and Splitting are forced to be FIXED" << std::endl;
        std::cout << "            Be aware that may lead to unwished results if you are not careful about the used model" << std::endl;
        io_calls.fill_param(&Inc_in, "Empty", "Fix", 0, Inc_in.priors.col(0), 0, 1); // Note that inputs_MS_global.modes_common.row(0) is not used... just dummy
        io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, Snlm_in.priors.col(0), 0, 1); // Note that inputs_MS_global.modes_common.row(0) is not used... just dummy
	}

}
if((bool_a1cosi == 1) && (bool_a1sini ==1)){
	if (all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic"){
		std::cout << "Warning: We cannot use " << all_in.model_fullname << " with variables sqrt(splitting_a1).cosi and sqrt(splitting_a1).sini..." << std::endl;
		std::cout << "         Use Inclination and Splitting_a1 instead for the model this model"  << std::endl;
		std::cout << "	       Alternatively, you can use 'model_MS_Global_a1etaa3_HarveyLike with sqrt(splitting_a1).cosi and sqrt(splitting_a1).sini keywords'" <<std::endl;
		std::cout << "         The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
    std::cout << "    NOTICE: sqrt(splitting_a1).sini and sqrt(splitting_a1).cosi superseed inclination and splitting" << std::endl;
    std::cout << "            Slots for the variables Inclination and Splitting are therefore forced to be FIXED" << std::endl;
    std::cout << "            Be aware that may lead to unwished results if you are not careful about the used model" << std::endl;
   	
   	io_calls.fill_param(&Inc_in, "Empty", "Fix", 0, Inc_in.priors.col(0), 0, 1); // Note that inputs_MS_global.modes_common.row(0) is not used... just dummy   	
   	io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, Snlm_in.priors.col(0), 0, 1); // "Splitting_a1" default values are erased
}
*/
	// ----------------------------------------------------
	// ---------------- Handling noise --------------------
	// ----------------------------------------------------
	io_calls.initialise_param(&Noise_in, 10, Nmax_prior_params, -1, -1);
	set_noise_params(&Noise_in, inputs_MS_global.noise_s2, inputs_MS_global.noise_params); 	// Set defaults
	// -------------------------------------------------------------------------------------
	// ------------- Sticking everything together in a Input_Data structure --------------
	// -------------------------------------------------------------------------------------
	
	//if(model_type ==-1 && bias_type==-1 && (all_in.model_fullname != "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4")){ // If it is a model that does not include any model_type (or bias_type) value, then it has no bias function 
	plength.resize(11);
	plength[0]=h_inputs.size(); plength[1]=lmax                  ; plength[2]=Nf_el[0];
	plength[3]=Nf_el[1]       ; plength[4]=Nf_el[2]		   ; plength[5]=Nf_el[3];
	plength[6]=Snlm_in.inputs.size(); plength[7]=width_in.inputs.size() ; plength[8]=Noise_in.inputs.size(); 
	plength[9]=Inc_in.inputs.size();
 	if (model_type ==-1 && bias_type ==-1){
		plength[10]=3; // This is trunc_c, do_amp and sigma_limit;
	} else{
		plength[10]=6; // This is trunc_c, do_amp and sigma_limit, model_type, bias_type, Nferr;
	}
	io_calls.initialise_param(&all_in, plength.sum(), Nmax_prior_params, plength, extra_priors); // The final two arguments are used for model_type or bias_type
	if ((model_type ==-1 && bias_type !=-1) || (model_type !=-1 && bias_type ==-1)){
		std::cout << "Error in io_asymptotic: model_type and bias_type must be defined for the specified model:" << all_in.model_fullname << std::endl;
		exit(EXIT_FAILURE);
	}
	// --- Put the Height or Amplitudes---
	p0=0;
	io_calls.add_param(&all_in, &height_in, p0); // height_in may contain either height or amplitudes depending on specified keywords
	// --- Put the Visibilities ---
	p0=all_in.plength[0];
	io_calls.add_param(&all_in, &Vis_in, p0);
	
	// --- Put the Frequencies ---
	p0=all_in.plength[0] + all_in.plength[1];	
	io_calls.add_param(&all_in, &freq_in, p0);

	// --- Put the Snlm (splittings and asymetry) ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5];
	io_calls.add_param(&all_in, &Snlm_in, p0);
	
	// --- Put the Width ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] +  all_in.plength[3]  + all_in.plength[4] + all_in.plength[5] + all_in.plength[6];
	if(do_width_Appourchaux == 0 || do_width_Appourchaux == 1){ // Most models
		io_calls.add_param(&all_in, &width_in, p0);
	} 
	if(do_width_Appourchaux == 2){// Case of: model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1
		width_in=set_width_App2016_params_v2(numax, err_numax, width_in);
		io_calls.add_param(&all_in, &width_in, p0);
	}
	// --- Put the Noise ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] +  all_in.plength[3]  + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7];
	io_calls.add_param(&all_in, &Noise_in, p0);
	// --- Put the Inclination ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8];
	io_calls.add_param(&all_in, &Inc_in, p0);
	
	// --- Add trunc_c that controls the truncation of the Lorentzian ---
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9];
	io_calls.fill_param(&all_in, "Truncation parameter", "Fix", trunc_c, inputs_MS_global.modes_common.row(0), p0, 1);
	if (all_in.inputs[p0] <= 0){
		std::cout << "Warning: trunc_c <= 0. This is forbidden. Setting default to 10000. (No truncation)" << std::endl;
		all_in.inputs[p0]=10000.; // In case of a non-sense value for c, we use Full-Lorentzian as default
	}
	// -- Add the Amplitude switch --
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9] + 1;
	io_calls.fill_param(&all_in, "Switch for fit of Amplitudes or Heights", "Fix", do_amp, inputs_MS_global.modes_common.row(0), p0,1);

	// -- Add the Limit on sigma... set to 10% of Dnu in the begining of this code --
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9] + 2;
	tmpXd.resize(4);
	tmpXd << -9999, -9999, -9999, -9999;
	io_calls.fill_param(&all_in, "Maximum limit on random values generated by N(0,sigma_m)", "Fix", sigma_limit, tmpXd, p0,1);
	// -- Add the hyper-priors for the bias function (if any)--
	if (model_type !=-1){
		p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9] + 3;
		tmpXd.resize(4);
		tmpXd << -9999, -9999, -9999, -9999;
		io_calls.fill_param(&all_in, "model type ", "Fix", model_type, tmpXd, p0,1);
	}
	if (bias_type !=-1){
		p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9] + 4;
		tmpXd.resize(4);
		tmpXd << -9999, -9999, -9999, -9999;
		io_calls.fill_param(&all_in, "bias type ", "Fix", bias_type, tmpXd, p0,1);
		io_calls.fill_param(&all_in, "Nferr ", "Fix", ferr.size(), tmpXd, p0+1,1); // Nferr
	}

	//
	if(verbose == 1){
		std::cout << " ----------------- Configuration summary -------------------" << std::endl;
		std::cout << "Model Name = " << all_in.model_fullname << std::endl;
		if(all_in.extra_priors[0] == 0){ 
			std::cout << "   freq_smoothness is set to 0 ==> NO smoothness condition on frequencies" << std::endl;
		} else{
			std::cout << "   freq_smoothness is set to 1 ==> APPLIES a smoothness condition on frequencies" << std::endl;
			std::cout << "   smoothness coeficient as specified by the user (of defined by default): " << all_in.extra_priors[1] << " microHz" << std::endl;
		}
		
		std::cout << "    Maximum ratio between a3 and a1: " << all_in.extra_priors[2] << std::endl;
			
		std::cout << " -----------------------------------------------------------" << std::endl;
		std::cout << " ---------- Configuration of the input vectors -------------" << std::endl;
		std::cout << "    Table of inputs " << std::endl;
		io_calls.show_param(all_in, 1);
		std::cout << " -----------------------------------------------------------" << std::endl;

		std::cout << "    The fit includes:" << std::endl;
		std::cout << "          - " << all_in.plength[0] << "  Heights parameters "<< std::endl; 
		std::cout << "          - " << all_in.plength[1] << "  Visibility parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[2] << "  l=0 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[3] << "  l=1 Frequencies parameters (Global parameters)" << std::endl; 
		std::cout << "          - " << all_in.plength[4] << "  l=2 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[5] << "  l=3 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[6] << "  Rotation parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[7] << "  Width parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[8] << "  Noise parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[9] << "  Inclination parameters" << std::endl; 
		std::cout << std::endl << "          - " << "Total Number of Parameters: " << all_in.plength.sum() << std::endl; 
		std::cout << " -----------------------------------------------------------" << std::endl;
	}
	
//	std::cout << "Exiting test " << std::endl;
//	exit(EXIT_SUCCESS);

return all_in;
}
