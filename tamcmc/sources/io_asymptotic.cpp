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
#include "colors.hpp"

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
	VectorXd aj_param_count(8);
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
	aj_param_count.setZero();
	all_in.model_fullname=" "; // Default is an empty string
	bool passed_model=false;
    for(int i=0; i<inputs_MS_global.common_names.size(); i++){
        if(inputs_MS_global.common_names[i] == "model_fullname" ){ // This defines if we assume S11=S22 or not (the executed model will be different)
        	all_in.model_fullname=inputs_MS_global.common_names_priors[i];
			if (all_in.model_fullname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4" || 
									all_in.model_fullname == "model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4"){
				std::cout << colors::red << " Models model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4 and model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4" << std::endl;
				std::cout << " have now a new name: model_RGB_asympt_aj_AppWidth_HarveyLike_v4 and model_RGB_asympt_aj_CteWidth_HarveyLike_v4" << std::endl;
				std::cout << " Update your model files accordingly." << colors::white  << std::endl;
				exit(EXIT_FAILURE);
			}
            if(all_in.model_fullname == "model_RGB_asympt_aj_AppWidth_HarveyLike_v4"){
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
				passed_model=true;	
            }
			if(all_in.model_fullname == "model_RGB_asympt_aj_CteWidth_HarveyLike_v4"){
				do_a11_eq_a12=1;
				do_avg_a1n=1;
				do_width_Appourchaux=1;
				if(numax != -9999){ 
					if (numax <= 0){ // Check that if there is a value, this one is a valid input for numax
						std::cout << "    The model: " << all_in.model_fullname << " is supposed to have numax for argument" << std::endl;
						std::cout << "    However, you provided a negative input. Please:" << std::endl;
						std::cout << "           Specify a valid (positive) numax." << std::endl;
						std::cout << "    The program will exit now. " << std::endl;
						exit(EXIT_SUCCESS);
					}
				}
				passed_model=true;	
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
    if(passed_model == false){
    	std::cout << colors::red << "Model Not recognized. Many models were removed in version >1.86.0. Check if your model is still supported." << std::endl;
    	std::cout << "    Supported models: " << std::endl;
		std::cout << "           - model_RGB_asympt_aj_CteWidth_HarveyLike_v4" << std::endl;
		std::cout << "           - model_RGB_asympt_aj_AppWidth_HarveyLike_v4" << colors::white << std::endl;
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
	int Nfix=0;
	if (all_in.model_fullname == "model_RGB_asympt_aj_AppWidth_HarveyLike_v4" || all_in.model_fullname =="model_RGB_asympt_aj_CteWidth_HarveyLike_v4"){
		ferr.resize(inputs_MS_global.hyper_priors.rows()); 
		fref.resize(inputs_MS_global.hyper_priors.rows());		
		// Generate initial set of ferr and fref. fref is basically a grid of nodes for the x-axis of the bias splines		
		// We use the hyper-prior section of the .model file in order to get the list of fref. It is up to the user to provide that after examination of the spectrum
		// MOD ON 24 OCT 2022 + Mod on 22 Sept 2023 (adding the possibility to fix the bias)
			for (int i=0; i<inputs_MS_global.hyper_priors.rows()-1;i++){ 
				// Check that we have values that are monotonically increasing
				if (inputs_MS_global.hyper_priors(i+1,0) < inputs_MS_global.hyper_priors(i,0)){
					std::cout << "Error in io_asymptotic, model (A) " << all_in.model_fullname << std::endl;
					std::cout << " You must ensure that values in the hyper prior section are all FIXED OR ALL non-zero and are monotonically increasing" << std::endl;
					exit(EXIT_FAILURE);
				} 
				// Counter to know how many Fix hyper param. If all in Fix ==> no bias correction. If one in Fix and other not ==> NOT ALLOWED
				if (inputs_MS_global.hyper_priors_names[i] == "Fix"){
					Nfix=Nfix+1;
				} 				
			}
			// Check if fixed params are either all or none
			if (Nfix != inputs_MS_global.hyper_priors_names.size()-1 && Nfix != 0){
				std::cout << "Error in io_asymptotic, model (B) " << all_in.model_fullname << std::endl;
				std::cout << " You must ensure that values in the hyper prior section are all FIXED OR ALL non-zero and are monotonically increasing" << std::endl;
				std::cout << " Nfix = " << Nfix << std::endl;
				std::cout << " inputs_MS_global.hyper_priors_names.size() = " << inputs_MS_global.hyper_priors_names.size() << std::endl;
				exit(EXIT_FAILURE);
			}

			for(int i=0; i<inputs_MS_global.hyper_priors.rows(); i++){
				fref[i]=inputs_MS_global.hyper_priors(i,0);
			}

			if (inputs_MS_global.hyper_priors.cols() <1 ){ // If there is not all of the prior information, we expect a single column
				std::cout << " inputs_MS_global.hyper_priors.rows() = " << inputs_MS_global.hyper_priors.cols() << std::endl;
				std::cout << " inputs_MS_global.hyper_priors = " << inputs_MS_global.hyper_priors << std::endl;
				std::cout <<  "Error in io_asymptotic, model " << all_in.model_fullname << std::endl;
				std::cout << " You must have at least 4 non-zero data points per line in the hyper_priors section in order to describe the bias frequencies used as hanchors" << std::endl;
				std::cout << " Or you just specify the anchor points (1 column entries) and default prior will be set" << std::endl;
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
	if (all_in.model_fullname == "model_RGB_asympt_aj_AppWidth_HarveyLike_v4" || all_in.model_fullname =="model_RGB_asympt_aj_CteWidth_HarveyLike_v4"){
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
			tmpXd << -inputs_MS_global.Dnu/2 , inputs_MS_global.Dnu/20, -9999, -9999; // default parameters for a Uniform prior
			io_calls.fill_param(&freq_in, "ferr_bias", "Uniform",  ferr[0], tmpXd, cpt, 0);
			cpt=cpt+1;
			// Core
			tmpXd << -inputs_MS_global.Dnu/20 , inputs_MS_global.Dnu/20, -9999, -9999; // default parameters for a Uniform prior
			for(int i=1; i<ferr.size()-1;i++){
				io_calls.fill_param(&freq_in, "ferr_bias", "Uniform",  ferr[i], tmpXd, cpt, 0);
				cpt=cpt+1;
			}

			// The LAST frequency (upper edge) is allowed to have larger excursion to the high frequency range
			tmpXd << -inputs_MS_global.Dnu/20 , inputs_MS_global.Dnu/2, -9999, -9999; // default parameters for a Uniform prior
			io_calls.fill_param(&freq_in, "ferr_bias", "Uniform",  ferr[ferr.size()-1], tmpXd, cpt, 0);
			cpt=cpt+1;
			status_model=true;
		} else{
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
	cpt=0;
	if(numax <=0){
		std::cout << "numax not provided. Input numax may be required by some models... Calculating numax ASSUMING PURE l=0 P MODES ONLY..." << std::endl;
		std::cout << "         ---- DUE TO MIXED MODES WE RECOMMEND YOU TO PROVIDE NUMAX AS AN ARGUMENT ---" << std::endl;
		if (all_in.model_fullname == "model_RGB_asympt_aj_AppWidth_HarveyLike_v4" ){
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
        io_calls.initialise_param(&Snlm_in, 10, Nmax_prior_params, -1, -1); // Add (a1_env,a1_core,a2_env,a2_core, a3_env, a4_env,a5_env,a6_env) +  1 val (bool 0/1) for eta0 + 1 asymetry
    }  
 
	// -------------- Set Extra_priors ----------------	
	extra_priors.resize(5); 
	extra_priors[0]=1; // By default, we apply a smoothness condition
	extra_priors[1]=2.; // By default, the smoothness coeficient is 2 microHz
	extra_priors[2]=0.2; // By default a3/a1<=0.2
	extra_priors[3]=0; // Switch to control whether a prior imposes Sum(Hnlm)_{m=-l, m=+l}=1. Default: 0 (none). >0 values are model_dependent
	extra_priors[4]=-1; // By default the model switch is not defined
	if (all_in.model_fullname == "model_RGB_asympt_aj_AppWidth_HarveyLike_v4"){ // The spline is handling inside the model, so there is no new priors for it. Prior MUST differ from v3
		extra_priors[4]=3;
	}
	if (all_in.model_fullname == "model_RGB_asympt_aj_CteWidth_HarveyLike_v4"){ // The spline is handling inside the model, so there is no new priors for it. Prior MUST differ from v3
		extra_priors[4]=3;
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
				if(Nfix != inputs_MS_global.hyper_priors_names.size()-1){
					bias_type=inputs_MS_global.modes_common(i,0); // This value is added to the input vector at the end of this function.
				} else{
					bias_type=0; // Use this value to signify that no bias will be used in the code (ALL FIXED)

				}
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
		/* 2nd Order Polynomial for p modes: curvature
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
		*/
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
		aj_param_count=aj_param_count + settings_aj_splittings_RGB(i, inputs_MS_global, &Snlm_in);

		if(inputs_MS_global.common_names[i] == "asymetry" || inputs_MS_global.common_names[i] == "Asymetry"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("asymetry", "Fix_Auto", "", "" );
			}
			p0=9;
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
				std::cout << colors::yellow << "Warning: lmax=" << lmax << " but keyword 'visibility_l1' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << colors::white << std::endl;
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
				std::cout << colors::yellow << "Warning: lmax=" << lmax << " but keyword 'visibility_l2' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << colors::white << std::endl;
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
				std::cout << colors::yellow << "Warning: lmax=" << lmax << " but keyword 'visibility_l3' detected" << std::endl;
				std::cout << "         This visibilitiy input will be ignored" << std::endl;
				std::cout << "         Proceeding..." << colors::white << std::endl;
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
			std::cout << colors::yellow << " WARNING : keyword sqrt(splitting_a1).cosi found but is incompatible with models in io_asymptotic.cpp" << std::endl;
			std::cout << "           The keyword will be ignored. Please be sure to have the keywords 'rot_core' and rot_env' to describe the rotation" << colors::white << std::endl;
		}
		
		if(inputs_MS_global.common_names[i] == "sqrt(splitting_a1).sini"){
			std::cout << colors::yellow << " WARNING : keyword sqrt(splitting_a1).sini found but is incompatible with models in io_asymptotic.cpp" << std::endl;
			std::cout << "           The keyword will be ignored. Please be sure to have the keywords 'rot_core' and rot_env' to describe the rotation" << colors::white << std::endl;
		}
	}

	if (aj_param_count.sum() !=8){
		std::cout << colors::red << "aj_param_count = " << aj_param_count.transpose() << std::endl;
		std::cout << " Invalid number of constraints for the model." << std::endl;
		std::cout << " Please set a1_env (or rot_env), a1_core (rot_core), a2_core, a2_env, a3_env, a4_env, a5_env, a6_env (8 parameters) for the considered model" << colors::white << std::endl;
		exit(EXIT_SUCCESS); 	
	} else{
		aj_param_count.setZero();
	}
	// ----- Required changes to ensure that sizes of l=1 arrays are fine ----
	Nf_el[1]=Nmixedmodes_params;
	// ---------------------------

	// ----------------------------------------------------
	// ---------------- Handling noise --------------------
	// ----------------------------------------------------
	io_calls.initialise_param(&Noise_in, 10, Nmax_prior_params, -1, -1);
	set_noise_params(&Noise_in, inputs_MS_global.noise_s2, inputs_MS_global.noise_params); 	// Set defaults
	// -------------------------------------------------------------------------------------
	// ------------- Sticking everything together in a Input_Data structure --------------
	// -------------------------------------------------------------------------------------
	
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
	if (bias_type == 0 && Nfix != inputs_MS_global.hyper_priors_names.size()-1 ){
		std::cout << "  Warning : Inconcistency detected between the requested bias_type = 0" << std::endl;
		std::cout << "            and the hyper priors being not set to FIX. In this scenario, " << std::endl;
		std::cout << "            bias_type superseed the hyper priors ==> Setting hyper priors to FIX and 0" << std::endl;
		tmpXd.resize(4);
		tmpXd << -9999, -9999, -9999, -9999;
		std::vector<int> positions=io_calls.lookupIndices("ferr_bias", all_in.inputs_names);
		for (p0=0; p0<positions.size();p0++){ // Fix to 0 all ferr_bias.
			io_calls.fill_param(&all_in, "ferr_bias", "Fix", 0, tmpXd, positions[p0],1);
		}
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
	
	//std::cout << "Exiting test " << std::endl;
	//exit(EXIT_SUCCESS);

return all_in;
}


// This function handles all of the aj related options. 
// It is designed for aj models of the RGB phase
VectorXd settings_aj_splittings_RGB(const int i, const MCMC_files inputs_MS_global, Input_Data* Snlm_in){
	int p0;
	VectorXd aj_param_count(8); // counter for knowing where we passed
	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure

	aj_param_count.setZero();
	if(inputs_MS_global.common_names[i] == "rot_env" || inputs_MS_global.common_names[i] == "Rot_env" || inputs_MS_global.common_names[i] == "a1_env"){  // This is a valid keyword only for aj models
		aj_param_count[0]=aj_param_count[0]+1;
		p0=0; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "rot_env", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Fix_Auto requested for rot_env ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "rot_core"  || inputs_MS_global.common_names[i] == "Rot_core" || inputs_MS_global.common_names[i] == "a1_core"){ // This is a valid keyword only for aj models
		aj_param_count[0]=aj_param_count[1]+1;
		p0=1; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "rot_core", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		}
		else{
			std::cout << "    Fix_Auto requested for rot_core ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "a2_core"){  // This is a valid keyword only for aj models
		aj_param_count[1]=aj_param_count[2]+1;
		p0=2; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "a2_core", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Fix_Auto requested for a2_core ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "a2_env"){ // This is a valid keyword only for aj models
		aj_param_count[1]=aj_param_count[3]+1;
		p0=3; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "a2_env", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Fix_Auto requested for a2_env ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "a3_env"){  // This is a valid keyword only for aj models
		aj_param_count[2]=aj_param_count[4]+1;
		p0=4; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "a3_env", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Fix_Auto requested for a3_env ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "a4_env"){ // This is a valid keyword only for aj models
		aj_param_count[2]=aj_param_count[5]+1;
		p0=5; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "a4_env", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Fix_Auto requested for a4_env ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "a5_env"){  // This is a valid keyword only for aj models
		aj_param_count[3]=aj_param_count[6]+1;
		p0=6; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "a5_env", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Fix_Auto requested for a5_env ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "a6_env"){ // This is a valid keyword only for aj models
		aj_param_count[3]=aj_param_count[7]+1;
		p0=7; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
			io_calls.fill_param(Snlm_in, "a6_env", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Fix_Auto requested for a6_env ... This is not allowed. Please use an explicit prior " << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	if(inputs_MS_global.common_names[i] == "asphericity_eta"|| inputs_MS_global.common_names[i] == "Asphericity_eta"){  
		std::cout << colors::yellow << " Warning: Update on the code (07/12/2021) does not allow to use Asphericity_eta as a keyword in io_asymptotic models" << std::endl;
		std::cout << "          The calculation of eta0 is now made if eta0_switch = 1. Otherwise (eta0_switch = 0, it is set to 0 and might be included in a2 coefficients)" << std::endl;
		std::cout << "          Consequently, this parameter will be igno" << colors::white << std::endl;
		Snlm_in->inputs_names[8]="eta0_switch";
		Snlm_in->priors_names[8]="Fix";
		Snlm_in->relax[8]=0;
		Snlm_in->inputs[8]=0;
	}
	if(inputs_MS_global.common_names[i] == "eta0_switch"){ // This is a valid keyword only for aj models
		p0=8; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
		if(inputs_MS_global.common_names_priors[i] == "Fix"){ 
			io_calls.fill_param(Snlm_in, "eta0_switch", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
		} else{
			std::cout << "    Error: eta0_switch must be a boolean 0 or 1: eta0_switch=0 means no computation of the a2_CF. eta0_switch means a2_CF is computed and included." << std::endl;
			std::cout << "           Due to the fact that this is RGB models with a2_core and a2_env, we recommend the default eta0_switch=0" << std::endl;
			exit(EXIT_SUCCESS);
		}
	}
	return aj_param_count;
}
