/*
 * io_ms_global.cpp
 *
 * Contains all kind of functions
 * used to arrange/handle the strings for the MS_Global models
 * 
 *  Created on: 20 Jun 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "data.h" // contains the structure Data
//#include "string_handler.h"
#include "io_local.h"
#include "io_models.h"
#include "noise_models.h"
#include "function_rot.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


MCMC_files read_MCMC_file_local(const std::string cfg_model_file, const int slice_ind,  const bool verbose){

	int i, out, nl, el, cpt, range_counter;
	std::vector<int> pos;
	std::string line0, char0, char1, char2;
	std::vector<std::string> word, tmp;
	std::ifstream cfg_session;
    MatrixXd tmpXd;

	MCMC_files i_local;

	i_local.numax=-9999; // Initialize the optional variable numax
	
	range_counter=0;
    cpt=0;
    i=0;
    out=0;
    nl=200; // maximum number of lines
    el=4; // maximum degree of the modes
    cfg_session.open(cfg_model_file.c_str());
    if (cfg_session.is_open()) {

	char0="#"; 
	std::getline(cfg_session, line0);
	if(verbose == 1) {std::cout << "  - Global parameters:" << std::endl;}
	i_local.freq_range.resize(2);
	i_local.els.resize(400);
	i_local.freqs_ref.resize(400);
	while((out < 3) && (!cfg_session.eof())){

		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1));
		char1=line0.substr(1, 1);
		if (char0 == "#" && char1 == "K"){
			word=strsplit(line0, "= \t");
			i_local.ID=strtrim(strtrim(word[1]));
			if(verbose == 1) {std::cout << "           ID=" << i_local.ID << std::endl;}
		}
		if (char0 == "!" && char1 == "n"){
			word=strsplit(line0, " ");
			i_local.numax=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           numax =" << i_local.numax << std::endl;}		
		}
		if (char0 == "!" && char1 != "!" && char1 != "n"){
			word=strsplit(line0, " ");
			i_local.Dnu=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           Dnu =" << i_local.Dnu << std::endl;}
		}
		if (char0 == "!" && char1 == "!"){
			word=strsplit(line0, " ");
			i_local.C_l=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           C_l =" << i_local.C_l << std::endl;}
		}
		if (char0 == "*"){
			if (range_counter == slice_ind){
				word=strsplit(line0, " ");
				i_local.freq_range[0]=str_to_dbl(word[1]);
				i_local.freq_range[1]=str_to_dbl(word[2]);
				if(verbose == 1) {std::cout << "       Analysed freq_range = [" << i_local.freq_range[0] << " , " << i_local.freq_range[1] << "]" << std::endl;}
			} else{
				word=strsplit(line0, " ");
				if(verbose == 1) {std::cout << "       Skipped freq_range = [" << str_to_dbl(word[1]) << " , " << str_to_dbl(word[2]) << "]  (skipped as we do not analyse this slice now...)" << std::endl;}
			}
			range_counter=range_counter + 1;	
		}
		if ((char0 != "#") && (char0 != "!") && (char0 != "*")){
			word=strsplit(line0, " ");
			if(strtrim(word[0]) == "p" || strtrim(word[0]) == "g" || strtrim(word[0]) == "co"){
				i_local.param_type.push_back(strtrim(word[0]));
				i_local.els[cpt]=str_to_int(word[1]);
                i_local.freqs_ref[cpt]=str_to_dbl(word[2]);
                if (word.size() >=4){
                    i_local.relax_freq.push_back(str_to_bool(word[3]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    i_local.relax_freq.push_back(1);  // By default, we fit frequencies
                }
                if (word.size() >=5){
                    //i_local.relax_gamma.push_back(str_to_bool(word[4]));
                    i_local.relax_H.push_back(str_to_bool(word[4]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the Width of the mode at frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    //i_local.relax_gamma.push_back(1);
                    i_local.relax_H.push_back(1);
                }
                if (word.size() >=6){
                    //i_local.relax_H.push_back(str_to_bool(word[5]));
                    i_local.relax_gamma.push_back(str_to_bool(word[5]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the Height of the mode at frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    //i_local.relax_H.push_back(1);
                    i_local.relax_gamma.push_back(1);
                }
                cpt=cpt+1;
			} else{
				std::cout << "Problem with the keyword that identifies the type of mode" << std::endl;
				std::cout << "The type of mode should be either 'p', 'g' or 'co'" << std::endl;
				std::cout << "THe program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
			
		}
		if (char0 == "#"){
			out=out+1;
		}	
		i=i+1;
		char0="";
		char1="";
		std::getline(cfg_session, line0);
	}
	i_local.freqs_ref.conservativeResize(i_local.relax_freq.size());
	i_local.els.conservativeResize(i_local.relax_freq.size());
	if(verbose == 1) {
		std::cout << " - List of relaxed parameters:" << std::endl;
		std::cout << "          l    nu       relax(nu)   relax(W)  relax(H)" << std::endl;
		for(int k=0; k<i_local.freqs_ref.size();k++){
			std::cout << "              " << i_local.els[k] << "       " << i_local.freqs_ref[k] << "       " << i_local.relax_freq[k] << "  "
				  << i_local.relax_gamma[k] << "      " << i_local.relax_H[k] << std::endl;
		}
	}
    
	// -------------------------------------
// Modified on 18 Aug 2022 to handle an array of values for the hyper priors
// This was specifically redesigned to handle spline fitting with hyper parameters
// The assumed structure is: 
//   value for x-axis, prior type, prior_parameters
	i=0;
	cpt=0;
	if(verbose == 1) {std::cout << " - Hyper priors:" << std::endl;}
	while ((out < 4) && !cfg_session.eof()){ // the priors, until we reach the next # symbol
			std::getline(cfg_session, line0);
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				word=strsplit(line0, " \t");
				if (cpt == 0){ // We need to detect the number of columns of the hyper_priors
					if(word.size() == 1){ // Case where no prior information is actually. It is more a compatibility mode with old code.				
						i_local.hyper_priors.resize(50,1);
					}
					if(word.size() ==2){
						std::cout << "  Error in the definition of the hyper priors" << std::endl;
						std::cout << "  only two columns detected while it should be either one column" << std::endl;
						std::cout << "  or at least columns (for a Fix hyper prior)" << std::endl;
						std::cout << "  Check your model file" << std::endl;
						exit(EXIT_FAILURE);
					}
					if(word.size() >2){ // We have at least 3 elements: The hyper_prior value in the x-axis, prior_type and then a hyper_prior parameter
						i_local.hyper_priors.resize(50,word.size()-1);
					}
				}
				// Gather all numerical values in hyper_priors
				i_local.hyper_priors(i,0)=str_to_dbl(word[0]); 
				for(int j=2; j<word.size();j++){
					i_local.hyper_priors(i,j-1)=str_to_dbl(word[j]);
				}
				// Gather all the prior types on hyper_priors_type
				i_local.hyper_priors_names.push_back(word[1]);
				//i_local.hyper_priors.row(i)=str_to_Xdarr(line0, " \t");
				cpt=cpt+1;
			} else{
				 out=out+1;
			}
			i=i+1;
	  }
	  i_local.hyper_priors.conservativeResize(cpt, word.size());
	  if(verbose == 1) {
		std::cout << i_local.hyper_priors.transpose() << std::endl;
	  }


	i=0;
	cpt=0;
	i_local.eigen_params.resize(200,6);
	std::getline(cfg_session, line0);
	if(verbose == 1) {std::cout << " - Initial guesses and frequency priors:" << std::endl;}
	while ((out < 5) && !cfg_session.eof() ){ //the priors, until we reach the 5th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				word=strsplit(line0, " \t");
				i_local.eigen_params.row(i)=str_to_Xdarr(line0, " \t");
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}
	i_local.eigen_params.conservativeResize(cpt, 6);
	if(verbose == 1) {
		std::cout << i_local.eigen_params << std::endl;
	}
	// -------------------------------------

	i=0;
	cpt=0;
	i_local.noise_params.resize(10);
	while ( (out < 6) && !cfg_session.eof() ){  // the priors, until we reach the 6th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){		
				word=strsplit(line0," \t");
          			for(int j=0; j<word.size(); j++){
          			 	i_local.noise_params[cpt+j]=str_to_dbl(word[j]);
          			}
				cpt=cpt+word.size();
		 	} else{
				 out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	  }
    // Ensure a compatibility with the 3 Harvey profile + White noise standard
      if(cpt < 10){
        tmpXd=i_local.noise_params.segment(0,cpt);
        i_local.noise_params.setConstant(-1);
        i_local.noise_params.segment(10-cpt, cpt)=tmpXd;
      }
	if(verbose == 1) {
		std::cout << " - Noise inputs:" << std::endl;
		std::cout << i_local.noise_params.transpose() << std::endl;
	}
	// -------------------------------------	

	i=0;
	cpt=0;
	i_local.noise_s2.resize(10, 3);
	while ((out < 7) && !cfg_session.eof() ){ //the priors, until we reach the 7th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				word=strsplit(line0," \t");
				i_local.noise_s2.row(i)=str_to_Xdarr(line0, " \t");
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}
    if(i < 11){
        tmpXd=i_local.noise_s2.block(0,0,cpt,i_local.noise_s2.cols());
        i_local.noise_s2.setConstant(-1);
        i_local.noise_s2.block(10-cpt, 0, tmpXd.rows(), tmpXd.cols())=tmpXd;
    }
    if(verbose == 1) {
		std::cout << " - Noise information extracted during step s2:" << std::endl;
		std::cout << i_local.noise_s2 << std::endl;
    }
	// -------------------------------------	
	
	// -------- The common parameters follow ------
	VectorXd a;
	i=0;
	cpt=0;
	i_local.modes_common.resize(20,5);
	i_local.modes_common.setConstant(-9999); // up to 10 variables and 4 prior parameters
	while ( (out < 9) && !cfg_session.eof() ){ // the initial values for the common parameters + priors, until we reach the 9th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			word=strsplit(line0," \t");
			if (char0 != "#"){
				word=strsplit(line0," \t");
				i_local.common_names.push_back(strtrim(word[0]));
				i_local.common_names_priors.push_back(strtrim(word[1]));
 				word.erase(word.begin()); // erase the slot containing the keyword name
				word.erase(word.begin()); // erase the slot containing the type of prior/switch
				a=arrstr_to_Xdarrdbl(word);
				
				for(int k=0; k<a.size();k++){
					i_local.modes_common(cpt, k)=a[k];
				}
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}

	i_local.modes_common.conservativeResize(i_local.common_names.size(), 5);

	if(verbose == 1) {
		std::cout << " - Common parameters for modes:" << std::endl;
		for(i=0; i<i_local.common_names.size();i++){
			std::cout << "  " << i_local.common_names[i] << "  ";
			std::cout << "  " << i_local.common_names_priors[i] << "  ";
			std::cout << "  " << i_local.modes_common.row(i) << std::endl;
		}
	}
    
	// -------------------------------------
   cfg_session.close();
   } else {
   		std::cout << "Unable to open the file: " << cfg_model_file << std::endl;
   		std::cout << "Check that the file exist and that the path is correct" << std::endl;
   		std::cout << "Cannot proceed" << std::endl;
   		std::cout << "The program will exit now" << std::endl;
   		exit(EXIT_FAILURE);
   }
      
   return i_local;
}

Input_Data build_init_local(const MCMC_files inputs_local, const bool verbose, const double resol){

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

	double rho=pow(inputs_local.Dnu/Dnu_sun,2.) * rho_sun;
	double Dnl=0.75, trunc_c=-1;
	//double numax=inputs_local.numax;
	
	// All Default booleans
	bool do_a11_eq_a12=1, do_avg_a1n=1, do_amp=0;
	bool bool_a1sini=0, bool_a1cosi=0;
	int lmax, en, Ntot, p0, cpt, ind, pos_prior_height=-1;
	//uint8_t do_width_Appourchaux=0; // We need more than a boolean here, but no need to use a 64 bit signed int
	double tol=1e-2, tmp;
	VectorXi pos_el, pos_relax0, els_eigen, Nf_el(4), plength;
	VectorXd ratios_l, extra_priors, tmpXd;
	MatrixXd tmpXd_mat;
	std::vector<int> pos_relax, pos_OK;
	std::vector<double> f0_inputs, h0_inputs, w0_inputs, f0_priors_min, f0_priors_max, f_el;
	std::vector<double> f1_inputs, h1_inputs, w1_inputs, f1_priors_min, f1_priors_max;
	std::vector<double> f2_inputs, h2_inputs, w2_inputs, f2_priors_min, f2_priors_max;
	std::vector<double> f3_inputs, h3_inputs, w3_inputs, f3_priors_min, f3_priors_max;
	std::vector<bool> f0_relax, h0_relax, w0_relax; 
	std::vector<bool> f1_relax, h1_relax, w1_relax; 
	std::vector<bool> f2_relax, h2_relax, w2_relax; 
	std::vector<bool> f3_relax, h3_relax, w3_relax; 
	std::vector<int> rf_el, rw_el, rh_el;
	
	std::string tmpstr_h, tmpstr;
	
	Input_Data Snlm_in, Inc_in, Noise_in, freq_in, height_in, width_in;//, Vis_in; //, width_App2016_params; // This is by block, each category of parameters		
	Input_Data all_in; // The final structure of parameters, using the standards of my code
	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure

	// Flatening and ordering of all the inputs/relax variables
	lmax=inputs_local.els.maxCoeff();
	
	// -- Initialisation of structures --
    // --- Look for common instruction That must be run before the setup ---------
	all_in.model_fullname=" "; // Default is an empty string
    for(int i=0; i<inputs_local.common_names.size(); i++){
        if(inputs_local.common_names[i] == "model_fullname" ){ // This defines if we assume S11=S22 or not (the executed model will be different)
        	all_in.model_fullname=inputs_local.common_names_priors[i];
            if(all_in.model_fullname == "model_MS_local_basic"){
            	//Previously corresponding to average_a1nl     bool    1    1 
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            }
            if(all_in.model_fullname == "model_MS_local_Hnlm"){
            	//Previously corresponding to average_a1nl     bool    1    1 
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            }
        }
        if(inputs_local.common_names[i] == "fit_squareAmplitude_instead_Height" ){ 
        		do_amp=1;
        	if(inputs_local.common_names_priors[i] != "bool"){
				fatalerror_msg_io_local("fit_squareAmplitude_instead_Height", "bool", "[0/1]  ", "1 " );
			} else{
				do_amp=inputs_local.modes_common(i,0);
				std::cout << "Using do_amp = " << do_amp << std::endl;
			}
        }
    }
    if(all_in.model_fullname == " "){
    	std::cout << "Model name empty. Cannot proceed. Check that the .model file contains the model_fullname variable." << std::endl;
    	exit(EXIT_FAILURE);
    }
 	//io_calls.initialise_param(&Vis_in, lmax, Nmax_prior_params, -1, -1); // NOT USED FOR LOCAL FIT
	io_calls.initialise_param(&Inc_in, 1, Nmax_prior_params, -1, -1);	
	
	// -----------------------------------------------------------------
	// ------------ Handling Frequencies/Widths/Heights ----------------
	// -----------------------------------------------------------------
	Nf_el.setZero();
	for(int el=0; el<lmax+1; el++){
		els_eigen.resize(inputs_local.eigen_params.col(0).size());
		for(int i=0; i<els_eigen.size();i++){
			els_eigen[i]=inputs_local.eigen_params(i,0);
		}
		// --- sublist of relax taken from the first frequency list at the given el----
		pos_relax0=where_int(inputs_local.els, el); // Get all the positions for the modes of degree el in the relax tab
		if (pos_relax0[0] != -1){
			for(int i=0; i<pos_relax0.size(); i++){ // fill temporary variables for the given el 
				f_el.push_back(inputs_local.freqs_ref[pos_relax0[i]]);
				rf_el.push_back(inputs_local.relax_freq[pos_relax0[i]]);
				rw_el.push_back(inputs_local.relax_gamma[pos_relax0[i]]);
				rh_el.push_back(inputs_local.relax_H[pos_relax0[i]]);
			}
			pos_el=where_int(els_eigen, el); // find where we have modes of degree el in the eigen vectors
			Nf_el[el]=pos_el.size(); // Get the Number of frequencies of a given el... this will be used in plength
			for(int en=0; en<pos_el.size(); en++){
				if(el ==0){
					f0_inputs.push_back(inputs_local.eigen_params(pos_el[en],1));
					f0_priors_min.push_back(inputs_local.eigen_params(pos_el[en],2));
					f0_priors_max.push_back(inputs_local.eigen_params(pos_el[en],3));
					w0_inputs.push_back(inputs_local.eigen_params(pos_el[en],4)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE WIDTH ARE INTERPOLATED USING l=0
					h0_inputs.push_back(inputs_local.eigen_params(pos_el[en],5)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE HEIGHT ARE SCALED USING VISBILITIES				
				}
				if(el ==1){
					f1_inputs.push_back(inputs_local.eigen_params(pos_el[en],1));
					f1_priors_min.push_back(inputs_local.eigen_params(pos_el[en],2));
					f1_priors_max.push_back(inputs_local.eigen_params(pos_el[en],3));
					w1_inputs.push_back(inputs_local.eigen_params(pos_el[en],4)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE WIDTH ARE INTERPOLATED USING l=0
					h1_inputs.push_back(inputs_local.eigen_params(pos_el[en],5)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE HEIGHT ARE SCALED USING VISBILITIES				
				}
				if(el ==2){
					f2_inputs.push_back(inputs_local.eigen_params(pos_el[en],1));
					f2_priors_min.push_back(inputs_local.eigen_params(pos_el[en],2));
					f2_priors_max.push_back(inputs_local.eigen_params(pos_el[en],3));
					w2_inputs.push_back(inputs_local.eigen_params(pos_el[en],4)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE WIDTH ARE INTERPOLATED USING l=0
					h2_inputs.push_back(inputs_local.eigen_params(pos_el[en],5)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE HEIGHT ARE SCALED USING VISBILITIES				
				}
				if(el ==3){
					f3_inputs.push_back(inputs_local.eigen_params(pos_el[en],1));
					f3_priors_min.push_back(inputs_local.eigen_params(pos_el[en],2));
					f3_priors_max.push_back(inputs_local.eigen_params(pos_el[en],3));
					w3_inputs.push_back(inputs_local.eigen_params(pos_el[en],4)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE WIDTH ARE INTERPOLATED USING l=0
					h3_inputs.push_back(inputs_local.eigen_params(pos_el[en],5)); // COLUMNS for el>0 ARE IGNORED BECAUSE THE HEIGHT ARE SCALED USING VISBILITIES				
				}			
				pos_relax=where_dbl(f_el, inputs_local.eigen_params(pos_el[en],1), tol); // Get the (unique) position in the relax tab for the seeked frequency

				if(pos_relax[0] != -1 && pos_relax.size() == 1){ // We impose uniqueness of the solution, within the tolerance
					if(el ==0){
						f0_relax.push_back(rf_el[pos_relax[0]]);
						w0_relax.push_back(rw_el[pos_relax[0]]);
						h0_relax.push_back(rh_el[pos_relax[0]]);
					}
					if(el ==1){
						f1_relax.push_back(rf_el[pos_relax[0]]);
						w1_relax.push_back(rw_el[pos_relax[0]]);
						h1_relax.push_back(rh_el[pos_relax[0]]);
					}
					if(el ==2){
						f2_relax.push_back(rf_el[pos_relax[0]]);
						w2_relax.push_back(rw_el[pos_relax[0]]);
						h2_relax.push_back(rh_el[pos_relax[0]]);
					}
					if(el ==3){
						f3_relax.push_back(rf_el[pos_relax[0]]);
						w3_relax.push_back(rw_el[pos_relax[0]]);
						h3_relax.push_back(rh_el[pos_relax[0]]);
					}
				} else{
					std::cout << "Error when preparing the vector of parameters using the MCMC file" << std::endl;
					std::cout << "The uniqueness of the frequency " << inputs_local.eigen_params(pos_el[en],1) << " is not respected" << std::endl;
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
		} 
		// Empty the temporary variables for the relax and the freqs_ref at a given el
		f_el.resize(0);
		rf_el.resize(0);
		rw_el.resize(0);
		rh_el.resize(0);
	}
	// -------------------------------------------------------------------------------------------		
	// ------------------- Filtering according to specified the frequency range ------------------
	// -------------------------------------------------------------------------------------------
	
	if (Nf_el[0] != 0){
		f_el=f0_inputs;
		f0_inputs=filter_range(f0_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h0_inputs=filter_range(h0_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w0_inputs=filter_range(w0_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f0_relax=filter_range(f0_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h0_relax=filter_range(h0_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w0_relax=filter_range(w0_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f0_priors_min=filter_range(f0_priors_min, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f0_priors_max=filter_range(f0_priors_max, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
	}
	if (Nf_el[1] != 0){
		f_el=f1_inputs;
		f1_inputs=filter_range(f1_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h1_inputs=filter_range(h1_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w1_inputs=filter_range(w1_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f1_relax=filter_range(f1_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h1_relax=filter_range(h1_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w1_relax=filter_range(w1_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f1_priors_min=filter_range(f1_priors_min, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f1_priors_max=filter_range(f1_priors_max, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
	}
	if (Nf_el[2] != 0){
		f_el=f2_inputs;
		f2_inputs=filter_range(f2_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h2_inputs=filter_range(h2_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w2_inputs=filter_range(w2_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f2_relax=filter_range(f2_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h2_relax=filter_range(h2_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w2_relax=filter_range(w2_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f2_priors_min=filter_range(f2_priors_min, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f2_priors_max=filter_range(f2_priors_max, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
	}
	if (Nf_el[3] != 0){
		f_el=f3_inputs;
		f3_inputs=filter_range(f3_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h3_inputs=filter_range(h3_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w3_inputs=filter_range(w3_inputs, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f3_relax=filter_range(f3_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		h3_relax=filter_range(h3_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		w3_relax=filter_range(w3_relax, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f3_priors_min=filter_range(f3_priors_min, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
		f3_priors_max=filter_range(f3_priors_max, f_el, inputs_local.freq_range[0], inputs_local.freq_range[1]);
	}
	// Resizing the pointers for the number of modes
	Nf_el[0]=f0_inputs.size();
	Nf_el[1]=f1_inputs.size();
	Nf_el[2]=f2_inputs.size();
	Nf_el[3]=f3_inputs.size();

	if ((h0_inputs.size() + h1_inputs.size() + h2_inputs.size() + h3_inputs.size()) == 0){
		std::cout << " -------------------------------------------------------------------" << std::endl;
		std::cout << " No parameters found in the specified frequency range!" << std::endl;
		std::cout << " Cannot perform a fit: The frequency range(s) may not include any frequencies listed during the initial peak bagging" << std::endl;
		std::cout << " The program must exit now " << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// ------------------------------------------------------------------------------------------
	// ------------------------------- Handling the Common parameters ---------------------------
	// ------------------------------------------------------------------------------------------
	if( do_amp){
		std::cout << "   ===> Requested to fit squared amplitudes instead of Height... Converting height inputs into A_squared = pi*Height*Width..." << std::endl;
		tmpstr_h="Amplitude_l";
		for(int i=0; i<h0_inputs.size(); i++){
			h0_inputs[i]=pi*w0_inputs[i]*h0_inputs[i]; 
		}
		for(int i=0; i<h1_inputs.size(); i++){
 			h1_inputs[i]=pi*w1_inputs[i]*h1_inputs[i]; 
		}
		for(int i=0; i<h2_inputs.size(); i++){		
			h2_inputs[i]=pi*w2_inputs[i]*h2_inputs[i]; 
		}
		for(int i=0; i<h3_inputs.size(); i++){		
			h3_inputs[i]=pi*w3_inputs[i]*h3_inputs[i]; 
		}
	} else{
		tmpstr_h="Height_l";
	}
	// ------------------------ 		
    // Set default value of priors for Height Width and frequency
	tmpXd.resize(4);
	tmpXd << Hmin, Hmax, -9999., -9999.; // default hmin and hmax for the Jeffreys prior

	if(all_in.model_fullname == "model_MS_local_Hnlm"){
		io_calls.initialise_param(&height_in, Nf_el[0] + Nf_el[1]*2 + Nf_el[2]*3 + Nf_el[3]*4, Nmax_prior_params, -1, -1);
		p0=0;
		io_calls.fill_param_vect(&height_in, h0_inputs, h0_relax, tmpstr_h, "Jeffreys", tmpXd, p0, 0, 0);
		// NOTE; IN THAT CASE WE WAIT THE INCLINATION INFORMATION TO FILL Hnlm for l>0
	} else{
		io_calls.initialise_param(&height_in, Nf_el.sum(), Nmax_prior_params, -1, -1);
		p0=0;
		io_calls.fill_param_vect(&height_in, h0_inputs, h0_relax, tmpstr_h, "Jeffreys", tmpXd, p0, 0, 0);
		p0=h0_inputs.size();
		io_calls.fill_param_vect(&height_in, h1_inputs, h1_relax, tmpstr_h, "Jeffreys", tmpXd, p0, 0, 0);
		p0=h0_inputs.size()+h1_inputs.size();
		io_calls.fill_param_vect(&height_in, h2_inputs, h2_relax, tmpstr_h, "Jeffreys", tmpXd, p0, 0, 0);
		p0=h0_inputs.size()+h1_inputs.size()+h2_inputs.size();
		io_calls.fill_param_vect(&height_in, h3_inputs, h3_relax, tmpstr_h, "Jeffreys", tmpXd, p0, 0, 0);
	}
	
	io_calls.initialise_param(&width_in, Nf_el.sum(), Nmax_prior_params, -1, -1);
	io_calls.initialise_param(&freq_in, Nf_el.sum(), Nmax_prior_params, -1, -1); 

	// --- Default setup for widths ---
	tmpXd.resize(4);
	if (inputs_local.Dnu > 0){
    	tmpXd << resol, inputs_local.Dnu/3., -9999., -9999.;
	} else{
		std::cout << "Warning: No large separation provided. The default maximum width will be set to 20 microHz " << std::endl;
		std::cout << "         If you want to set this yourself, use the 'width' keyword for the section '# Controls and priors for common parameters'" << std::endl;
		std::cout << "         of the .model file" << std::endl;
		tmpXd << resol, 20., -9999., -9999.;
	}
	p0=0;
	io_calls.fill_param_vect(&width_in, w0_inputs, w0_relax, "Width_l", "Jeffreys", tmpXd, p0, 0, 0);
	p0=w0_inputs.size();
	io_calls.fill_param_vect(&width_in, w1_inputs, w1_relax, "Width_l", "Jeffreys", tmpXd, p0, 0, 0);
	p0=w0_inputs.size()+w1_inputs.size();
	io_calls.fill_param_vect(&width_in, w2_inputs, w2_relax, "Width_l", "Jeffreys", tmpXd, p0, 0, 0);
	p0=w0_inputs.size()+w1_inputs.size()+w2_inputs.size();
	io_calls.fill_param_vect(&width_in, w3_inputs, w3_relax,"Width_l", "Jeffreys", tmpXd, p0, 0, 0);
	
	// --- Default setup for frequencies ---
	p0=0;
	for(int i=0; i<f0_inputs.size(); i++){
		if(f0_relax[i]){
			tmpXd << f0_priors_min[i], f0_priors_max[i], 0.01*std::abs(f0_priors_max[i]-f0_priors_min[i]), 0.01*std::abs(f0_priors_max[i]-f0_priors_min[i]); // default parameters for a GUG prior on frequencies
			io_calls.fill_param(&freq_in, "Frequency_l", "GUG", f0_inputs[i], tmpXd, i+p0, 0);	
		} else{
			io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f0_inputs[i], tmpXd, i+p0, 0);			
		}
	}

	p0=f0_inputs.size();
	for(int i=0; i<f1_inputs.size(); i++){
		if(f1_relax[i]){
			tmpXd << f1_priors_min[i], f1_priors_max[i], 0.01*std::abs(f1_priors_max[i]-f1_priors_min[i]), 0.01*std::abs(f1_priors_max[i]-f1_priors_min[i]); // default parameters for a GUG prior on frequencies
			io_calls.fill_param(&freq_in, "Frequency_l", "GUG", f1_inputs[i], tmpXd, i+p0, 0);	
		} else{
			io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f1_inputs[i], tmpXd, i+p0, 0);			
		}
	}

	p0=f0_inputs.size()+f1_inputs.size();
	for(int i=0; i<f2_inputs.size(); i++){
		if(f2_relax[i]){
			tmpXd << f2_priors_min[i], f2_priors_max[i], 0.01*std::abs(f2_priors_max[i]-f2_priors_min[i]), 0.01*std::abs(f2_priors_max[i]-f2_priors_min[i]); // default parameters for a GUG prior on frequencies
			io_calls.fill_param(&freq_in, "Frequency_l", "GUG", f2_inputs[i], tmpXd, i+p0, 0);	
		} else{
			io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f2_inputs[i], tmpXd, i+p0, 0);			
		}
	}

	p0=f0_inputs.size()+f1_inputs.size()+f2_inputs.size();
	for(int i=0; i<f3_inputs.size(); i++){
		if(f3_relax[i]){
			tmpXd << f3_priors_min[i], f3_priors_max[i], 0.01*std::abs(f3_priors_max[i]-f3_priors_min[i]), 0.01*std::abs(f3_priors_max[i]-f3_priors_min[i]); // default parameters for a GUG prior on frequencies
			io_calls.fill_param(&freq_in, "Frequency_l", "GUG", f3_inputs[i], tmpXd, i+p0, 0);	
		} else{
			io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f3_inputs[i], tmpXd, i+p0, 0);			
		}
	}

	// ----- Switch between the models that handle averaging over n,l or both -----
   if(do_a11_eq_a12 == 1 && do_avg_a1n == 1){
        io_calls.initialise_param(&Snlm_in, 6, Nmax_prior_params, -1, -1);
    }  
    if(do_a11_eq_a12 == 0 && do_avg_a1n == 1){
        io_calls.initialise_param(&Snlm_in, 7, Nmax_prior_params, -1, -1);
    }
    if(do_a11_eq_a12 == 1 && do_avg_a1n == 0){
		if (Nf_el[1] == Nf_el[2]){
			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1], Nmax_prior_params, -1, -1);
		} else {
			std::cout << "When considering a11=a22"<< std::endl;
			std::cout <<" You must have as many l=1 than l=2" << std::endl;
			std::cout <<" Here we have: Nf(l=1) = " << Nf_el[1] << " and Nf(l=2) = " << Nf_el[2] << std::endl;
			std::cout <<" Check the .model file" << std::endl;
			std::cout <<" The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
    }
    if(do_a11_eq_a12 == 0 && do_avg_a1n == 0){
    	io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2], Nmax_prior_params, -1, -1);
    }
      
	// -------------- Set Extra_priors ----------------	
	extra_priors.resize(4);
	extra_priors[0]=0; // Empty slot
	extra_priors[1]=0; // Empty slot
	extra_priors[2]=0.2; // By default a3/a1<=1
	extra_priors[3]=0; // Default is not imposing Sun(Hnlm)_{l=-m, l=+m} = 1
	// ------------------------------------------------
	
	for(int i=0; i<inputs_local.common_names.size(); i++){
		// --- Common parameters than can be run during setup ---
		if(inputs_local.common_names[i] == "freq_smoothness" || inputs_local.common_names[i] == "Freq_smoothness"){
			std::cout << "    freq_smoothness keyword found but the analysis performs a local fit (use of io_local.cpp)" << std::endl;
			std::cout << "    For a local fit, freq_smoothness is irrelevant.  The keyword will be skipped." << std::endl;
			std::cout << "    If you intended to perform a fit on a large range, please consider using io_MS_Global instead of io_local" << std::endl;
			std::cout << "    Pursuing..." << std::endl; 
		}
		if(inputs_local.common_names[i] == "Visibility_l1" || inputs_local.common_names[i] == "Visibility_l2" || inputs_local.common_names[i] == "Visibility_l3"){
			std::cout << "    Visibility_lx (x E [1,2,3]) keyword found but the analysis performs a local fit (use of io_local.cpp)" << std::endl;
			std::cout << "    For a local fit, Heights of modes are directly fitted.  The keyword will be skipped." << std::endl;
			std::cout << "    If you intended to perform a fit on a large range, please consider using io_MS_Global instead of io_local" << std::endl;
			std::cout << "    Pursuing..." << std::endl; 
		}
		if(inputs_local.common_names[i] == "trunc_c"){
			if(inputs_local.common_names_priors[i] != "Fix"){
				fatalerror_msg_io_local("trunc_c", "Fix", "[Truncation parameter]", "20" );
			} else{
				trunc_c=inputs_local.modes_common(i,0); // This value is added to the input vector at the end of this function.
			}
		}	

		// --- Frequencies ---
//  ----- WE REMOVE THE POSSIBILITY OF CHANGING THE DEFAULT PRIRO ON THE FREQUENCY FROM GUG TO UNIFORM ------
//  ----- THIS WOULD ADD TOO MANY LINES OF CODES WITH THE CURRENT CODING CHOICE (f0, f1, etc... all a independent vector subject to high code redundancy) -----
// ------- THE CODING CHOICE HAS TO BE IMPROVED BEFORE THINKING OF THAT ----
/*		if(inputs_local.common_names[i] == "Frequency" || inputs_local.common_names[i] == "frequency"){ 
			if(inputs_local.common_names_priors[i] == "GUG" || inputs_local.common_names_priors[i] == "Uniform"){
				for(int p0=0; p0<f_inputs.size(); p0++){
					if(f_relax[p0]){
						if(inputs_local.common_names_priors[i] == "GUG"){
							tmpXd << f_priors_min[p0], f_priors_max[p0], inputs_local.modes_common(i,3), inputs_local.modes_common(i,4);
						} else{
							tmpXd << f_priors_min[p0], f_priors_max[p0], -9999, -9999; 						
						}
						io_calls.fill_param(&freq_in, "Frequency_l", inputs_local.common_names_priors[i], f_inputs[p0], tmpXd, p0, 0);	
					} else{
						io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f_inputs[p0], tmpXd, p0, 0);			
					}
				}
			} else{
				fatalerror_msg_io_local(inputs_local.common_names[i], "GUG or Uniform", "", "" );
			}
		}
*/		
		std::cout << "i=[" << i << "]" << std::endl;
		// --- Height or Amplitude ---
		if(inputs_local.common_names[i] == "height" || inputs_local.common_names[i] == "Height" || 
		   inputs_local.common_names[i] == "amplitude" || inputs_local.common_names[i] == "Amplitude"){

			if(inputs_local.common_names_priors[i] == "Fix_Auto"){
				//In auto mode, A Jeffrey's prior is used, with an upper bound given by the input value multiplied by a factor X=inputs_local.modes_common(i,1) given by the user (typically 3 is fine). 
				//The lower bound is fixed as Y=inputs_local.modes_common(i,0) times lower than the actual input.
				//The reason is that for RGBs, a single common value of the upper bound for the Jeffrey prior lead to a too large dynamic range for low SNR modes (e.g. low freq or l=3).
				//In the .model file, the new prior can be set by setting e.g. (here X=3 and Y=10):

					//fatalerror_msg_io_local(inputs_local.common_names[i], "Fix_Auto", "", "" );
					tmpXd_mat.resize(h0_inputs.size(),4);
					for(int k=0; k<h0_inputs.size(); k++){
						if (inputs_local.common_names[i] == "amplitude" || inputs_local.common_names[i] == "Amplitude"){
							tmpXd_mat.row(k) << pi*inputs_local.Dnu/3.*h0_inputs[k]/ inputs_local.modes_common(i,0), pi*inputs_local.Dnu/3.*h0_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
						} else{
							tmpXd_mat.row(k) << h0_inputs[k]/ inputs_local.modes_common(i,0), h0_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
						}
					}
					if(all_in.model_fullname == "model_MS_local_Hnlm"){
						pos_prior_height=i;
						io_calls.fill_param_vect2(&height_in, h0_inputs, h0_relax, tmpstr_h, "Jeffreys", tmpXd_mat, p0, 0, 0);
					} else{
						p0=0;
						io_calls.fill_param_vect2(&height_in, h0_inputs, h0_relax, tmpstr_h, "Jeffreys", tmpXd_mat, p0, 0, 0);
						p0=h0_inputs.size();
						tmpXd_mat.resize(h1_inputs.size(),4);
						for(int k=0; k<h1_inputs.size(); k++){
							if (inputs_local.common_names[i] == "amplitude" || inputs_local.common_names[i] == "Amplitude"){
								tmpXd_mat.row(k) << pi*inputs_local.Dnu/3.*h1_inputs[k]/ inputs_local.modes_common(i,0), pi*inputs_local.Dnu/3.*h1_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
							} else{
								tmpXd_mat.row(k) << h1_inputs[k]/ inputs_local.modes_common(i,0), h1_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
							}
						}
						io_calls.fill_param_vect2(&height_in, h1_inputs, h1_relax, tmpstr_h,  "Jeffreys", tmpXd_mat, p0, 0, 0);
						p0=h0_inputs.size()+h1_inputs.size();
						tmpXd_mat.resize(h2_inputs.size(),4);
						for(int k=0; k<h2_inputs.size(); k++){
							if (inputs_local.common_names[i] == "amplitude" || inputs_local.common_names[i] == "Amplitude"){
								tmpXd_mat.row(k) << pi*inputs_local.Dnu/3.*h2_inputs[k]/ inputs_local.modes_common(i,0), pi*inputs_local.Dnu/3.*h2_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
							} else{
								tmpXd_mat.row(k) << h2_inputs[k]/ inputs_local.modes_common(i,0), h2_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
							}
						}
						io_calls.fill_param_vect2(&height_in, h2_inputs, h2_relax, tmpstr_h,  "Jeffreys", tmpXd_mat, p0, 0, 0);
						p0=h0_inputs.size()+h1_inputs.size()+h2_inputs.size();
						tmpXd_mat.resize(h3_inputs.size(),4);
						for(int k=0; k<h3_inputs.size(); k++){
							if (inputs_local.common_names[i] == "amplitude" || inputs_local.common_names[i] == "Amplitude"){
								tmpXd_mat.row(k) << pi*inputs_local.Dnu/3.*h3_inputs[k]/ inputs_local.modes_common(i,0), pi*inputs_local.Dnu/3.*h3_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
							} else{
								tmpXd_mat.row(k) << h3_inputs[k]/ inputs_local.modes_common(i,0), h3_inputs[k]* inputs_local.modes_common(i,1) , -9999, -9999;
							}
						}						
						io_calls.fill_param_vect2(&height_in, h3_inputs, h3_relax, tmpstr_h,  "Jeffreys", tmpXd_mat, p0, 0, 0);
					}
			} else{
	
				if(all_in.model_fullname == "model_MS_local_Hnlm"){
					pos_prior_height=i;
					io_calls.fill_param_vect(&height_in, h0_inputs, h0_relax, tmpstr_h, inputs_local.common_names_priors[i], inputs_local.modes_common.row(i), p0, 0, 0);
				} else{
					p0=0;
					io_calls.fill_param_vect(&height_in, h0_inputs, h0_relax, tmpstr_h, inputs_local.common_names_priors[i], inputs_local.modes_common.row(i), p0, 1, 1);
					p0=h0_inputs.size();
					io_calls.fill_param_vect(&height_in, h1_inputs, h1_relax, tmpstr_h, inputs_local.common_names_priors[i], inputs_local.modes_common.row(i), p0, 1, 1);
					p0=h0_inputs.size()+h1_inputs.size();
					io_calls.fill_param_vect(&height_in, h2_inputs, h2_relax, tmpstr_h, inputs_local.common_names_priors[i], inputs_local.modes_common.row(i), p0, 1, 1);
					p0=h0_inputs.size()+h1_inputs.size()+h2_inputs.size();
					io_calls.fill_param_vect(&height_in, h3_inputs, h3_relax, tmpstr_h, inputs_local.common_names_priors[i], inputs_local.modes_common.row(i), p0, 1, 1);
				}
			}
		}
		// -- Mode Width ---
		if(inputs_local.common_names[i] == "width" || inputs_local.common_names[i] == "Width"){
				if(inputs_local.common_names_priors[i] == "Fix_Auto"){
                   std::cout << "Fix_Auto requested for Widths... For all free Widths, the prior will be with this syntax:" << std::endl;
                   tmpstr="Jeffreys";
                    if(inputs_local.Dnu > 0){
                    	std::cout << "          " << std::left << std::setw(15) << tmpstr << " [Spectrum Resolution]   [Deltanu / 3]   -9999    -9999" << std::endl;
	                    tmpXd << resol, inputs_local.Dnu/3, -9999, -9999;
                    } else{
                    	std::cout << "          " << std::left << std::setw(15) << tmpstr << " [Spectrum Resolution]   20   -9999    -9999" << std::endl;                   	
                    	tmpXd << resol, 20, -9999, -9999;

                    }
                    std::cout << "          " << "Resolution: " << resol << std::endl;
                } else{
                	std::cout << "Using user-defined width prior" << std::endl;
                    tmpstr=inputs_local.common_names_priors[i];
                    tmpXd=inputs_local.modes_common.row(i);
                }
            p0=0;
            io_calls.fill_param_vect(&width_in, w0_inputs, w0_relax, "Width_l", tmpstr, tmpXd, p0, 0, 1);
            p0=w0_inputs.size();
            io_calls.fill_param_vect(&width_in, w1_inputs, w1_relax, "Width_l", tmpstr, tmpXd, p0, 0, 1);
            p0=w0_inputs.size()+w1_inputs.size();
            io_calls.fill_param_vect(&width_in, w2_inputs, w2_relax, "Width_l", tmpstr, tmpXd, p0, 0, 1);
            p0=w0_inputs.size()+w1_inputs.size()+w2_inputs.size();
            io_calls.fill_param_vect(&width_in, w3_inputs, w3_relax, "Width_l", tmpstr, tmpXd, p0, 0, 1);
            //exit(EXIT_SUCCESS);
		} 	
		// --- Splittings and asymetry ---
		if(inputs_local.common_names[i] == "splitting_a1" || inputs_local.common_names[i] == "Splitting_a1"){
			if(inputs_local.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_local("splitting_a1", "Fix_Auto", "", "" );
			}
			p0=0;
			io_calls.fill_param(&Snlm_in, "Splitting_a1", inputs_local.common_names_priors[i], inputs_local.modes_common(i,0), inputs_local.modes_common.row(i), p0, 1);	
        	if(do_a11_eq_a12 == 0 && do_avg_a1n == 1){ // Case a1(l), we put the same prior for a1(2) than a1(1)
		    	std::cout << "do_a11_eq_a12 == 0 && do_avg_a1n == 1" << std::endl;
		    	p0=6;
        	    io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_local.modes_common.row(i), p0, 1);
        	}
        	if(do_a11_eq_a12 == 1 && do_avg_a1n == 0){ // Case a1(n), we put the same prior for a1(n)
				std::cout << "do_a11_eq_a12 == 1 && do_avg_a1n == 0" << std::endl;
				for(int kk=0; kk<Nf_el[1]; kk++){
					p0=6+kk;
					io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_local.modes_common.row(i), p0, 1); 
				}
				p0=0;
				io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, inputs_local.modes_common.row(i), p0, 1);
        	}
        	if(do_a11_eq_a12 == 0 && do_avg_a1n == 0){ // Case a1(n,l), we put the same prior for a1(n,l)
				std::cout << "do_a11_eq_a12 == 0 && do_avg_a1n == 0" << std::endl;
				for(int kk=0; kk<Nf_el[1]+Nf_el[2]; kk++){
					p0=6+kk;
					io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_local.modes_common.row(i), p0, 1); 
				}
				p0=0;
				io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, inputs_local.modes_common.row(i), p0, 1); // Erase values of the default Splitting_a1 block
        	}
		}
		if(inputs_local.common_names[i] == "asphericity_eta"|| inputs_local.common_names[i] == "Asphericity_eta"){  
			Snlm_in.inputs_names[1]="Asphericity_eta0";
			if(inputs_local.common_names_priors[i] == "Fix_Auto"){ // Case where the centrifugal force is added but FIXED.
				Snlm_in.priors_names[1]="Fix"; //In all cases the centrifugal force is fixed
				Snlm_in.relax[1]=0;
				if(inputs_local.modes_common(i,0) == 1){
  					if(inputs_local.Dnu > 0){
						Snlm_in.inputs[1]=3./(4.*M_PI*rho*G); // rho is calculated using Dnu
					} else{
						std::cout << "Centrifugal force set to 0 because the user did not provide Dnu in the model file" << std::endl;
						Snlm_in.inputs[1]=0; // If Dnu was not provided, we must fix the centrifugal distorsion constant part to 0.
					}
					std::cout << " -------------" << std::endl;
					std::cout << "Asphericity_eta given with Fix_Auto ==> Setting the asphericity to the centrifugal force" << std::endl;
					std::cout << "      eta0=" << Snlm_in.inputs[1] << std::endl;
					std::cout << " -------------" << std::endl;
				} else{
					Snlm_in.inputs[1]=0;
				}
			} else{ // Case where the centrifugal force is user-defined: Could be free or fixed
				p0=1;
				std::cout << " WARNING: Asphericity_eta is set to a non-Fix_Auto value: It means that the asphericity parameter is not set to the Centrifugal force automatically "<< std::endl;
				std::cout << "          If this is not your intend, please set Asphericity_eta to Fix_Auto" << std::endl;
				io_calls.fill_param(&Snlm_in, "Asphericity_eta", inputs_local.common_names_priors[i], inputs_local.modes_common(i,0), inputs_local.modes_common.row(i), p0, 1);	
			}
		}

		if(inputs_local.common_names[i] == "splitting_a3" || inputs_local.common_names[i] == "Splitting_a3"){
			if(inputs_local.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_local("splitting_a3", "Fix_Auto", "", "" );
			}
			p0=2;
			io_calls.fill_param(&Snlm_in, "Splitting_a3", inputs_local.common_names_priors[i], inputs_local.modes_common(i,0), inputs_local.modes_common.row(i), p0, 1);
		}

		if(inputs_local.common_names[i] == "asymetry" || inputs_local.common_names[i] == "Asymetry"){
			if(inputs_local.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_local("asymetry", "Fix_Auto", "", "" );
			}
			p0=5;
			io_calls.fill_param(&Snlm_in, "Lorentzian_asymetry", inputs_local.common_names_priors[i], inputs_local.modes_common(i,0), inputs_local.modes_common.row(i), p0, 1);
		}
		// --- Finally the inclination ---
		if(inputs_local.common_names[i] == "inclination" || inputs_local.common_names[i] == "Inclination"){
			if(inputs_local.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_local("inclination", "Fix_Auto", "", "" );
			}
			if( inputs_local.modes_common(i,0) >= 90){ // Avoid some nan when computing the prior
				    tmp=89.99999; 
			} else{
				tmp=inputs_local.modes_common(i,0);
			}
			p0=0;
			io_calls.fill_param(&Inc_in, "Inclination", inputs_local.common_names_priors[i], tmp, inputs_local.modes_common.row(i), p0, 1);
		}
		
		if(inputs_local.common_names[i] == "sqrt(splitting_a1).cosi"){
			if(inputs_local.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_local("sqrt(splitting_a1).cosi", "Fix_Auto", "", "" );
			}
			p0=3;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", inputs_local.common_names_priors[i], inputs_local.modes_common(i,0), inputs_local.modes_common.row(i), p0, 1);

            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0  and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1cosi=1;
		}
		
		if(inputs_local.common_names[i] == "sqrt(splitting_a1).sini"){
			if(inputs_local.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_local("sqrt(splitting_a1).sini", "Fix_Auto", "", "" );
			}
			p0=4;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", inputs_local.common_names_priors[i], inputs_local.modes_common(i,0), inputs_local.modes_common.row(i), p0, 1);
			
            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0 and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1sini=1;
		}
	}

if(bool_a1cosi != bool_a1sini){ // Case when one of the projected splitting quantities is missing ==> Problem
	std::cout << "Warning: Both 'sqrt(splitting_a1).sini' and 'sqrt(splitting_a1).cosi' keywords must appear" << std::endl;
	std::cout << "         It is forbidden to use only one of them" << std::endl;
	std::cout << "         Edit the .MODEL file accordingly" << std::endl;
	std::cout << "The program will exit now" << std::endl;
	exit(EXIT_FAILURE);
}
if((bool_a1cosi == 0) && (bool_a1sini == 0)){ // Case where Inclination and Splitting_a1 are supposed to be used 
	if( all_in.model_fullname == "model_MS_local_basic"){
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
        io_calls.fill_param(&Inc_in, "Empty", "Fix", 0, Inc_in.priors.col(0), 0, 1); // Note that inputs_local.modes_common.row(0) is not used... just dummy
        io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, Snlm_in.priors.col(0), 0, 1); // Note that inputs_local.modes_common.row(0) is not used... just dummy
	}

	if(all_in.model_fullname == "model_MS_local_Hnlm"){	
	
		tmpXd.resize(4);
		tmpXd << -9999, -9999, -9999., -9999.;

		tmp=Inc_in.inputs[0]; // Recover the initial stellar inclination as defined by the user in the .model file

		// Visibilities are not used in the is model ==> deactivated
        io_calls.fill_param(&Inc_in, "Empty", "Fix", 0, Inc_in.priors.col(0), 0, 1); // Note that inputs_local.modes_common.row(0) is not used... just dummy
		ind=h0_inputs.size();
		if(pos_prior_height <= 0){
			tmpXd << Hmin, Hmax, -9999., -9999.;
			tmpstr="Jeffreys";
		} else{
			tmpXd=inputs_local.modes_common.row(pos_prior_height);
			tmpstr=inputs_local.common_names_priors[pos_prior_height];
		}
		for(int el=1; el<=lmax; el++){
			ratios_l=amplitude_ratio(el, tmp); 
			for(int en=0; en<Nf_el[el]; en++){
				for(int em=0; em<=el; em++){
					switch(el){
						case 1:
							io_calls.fill_param(&height_in, "H(" + int_to_str(en) + "," + int_to_str(el) + "," + int_to_str(em) + ")", tmpstr, h1_inputs[en]*ratios_l[el+em], tmpXd, ind, 0);
							ind=ind+1;
							break;
						case 2:
							io_calls.fill_param(&height_in, "H(" + int_to_str(en) + "," + int_to_str(el) + "," + int_to_str(em) + ")", tmpstr, h2_inputs[en]*ratios_l[el+em], tmpXd, ind, 0);
							ind=ind+1;
							break;
						case 3:
							io_calls.fill_param(&height_in, "H(" + int_to_str(en) + "," + int_to_str(el) + "," + int_to_str(em) + ")", tmpstr, h3_inputs[en]*ratios_l[el+em], tmpXd, ind, 0);
							ind=ind+1;
							break;
					}
				}
			}
		}
		extra_priors[3]=2; // Impose Sum(H(nlm))_{l=-m,l=+m} =1 for that model (case == 2 of priors_local())
		//exit(EXIT_SUCCESS);
	}
}
if((bool_a1cosi == 1) && (bool_a1sini ==1)){
    std::cout << "    NOTICE: sqrt(splitting_a1).sini and sqrt(splitting_a1).cosi superseed inclination and splitting" << std::endl;
    std::cout << "            Slots for the variables Inclination and Splitting are therefore forced to be FIXED" << std::endl;
    std::cout << "            Be aware that may lead to unwished results if you are not careful about the used model" << std::endl;
   	
   	io_calls.fill_param(&Inc_in, "Empty", "Fix", 0, Inc_in.priors.col(0), 0, 1); // Note that inputs_local.modes_common.row(0) is not used... just dummy   	
   	io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, Snlm_in.priors.col(0), 0, 1); // "Splitting_a1" default values are erased
}
	// ----------------------------------------------------
	// ---------------- Handling noise --------------------
	// ----------------------------------------------------
	io_calls.initialise_param(&Noise_in, 1, Nmax_prior_params, -1, -1);
	set_noise_params_local(&Noise_in, inputs_local.noise_s2, inputs_local.noise_params, inputs_local.freq_range); 	// Set defaults

	// -------------------------------------------------------------------------------------
	// ------------- Sticking everything together in a Input_Data structure --------------
	// -------------------------------------------------------------------------------------
	
	std::cout << "Setting plength..." << std::endl;
	plength.resize(11);
	if(all_in.model_fullname == "model_MS_local_Hnlm"){
		plength[0]=h0_inputs.size()+2*h1_inputs.size()+3*h2_inputs.size()+4*h3_inputs.size(); 
	} else{
		plength[0]=h0_inputs.size()+h1_inputs.size()+h2_inputs.size()+h3_inputs.size(); 		
	}
	plength[1]=0; // Local models do not have visibilities // lmax                  ; 
	plength[2]=Nf_el[0];            plength[3]=Nf_el[1]       ; plength[4]=Nf_el[2]		   ; plength[5]=Nf_el[3];
	plength[6]=Snlm_in.inputs.size(); 
	plength[7]=w0_inputs.size()+w1_inputs.size()+w2_inputs.size()+w3_inputs.size() ; 
	plength[8]=Noise_in.inputs.size(); 
	plength[9]=Inc_in.inputs.size();
	plength[10]=2; // This is trunc_c and do_amp;
	
	std::cout << "init all_in..." << std::endl;	
	io_calls.initialise_param(&all_in, plength.sum(), Nmax_prior_params, plength, extra_priors);
	
	// --- Put the Height or Amplitudes---
	std::cout << "Setting heights..." << std::endl;
	io_calls.show_param(height_in,1);
	std::cout << " ----" << std::endl;
	p0=0;
	io_calls.add_param(&all_in, &height_in, p0); // height_in may contain either height or amplitudes depending on specified keywords
	// --- Put the Visibilities --- // NOT USED FOR LOCAL FIT
	//std::cout << "Setting visibilities..." << std::endl;	
	//p0=all_in.plength[0];
	//io_calls.add_param(&all_in, &Vis_in, p0);
	
	// --- Put the Frequencies ---
	//std::cout << "Setting freqs..." << std::endl;
	p0=all_in.plength[0] + all_in.plength[1];	
	io_calls.add_param(&all_in, &freq_in, p0);

	// --- Put the Snlm (splittings and asymetry) ---
	//std::cout << "Setting asym..." << std::endl;
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5];
	io_calls.add_param(&all_in, &Snlm_in, p0);
	
	// --- Put the Width ---
	//std::cout << "Setting widths..." << std::endl;
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] +  all_in.plength[3]  + all_in.plength[4] + all_in.plength[5] + all_in.plength[6];
	io_calls.add_param(&all_in, &width_in, p0);
	
	// --- Put the Noise ---
	//std::cout << "Setting noise..." << std::endl;
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] +  all_in.plength[3]  + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7];
	io_calls.add_param(&all_in, &Noise_in, p0);
	
	// --- Put the Inclination ---
	//std::cout << "Setting inc..." << std::endl;	
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8];
	io_calls.add_param(&all_in, &Inc_in, p0);
	
	// --- Add trunc_c that controls the truncation of the Lorentzian ---
	//std::cout << "Setting trunc..." << std::endl;
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9];
	io_calls.fill_param(&all_in, "Truncation parameter", "Fix", trunc_c, inputs_local.modes_common.row(0), p0, 1);
	if (all_in.inputs[p0] <= 0){
		std::cout << "Warning: trunc_c <= 0. This is forbidden. Setting default to 10000. (No truncation)" << std::endl;
		all_in.inputs[p0]=10000.; // In case of a non-sense value for c, we use Full-Lorentzian as default
	}
	// -- Add the Amplitude switch --
	//std::cout << "Setting amp_switch..." << std::endl;
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9] + 1;
	io_calls.fill_param(&all_in, "Switch for fit of Amplitudes or Heights", "Fix", do_amp, inputs_local.modes_common.row(0), p0,1);
					
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
		std::cout << "          - " << all_in.plength[3] << "  l=1 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[4] << "  l=2 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[5] << "  l=3 Frequencies parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[6] << "  Splitting parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[7] << "  Width parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[8] << "  Noise parameters" << std::endl; 
		std::cout << "          - " << all_in.plength[9] << "  Inclination parameters" << std::endl; 
		std::cout << std::endl << "          - " << "Total Number of Parameters: " << all_in.plength.sum() << std::endl; 
		std::cout << " -----------------------------------------------------------" << std::endl;
	}
	
	//exit(EXIT_SUCCESS);
return all_in;
}

short int set_noise_params_local(Input_Data *Noise_in, const MatrixXd& noise_s2, const VectorXd& noise_params, const VectorXd& freq_range){
/*
 *
 * A function that prepares the Noise_in Data structure using 
 * inputs of the .model file regrouped inside the noise_s2 matrixXd
 *
*/
	int cpt=0;
	VectorXi pos;
	VectorXd noise_p;

	// On a local fit, we just need to approximate the noise by a constant... Later refinements may include a slope, 
	//but this might no be necessary
	(*Noise_in).inputs_names[0]="White_Noise_N0"; 
	(*Noise_in).priors_names[0]="Uniform"; 
	(*Noise_in).relax[0]=1; 
	
	int Nharvey;
	VectorXd y0(2), noise_vals;
	pos=where_dbl(noise_params, -2, 1e-6);

	if(pos[0] != -1){
		noise_p.resize(noise_params.size() -pos.size());
	}
	else{
		noise_p.resize(noise_params.size());	
	}
	noise_p.setConstant(-9999);

	y0 << 0, 0;
	
	for (int i=0; i<noise_params.size(); i++){
		std::cout << i << "  " << cpt << std::endl;
		if (noise_params[i] ==-1){
				noise_p[cpt]=0;
				cpt=cpt+1;
		}
		if (noise_params[i] >=0){
				noise_p[cpt] = noise_params[i];
				cpt=cpt+1;
		}
	}

	Nharvey=(noise_p.size()-1)/3;
	noise_vals=harvey_like(noise_p, freq_range, y0, Nharvey);

	(*Noise_in).inputs[0]=noise_vals.sum()/noise_vals.size();
	(*Noise_in).priors(0,0)=noise_vals.minCoeff()*0.5;
	(*Noise_in).priors(1,0)=noise_vals.maxCoeff()*1.5;	 
	
///	 --- For debug only ---
//	std::cout << "Stop in set_noise_params_local: Need to be checked" << std::endl;	
//	std::cout << "Nharvey = " << Nharvey << std::endl;
//	std::cout << "freq_range = " << freq_range << std::endl;
//	std::cout << "noise_params = " <<  noise_p << std::endl;
//	std::cout << "------" << std::endl;
//	std::cout << "noise_vals:" << noise_vals << std::endl;
//	std::cout << "(*Noise_in).inputs : " << (*Noise_in).inputs << std::endl;
//	std::cout << "(*Noise_in).priors : " << (*Noise_in).priors << std::endl;
	return 0;
}

short int fatalerror_msg_io_local(const std::string varname, const std::string param_type, const std::string syntax_vals, const std::string example_vals){
/*
* Function that handle error messages and warnings for io_local
*/
	if(param_type == "Fix_Auto"){
		std::cout << "         Warning: Fix_Auto is not implemented for that parameter" << std::endl;
		std::cout << "         		You must choose the prior yourself" << std::endl;
	} else{
		std::cout << "         Warning: " << varname << " should always be defined as '" << param_type << "'" << std::endl;
	}
	if(syntax_vals != ""){
		std::cout << "                  The syntax should be as follow: " << varname << "    " << param_type << "    " << syntax_vals << std::endl; 
	}
	if(example_vals !=""){
		std::cout << "                  Example: " << varname << "    " << param_type << "    " << example_vals << std::endl; 
	}
	std::cout << "         The program will exit now" << std::endl;
	exit(EXIT_FAILURE);
	return -1;
}

