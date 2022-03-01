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
#include "io_ms_global.h"
#include "io_models.h"
#include "function_rot.h"
#include "interpol.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


MCMC_files read_MCMC_file_MS_Global(const std::string cfg_model_file, const bool verbose){

	bool range_done=0;
	int i, out, nl, el, cpt;
	std::vector<int> pos;
	std::string line0, char0, char1, char2;
	std::vector<std::string> word, tmp;
	std::ifstream cfg_session;
    MatrixXd tmpXd;

	MCMC_files iMS_global;

	iMS_global.numax=-9999; // Initialize the optional variable numax
	
    cpt=0;
    i=0;
    out=0;
    nl=200; // maximum number of lines
    el=4; // maximum degree of the modes
    //cfg_session.open(modeling.cfg_model_file.c_str());
    cfg_session.open(cfg_model_file.c_str());
    if (cfg_session.is_open()) {

	char0="#"; 
	std::getline(cfg_session, line0);
	if(verbose == 1) {std::cout << "  - Global parameters:" << std::endl;}
	iMS_global.freq_range.resize(2);
	iMS_global.els.resize(400);
	iMS_global.freqs_ref.resize(400);
	while((out < 3) && (!cfg_session.eof())){

		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1));
		char1=line0.substr(1, 1);
		if (char0 == "#" && char1 == "K"){
			word=strsplit(line0, "= \t");
			iMS_global.ID=strtrim(strtrim(word[1]));
			if(verbose == 1) {std::cout << "           ID=" << iMS_global.ID << std::endl;}
		}
		if (char0 == "!" && char1 == "n"){
			word=strsplit(line0, " ");
			iMS_global.numax=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           numax =" << iMS_global.numax << std::endl;}		
		}
		if (char0 == "!" && char1 != "!" && char1 != "n"){
			word=strsplit(line0, " ");
			iMS_global.Dnu=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           Dnu =" << iMS_global.Dnu << std::endl;}
		}
		if (char0 == "!" && char1 == "!"){
			word=strsplit(line0, " ");
			iMS_global.C_l=str_to_dbl(word[1]);
			if(verbose == 1) {std::cout << "           C_l =" << iMS_global.C_l << std::endl;}
		}
		if (char0 == "*"){
			if (range_done == 0){
				word=strsplit(line0, " ");
				iMS_global.freq_range[0]=str_to_dbl(word[1]);
				iMS_global.freq_range[1]=str_to_dbl(word[2]);
				if(verbose == 1) {std::cout << "           freq_range = [" << iMS_global.freq_range[0] << " , " << iMS_global.freq_range[1] << "]" << std::endl;}
				range_done=1;
			} else{
				std::cout << "Error: Multiple range detected. This is not allowed with io_ms_global models and priors" << std::endl;
				std::cout << "       Check you .model file and correct the configuration accordingly" << std::endl;
				std::cout << "       The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		if ((char0 != "#") && (char0 != "!") && (char0 != "*")){
			word=strsplit(line0, " ");
			if(strtrim(word[0]) == "p" || strtrim(word[0]) == "g" || strtrim(word[0]) == "co"){
				iMS_global.param_type.push_back(strtrim(word[0]));
				iMS_global.els[cpt]=str_to_int(word[1]);
                iMS_global.freqs_ref[cpt]=str_to_dbl(word[2]);
                if (word.size() >=4){
                    iMS_global.relax_freq.push_back(str_to_bool(word[3]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    iMS_global.relax_freq.push_back(1);  // By default, we fit frequencies
                }
                if (word.size() >=5){
                    //iMS_global.relax_gamma.push_back(str_to_bool(word[4]));
                    iMS_global.relax_H.push_back(str_to_bool(word[4]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the Width of the mode at frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    //iMS_global.relax_gamma.push_back(1);
                    iMS_global.relax_H.push_back(1);
                }
                //std::cout << "H" << std::endl;
                if (word.size() >=6){
                    //iMS_global.relax_H.push_back(str_to_bool(word[5]));
                    iMS_global.relax_gamma.push_back(str_to_bool(word[5]));
                } else{
                    std::cout << "Warning: The .model file does not specify if the Height of the mode at frequency f=" << word[2] << " is fixed/free! ==> Using default condition (free parameter) " << std::endl;
                    //iMS_global.relax_H.push_back(1);
                    iMS_global.relax_gamma.push_back(1);
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
	iMS_global.freqs_ref.conservativeResize(iMS_global.relax_freq.size());
	iMS_global.els.conservativeResize(iMS_global.relax_freq.size());
	if(verbose == 1) {
		std::cout << " - List of relaxed parameters:" << std::endl;
		std::cout << "          l    nu       relax(nu)   relax(W)  relax(H)" << std::endl;
		for(int k=0; k<iMS_global.freqs_ref.size();k++){
			std::cout << "              " << iMS_global.els[k] << "       " << iMS_global.freqs_ref[k] << "       " << iMS_global.relax_freq[k] << "  "
				  << iMS_global.relax_gamma[k] << "      " << iMS_global.relax_H[k] << std::endl;
		}
	}
    
	// -------------------------------------
	
	i=0;
	cpt=0;
	iMS_global.hyper_priors.resize(10);
	if(verbose == 1) {std::cout << " - Hyper priors:" << std::endl;}
	while ((out < 4) && !cfg_session.eof()){ // the priors, until we reach the next # symbol
			std::getline(cfg_session, line0);
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				iMS_global.hyper_priors[i]=str_to_dbl(line0);
				cpt=cpt+1;
			} else{
				 out=out+1;
			}
			i=i+1;
	  }
	  iMS_global.hyper_priors.conservativeResize(cpt);
	  if(verbose == 1) {
		std::cout << iMS_global.hyper_priors.transpose() << std::endl;
	  }

	i=0;
	cpt=0;
	iMS_global.eigen_params.resize(200,6);
	std::getline(cfg_session, line0);
	if(verbose == 1) {std::cout << " - Initial guesses and frequency priors:" << std::endl;}
	while ((out < 5) && !cfg_session.eof() ){ //the priors, until we reach the 5th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				word=strsplit(line0, " \t");
				iMS_global.eigen_params.row(i)=str_to_Xdarr(line0, " \t");
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}
	iMS_global.eigen_params.conservativeResize(cpt, 6);
	if(verbose == 1) {
		std::cout << iMS_global.eigen_params << std::endl;
	}
	// -------------------------------------

	i=0;
	cpt=0;
	iMS_global.noise_params.resize(10);
	while ( (out < 6) && !cfg_session.eof() ){  // the priors, until we reach the 6th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){		
				word=strsplit(line0," \t");
          			for(int j=0; j<word.size(); j++){
          			 	iMS_global.noise_params[cpt+j]=str_to_dbl(word[j]);
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
        tmpXd=iMS_global.noise_params.segment(0,cpt);
        iMS_global.noise_params.setConstant(-1);
        iMS_global.noise_params.segment(10-cpt, cpt)=tmpXd;
      }
	if(verbose == 1) {
		std::cout << " - Noise inputs:" << std::endl;
		std::cout << iMS_global.noise_params.transpose() << std::endl;
	}
	// -------------------------------------	

	i=0;
	cpt=0;
	iMS_global.noise_s2.resize(10, 3);
	while ((out < 7) && !cfg_session.eof() ){ //the priors, until we reach the 7th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			if (char0 != "#"){
				word=strsplit(line0," \t");
				iMS_global.noise_s2.row(i)=str_to_Xdarr(line0, " \t");
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}
    if(i < 11){
        tmpXd=iMS_global.noise_s2.block(0,0,cpt,iMS_global.noise_s2.cols());
        iMS_global.noise_s2.setConstant(-1);
        iMS_global.noise_s2.block(10-cpt, 0, tmpXd.rows(), tmpXd.cols())=tmpXd;
    }
    if(verbose == 1) {
		std::cout << " - Noise information extracted during step s2:" << std::endl;
		std::cout << iMS_global.noise_s2 << std::endl;
    }
    //	exit(EXIT_SUCCESS);
	// -------------------------------------	
	
	// -------- The common parameters follow ------
	VectorXd a;
	i=0;
	cpt=0;
	iMS_global.modes_common.resize(50,5);
	iMS_global.modes_common.setConstant(-9999); // up to 10 variables and 4 prior parameters
	while ( (out < 9) && !cfg_session.eof() ){ // the initial values for the common parameters + priors, until we reach the 9th # symbol
			line0=strtrim(line0);
			char0=strtrim(line0.substr(0, 1));
			word=strsplit(line0," \t");
			if (char0 != "#"){
				word=strsplit(line0," \t");
				iMS_global.common_names.push_back(strtrim(word[0]));
				iMS_global.common_names_priors.push_back(strtrim(word[1]));
 				word.erase(word.begin()); // erase the slot containing the keyword name
				word.erase(word.begin()); // erase the slot containing the type of prior/switch
				a=arrstr_to_Xdarrdbl(word);
				
				for(int k=0; k<a.size();k++){
					iMS_global.modes_common(cpt, k)=a[k];
				}
				cpt=cpt+1;
		 	} else{
				out=out+1;
			}
			i=i+1;
			std::getline(cfg_session, line0);
	}

	iMS_global.modes_common.conservativeResize(iMS_global.common_names.size(), 5);

	if(verbose == 1) {
		std::cout << " - Common parameters for modes:" << std::endl;
		for(i=0; i<iMS_global.common_names.size();i++){
			std::cout << "  " << iMS_global.common_names[i] << "  ";
			std::cout << "  " << iMS_global.common_names_priors[i] << "  ";
			std::cout << "  " << iMS_global.modes_common.row(i) << std::endl;
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
    
   //exit(EXIT_SUCCESS);  
   return iMS_global;
}

Input_Data build_init_MS_Global(const MCMC_files inputs_MS_global, const bool verbose, const double resol){

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
	const std::vector<double> Vl{1, 1.5, 0.53, 0.08};
	double rho=pow(inputs_MS_global.Dnu/Dnu_sun,2.) * rho_sun;
	double Dnl=0.75, trunc_c=-1;
	double numax=inputs_MS_global.numax;
	

	// All Default booleans
	bool do_a11_eq_a12=1, do_avg_a1n=1, do_amp=0;
	bool bool_a1sini=0, bool_a1cosi=0;
	int lmax, en, ind, Ntot, p0, cpt, aj_switch, a2_param_count=0;
	uint8_t do_width_Appourchaux=0; // We need more than a boolean here, but no need to use a 64 bit signed int
	double tol=1e-2, tmp;
	VectorXi pos_el, pos_relax0, els_eigen, Nf_el(4), plength;
	VectorXd ratios_l, tmpXd, extra_priors, aj_param_count(6);
	std::vector<int> pos_relax;
	std::vector<double> f_inputs, h_inputs, w_inputs, f_priors_min, f_priors_max, f_el;;
	std::vector<bool> f_relax, h_relax, w_relax; 
	std::vector<int> rf_el, rw_el, rh_el;
	std::vector<std::string> tmpstr_vec;

	std::string tmpstr_h, tmpstr;
	
	Input_Data Snlm_in, Vis_in, Inc_in, Noise_in, freq_in, height_in, width_in; //, width_App2016_params; // This is by block, each category of parameters		
	Input_Data all_in; // The final structure of parameters, using the standards of my code
	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure

	// -------------- Set Extra_priors ----------------	
	extra_priors.resize(10);
	extra_priors[0]=1; // By default, we apply a smoothness condition
	extra_priors[1]=2.; // By default, the smoothness coeficient is 2 microHz
	extra_priors[2]=1e6; // aj/a1 No limit, j=1
	extra_priors[3]=0.50; // aj/a1 No limit, j=2
	extra_priors[4]=0.20; // aj/a1 No limit, j=3
	extra_priors[5]=0.15; // aj/a1 No limit, j=4
	extra_priors[6]=0.05; // aj/a1 No limit, j=5
	extra_priors[7]=0.05; // aj/a1 No limit, j=6
	extra_priors[8]=0; // Switch to control whether a prior imposes Sum(Hnlm)_{m=-l, m=+l}=1. Default: 0 (none). >0 values are model_dependent
	extra_priors[9]=-1;  // Specify rules to apply to specific models (e.g. models a1a2a3 must have extra_priors[5]=1)
	// ------------------------------------------------

	// Flatening and ordering of all the inputs/relax variables
	lmax=inputs_MS_global.els.maxCoeff();

	// -- Initialisation of structures --
    // --- Look for common instruction That must be run before the setup ---------
	all_in.model_fullname=" "; // Default is an empty string
	//all_in.prior_fullname="prior_MS_Global"; // Default set of prior
	aj_switch=0;
	aj_param_count.setZero();
    for(int i=0; i<inputs_MS_global.common_names.size(); i++){
        if(inputs_MS_global.common_names[i] == "model_fullname" ){ // This defines if we assume S11=S22 or not (the executed model will be different)
        	all_in.model_fullname=inputs_MS_global.common_names_priors[i];
            if(all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic" || all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic_v2" ||
               all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic_v3"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            }
            if(all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            	extra_priors[9]=0; // used in priors_calc.cpp 
            }
            if(all_in.model_fullname == "model_MS_Global_a1etaa3_Harvey1985"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;        
            	extra_priors[9]=0; // used in priors_calc.cpp     	
            }
            if(all_in.model_fullname == "model_MS_Global_a1a2a3_HarveyLike"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            	aj_switch=1; 
            	extra_priors[9]=1; // used in priors_calc.cpp           	
            }
/*            if(all_in.model_fullname == "model_MS_Global_a1acta3_HarveyLike"){ // OBSELETE
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            }
            if(all_in.model_fullname == "model_MS_Global_a1acta3_Harvey1985"){ // OBSELETE
         		do_a11_eq_a12=1;
           		do_avg_a1n=1;          	
           	}
*/
           	if(all_in.model_fullname == "model_MS_Global_a1l_etaa3_HarveyLike"){
           		//Previously corresponding to average_a1nl     bool    0    1 
           	    do_a11_eq_a12=0;
           		do_avg_a1n=1;
           		extra_priors[9]=2; // used in priors_calc.cpp 
           	}
        	if(all_in.model_fullname == "model_MS_Global_a1n_etaa3_HarveyLike"){
           		do_a11_eq_a12=1;
        		do_avg_a1n=0;
        		extra_priors[9]=3; // used in priors_calc.cpp 
            }
        	if(all_in.model_fullname == "model_MS_Global_a1nl_etaa3_HarveyLike"){
        		//Previously corresponding to average_a1nl     bool    0    0
        		do_a11_eq_a12=0;            		
        		do_avg_a1n=0;
        		extra_priors[9]=4; // used in priors_calc.cpp 
            }
        	if(all_in.model_fullname == "model_MS_Global_a1n_a2a3_HarveyLike"){
        		//Previously corresponding to average_a1nl     bool    0    0
        		do_a11_eq_a12=1;            		
        		do_avg_a1n=0;
        		aj_switch=2;
        		extra_priors[9]=5; // used in priors_calc.cpp 
            }
        	if(all_in.model_fullname == "model_MS_Global_a1l_a2a3_HarveyLike"){
        		//Previously corresponding to average_a1nl     bool    0    0
        		do_a11_eq_a12=0;            		
        		do_avg_a1n=1;
        		aj_switch=3;
        		extra_priors[9]=6; // used in priors_calc.cpp 
            }
        	if(all_in.model_fullname == "model_MS_Global_a1nl_a2a3_HarveyLike"){
        		//Previously corresponding to average_a1nl     bool    0    0
        		do_a11_eq_a12=0;            		
        		do_avg_a1n=0;
        		aj_switch=4;
        		extra_priors[9]=7; // used in priors_calc.cpp 
            }
        	if(all_in.model_fullname == "model_MS_Global_a1etaAlma3_HarveyLike"){
        		//Previously corresponding to average_a1nl     bool    0    0
        		do_a11_eq_a12=1;            		
        		do_avg_a1n=1;
        		aj_switch=5;
        		extra_priors[9]=8; // used in priors_calc.cpp 
            }
       	if(all_in.model_fullname == "model_MS_Global_aj_HarveyLike"){
        		//Previously corresponding to average_a1nl     bool    0    0
        		do_a11_eq_a12=1;            		
        		do_avg_a1n=1;
        		aj_switch=6; // Case for 2x(a1,a2,a3,a4,a5,a6) + 1 asymetry
        		extra_priors[9]=9; // used in priors_calc.cpp 
            }
            if(all_in.model_fullname == "model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1" || all_in.model_fullname == "model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2"){
            	do_a11_eq_a12=1;
            	do_avg_a1n=1;
            	if(all_in.model_fullname == "model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1"){ do_width_Appourchaux=1;}
            	if(all_in.model_fullname == "model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2"){ do_width_Appourchaux=2;}
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
	for(int el=0; el<lmax+1; el++){
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
	}
	
	// ------------------------------------------------------------------------------------------
	// ------------------------------- Handling the Common parameters ---------------------------
	// ------------------------------------------------------------------------------------------
	if( do_amp){
		std::cout << "   ===> Requested to fit squared amplitudes instead of Height... Converting height inputs into A_squared = pi*Height*Width..." << std::endl;
		tmpstr_h="Amplitude_l0";
		for(int i=0; i<h_inputs.size(); i++){
			h_inputs[i]=pi*w_inputs[i]*h_inputs[i]; 
		}
	} else{
		tmpstr_h="Height_l0";
	}
    // Set default value of priors for Height Width and frequency
	io_calls.initialise_param(&height_in, h_relax.size(), Nmax_prior_params, -1, -1);
	if (do_width_Appourchaux == 0){
		io_calls.initialise_param(&width_in, w_relax.size(), Nmax_prior_params, -1, -1); 
	} 
	if (do_width_Appourchaux == 1){
		io_calls.initialise_param(&width_in, 5, Nmax_prior_params, -1, -1);
	}
	if (do_width_Appourchaux == 2){
		io_calls.initialise_param(&width_in, 6, Nmax_prior_params, -1, -1);
	}
	if (do_width_Appourchaux != 0 && do_width_Appourchaux !=1 && do_width_Appourchaux !=2){
		std::cout << "Fatal error on io_ms_global.cpp: Allowed Appourchaux Widths models are only model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1 or" << std::endl;
		std::cout << "model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2. However the logical switches do_width_Appourchaux did not match any of those and are not consistent with non-Appourchaux width"<< std::endl;
		std::cout << "Serious debug required here. The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	io_calls.initialise_param(&freq_in, f_relax.size(), Nmax_prior_params, -1, -1); 


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

	// DEFAULT WIDTH
	tmpXd.resize(4);
	tmpXd << resol, inputs_MS_global.Dnu/3., -9999., -9999.;
	for(int i=0; i<w_inputs.size(); i++){
		if(w_relax[i]){
			io_calls.fill_param(&width_in, "Width_l0", "Jeffreys", h_inputs[i], tmpXd, i, 0);	
		} else{
			io_calls.fill_param(&width_in, "Width_l0", "Fix", h_inputs[i], tmpXd, i, 0);			
		}
	}

	
	// --- Default setup for frequencies ---
	for(int i=0; i<f_inputs.size(); i++){
		if(f_relax[i]){
			tmpXd << f_priors_min[i], f_priors_max[i], 0.01*inputs_MS_global.Dnu, 0.01*inputs_MS_global.Dnu; // default parameters for a GUG prior on frequencies
			io_calls.fill_param(&freq_in, "Frequency_l", "GUG", f_inputs[i], tmpXd, i, 0);	
		} else{
			io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f_inputs[i], tmpXd, i, 0);			
		}
	}
	// ----------- Calculate numax -----------
	// Flaten the output vector
	if(numax <=0){
		tmpXd.resize(Nf_el.sum());
		std::cout << "numax not provided. Input numax may be required by some models... Calculating numax..." << std::endl;
		//std::cout << height_in.inputs << std::endl;
		tmpXd.segment(0 , Nf_el[0])=height_in.inputs;
		cpt=Nf_el[0];
		std::cout << freq_in.inputs << std::endl;
		for(int el=1; el<=3; el++){
			for(int n=0; n<Nf_el[el]; n++){
				tmpXd[cpt+n]=std::abs(lin_interpol(freq_in.inputs.segment(0, Nf_el[0]), height_in.inputs, freq_in.inputs[Nf_el[0]+n]))*Vl[el]; // Modified on 1 Mar 2022 
				//std::cout << "l=" << el << "  n=" << n << std::endl;
			}
			cpt=cpt+Nf_el[el];
			//std::cout << "cpt[" << el <<	 "]=" << cpt << std::endl;
		}
		std::cout << "getting in getnumax..." << std::endl;
		numax=getnumax(freq_in.inputs , tmpXd); // We had to flatten the Height vector and put visibilities
		std::cout << "     numax: " << numax << std::endl;
	} else {
		std::cout << " Using provided numax: " << numax << std::endl;
	}
	std::cout << " ------------------" << std::endl;

	//exit(EXIT_SUCCESS);
 	// -------------------------------------

	// ----- Switch between the models that handle averaging over n,l or both -----
   if(do_a11_eq_a12 == 1 && do_avg_a1n == 1){
   		switch (aj_switch){
        	case 0:
        		io_calls.initialise_param(&Snlm_in, 6, Nmax_prior_params, -1, -1);
    			break;
    		case 1:
    			io_calls.initialise_param(&Snlm_in, 6 + 3, Nmax_prior_params, -1, -1); // Add 3 a2 coefficients in the a1a2a3 models
    			break;
    		case 2:
    			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax = Nf_el[0] = Nf_el[1] ... a2 coefficients in the a1n_a2a3 models
   				break;
    		case 3:
    			io_calls.initialise_param(&Snlm_in, 6 + lmax*3, Nmax_prior_params, -1, -1); // Add lmax*3 a2 coefficients in the a1l_a2a3 models
   				break;
 			case 4:
    			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax a2 coefficients in the a1nl_a2a3 models
   				break;
    		case 5:
    			io_calls.initialise_param(&Snlm_in, 6 + 3 + 2, Nmax_prior_params, -1, -1); // Add 3 epsilon_nl coefficients (a2_AR) and 2 theta parameters for the a1etaAlma3 models
    			break;
    		case 6:
    			io_calls.initialise_param(&Snlm_in, 12+1+1, Nmax_prior_params, -1, -1); // Add 2x(a1,a2,a3,a4,a5,a6) +  1 val (bool 0/1) for eta0 + 1 asymetry
    			tmpXd.resize(4);
    			tmpXd << -9999, -9999, -9999, -9999; 	
    			io_calls.fill_param(&Snlm_in, "eta0_switch", "Fix", 0, tmpXd, 12, 0); // set the defaut eta0 boolean to 0 (Not eta0 calculation)
    			break;		
   			default:
   				std::cout << "  Issues on the aj_switch configuration inside io_ms_global" << std::endl;
   				std::cout << "  aj_switch =" << aj_switch << " Has no preset rules" << std::endl;
   				std::cout << "  Debug required. The program will exit now" << std::endl;
   				exit(EXIT_SUCCESS);
    		}
    }  
    if(do_a11_eq_a12 == 0 && do_avg_a1n == 1){
   		switch (aj_switch){
        	case 0:
        		io_calls.initialise_param(&Snlm_in, 7, Nmax_prior_params, -1, -1);
    			break;
    		case 1:
    			io_calls.initialise_param(&Snlm_in, 7 + 3, Nmax_prior_params, -1, -1); // Add 3 a2 coefficients in the a1a2a3 models
    			break;
    		case 2:
    			io_calls.initialise_param(&Snlm_in, 7 + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax = Nf_el[0] = Nf_el[1] ... a2 coefficients in the a1n_a2a3 models
   				break;
    		case 3:
    			io_calls.initialise_param(&Snlm_in, 7 + lmax*3, Nmax_prior_params, -1, -1); // Add lmax*3 a2 coefficients in the a1l_a2a3 models
   				break;
 			case 4:
    			io_calls.initialise_param(&Snlm_in, 7 + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax a2 coefficients in the a1nl_a2a3 models
   				break;
   			default:
   				std::cout << "  Issues on the aj_switch configuration inside io_ms_global" << std::endl;
   				std::cout << "  aj_switch =" << aj_switch << " Has no preset rules in the case of do_a11_eq_a12 == 0 && do_avg_a1n == 1" << std::endl;
   				std::cout << "  Debug required. The program will exit now" << std::endl;
   				exit(EXIT_SUCCESS);
    		}

    }
    if(do_a11_eq_a12 == 1 && do_avg_a1n == 0){
		if (Nf_el[1] == Nf_el[2]){
   			switch (aj_switch){
        		case 0:
        			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1], Nmax_prior_params, -1, -1);
    				break;
    			case 1:
    				io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1] + 3, Nmax_prior_params, -1, -1); // Add 3 a2 coefficients in the a1a2a3 models
    				break;
    			case 2:
    				io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1] + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax = Nf_el[0] = Nf_el[1] ... a2 coefficients in the a1n_a2a3 models
   					break;
    			case 3:
    				io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1] + lmax*3, Nmax_prior_params, -1, -1); // Add lmax*3 a2 coefficients in the a1l_a2a3 models
   					break;
 				case 4:
    				io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1] + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax a2 coefficients in the a1nl_a2a3 models
   					break;
   				default:
   					std::cout << "  Issues on the aj_switch configuration inside io_ms_global" << std::endl;
   					std::cout << "  aj_switch =" << aj_switch << " Has no preset rules in the case of do_a11_eq_a12 == 1 && do_avg_a1n == 0" << std::endl;
   					std::cout << "  Debug required. The program will exit now" << std::endl;
   					exit(EXIT_SUCCESS);
    			}
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
   		switch (aj_switch){
       		case 0:
       			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2], Nmax_prior_params, -1, -1);
    			break;
    		case 1:
    			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2] + 3, Nmax_prior_params, -1, -1); // Add 3 a2 coefficients in the a1a2a3 models
    			break;
    		case 2:
    			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2] + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax = Nf_el[0] = Nf_el[1] ... a2 coefficients in the a1n_a2a3 models
   				break;
    		case 3:
    			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2] + lmax*3, Nmax_prior_params, -1, -1); // Add lmax*3 a2 coefficients in the a1l_a2a3 models
   				break;
 			case 4:
    			io_calls.initialise_param(&Snlm_in, 6 + Nf_el[1]+Nf_el[2] + Nf_el[0], Nmax_prior_params, -1, -1); // Add Nmax a2 coefficients in the a1nl_a2a3 models
   				break;
   			default:
   				std::cout << "  Issues on the aj_switch configuration inside io_ms_global" << std::endl;
   				std::cout << "  aj_switch =" << aj_switch << " Has no preset rules in the case of do_a11_eq_a12 == 0 && do_avg_a1n == 0" << std::endl;
   				std::cout << "  Debug required. The program will exit now" << std::endl;
   				exit(EXIT_SUCCESS);
    		} 
    }
 	
	for(int i=0; i<inputs_MS_global.common_names.size(); i++){	
		//std::cout << "5 - " << i << std::endl;
		//std::cout << "inputs_MS_global.common_names[i] = " << inputs_MS_global.common_names[i] << std::endl;

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
				try{
					trunc_c=str_to_dbl(inputs_MS_global.common_names_priors[i]); // Attempt to take the value of trunc_c from the 1 slot instead of the common part
				}
				catch(...){
					fatalerror_msg_io_MS_Global("trunc_c", "Fix", "[Truncation parameter]", "20" );
				}
			} else{
				trunc_c=inputs_MS_global.modes_common(i,0); // This value is added to the input vector at the end of this function.
			}
		}	

		// --- Frequencies ---
		if(inputs_MS_global.common_names[i] == "Frequency" || inputs_MS_global.common_names[i] == "frequency"){ 
			if(inputs_MS_global.common_names_priors[i] == "GUG" || inputs_MS_global.common_names_priors[i] == "Uniform"){
				for(int p0=0; p0<f_inputs.size(); p0++){
					if(f_relax[p0]){
						if(inputs_MS_global.common_names_priors[i] == "GUG"){
							tmpXd << f_priors_min[p0], f_priors_max[p0], inputs_MS_global.modes_common(i,3), inputs_MS_global.modes_common(i,4);
						} else{
							tmpXd << f_priors_min[p0], f_priors_max[p0], -9999, -9999; 						
						}
						io_calls.fill_param(&freq_in, "Frequency_l", inputs_MS_global.common_names_priors[i], f_inputs[p0], tmpXd, p0, 0);	
					} else{
						io_calls.fill_param(&freq_in, "Frequency_l", "Fix", f_inputs[p0], tmpXd, p0, 0);			
					}
				}
			} else{
				fatalerror_msg_io_MS_Global(inputs_MS_global.common_names[i], "GUG or Uniform", "", "" );
			}
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
				std::cout << "       Width argument found in the model file, but model name is: model_MS_Global_a1etaa3_AppWidth_HarveyLike. "  << std::endl;
				std::cout << "       This is incompatible. The argument will be ignored. Fit of the Width is performed using Appourchaux+2016 relation (5 parameters + numax)" << std::endl;
				std::cout << "       Pursuing..." <<std::endl;
		}	
		// --- Splittings and asymetry ---
		if(inputs_MS_global.common_names[i] == "splitting_a1" || inputs_MS_global.common_names[i] == "Splitting_a1"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("splitting_a1", "Fix_Auto", "", "" );
			}
			p0=0;
			io_calls.fill_param(&Snlm_in, "Splitting_a1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
        	if(do_a11_eq_a12 == 0 && do_avg_a1n == 1){ // Case a1(l), we put the same prior for a1(2) than a1(1)
		    	std::cout << "do_a11_eq_a12 == 0 && do_avg_a1n == 1" << std::endl;
		    	p0=6;
        	    io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_MS_global.modes_common.row(i), p0, 1);
        	}
        	if(do_a11_eq_a12 == 1 && do_avg_a1n == 0){ // Case a1(n), we put the same prior for a1(n)
				std::cout << "do_a11_eq_a12 == 1 && do_avg_a1n == 0" << std::endl;
				for(int kk=0; kk<Nf_el[1]; kk++){
					p0=6+kk;
					io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_MS_global.modes_common.row(i), p0, 1); 
				}
				p0=0;
				io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, inputs_MS_global.modes_common.row(i), p0, 1);
        	}
        	if(do_a11_eq_a12 == 0 && do_avg_a1n == 0){ // Case a1(n,l), we put the same prior for a1(n,l)
				std::cout << "do_a11_eq_a12 == 0 && do_avg_a1n == 0" << std::endl;
				for(int kk=0; kk<Nf_el[1]+Nf_el[2]; kk++){
					p0=6+kk;
					io_calls.fill_param(&Snlm_in, Snlm_in.inputs_names[0], Snlm_in.priors_names[0], Snlm_in.inputs[0], inputs_MS_global.modes_common.row(i), p0, 1); 
				}
				p0=0;
				io_calls.fill_param(&Snlm_in, "Empty", "Fix", 0, inputs_MS_global.modes_common.row(i), p0, 1); // Erase values of the default Splitting_a1 block
        	}
		}
		if(inputs_MS_global.common_names[i] == "asphericity_eta"|| inputs_MS_global.common_names[i] == "Asphericity_eta"){
			std::cout << " Warning: Update on the code (07/12/2021) does not allow to use Asphericity_eta as a keyword in io_ms_global models" << std::endl;
			std::cout << "          The calculation of eta is now made systematically within the models and no longer requires the eta parameter" << std::endl;
			std::cout << "          Consequently, this parameter will be ignored" << std::endl;
			Snlm_in.inputs_names[1]="Empty";
			Snlm_in.priors_names[1]="Fix";
			Snlm_in.inputs_names[1]="Asphericity_eta";
			Snlm_in.relax[1]=0;
			Snlm_in.inputs[1]=0;
			/*
			if (aj_switch == 0){
				if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){ // Case where the centrifugal force is added but FIXED.
					Snlm_in.priors_names[1]="Fix"; //In all cases the centrifugal force is fixed
					Snlm_in.relax[1]=0;
					if(inputs_MS_global.modes_common(i,0) == 1){
  						if (Snlm_in.inputs[0] != -9999){ 
							//eta=(4./3.)*!pi*Dnl*(a1_init*1d-6)^2/(rho*G)
							Snlm_in.inputs[1]=(4./3.)*pi*Dnl*pow(Snlm_in.inputs[0]*1e-6,2.)/(rho*G);
							std::cout << " -------------" << std::endl;
							std::cout << "Asphericity_eta given with Fix_Auto ==> Setting the asphericity to the centrifugal force" << std::endl;
							std::cout << "      eta=" << Snlm_in.inputs[1] << std::endl;
							std::cout << " -------------" << std::endl;
						} else{
							std::cout << "Warning: the keyword 'asphericity_eta' must appear after the keyword splitting_a1" << std::endl;
							std::cout << "         This because the splitting_a1 is used to define the initial value of asphericity" << std::endl;
							std::cout << "         Edit the .MCMC file accordingly" << std::endl;
							std::cout << "The program will exit now" << std::endl;
							exit(EXIT_FAILURE);
						}
					} else{
						Snlm_in.inputs[1]=0;
					}
				} else{ // Case where the centrifugal force is user-defined: Could be free or fixed
					p0=1;
					io_calls.fill_param(&Snlm_in, "Asphericity_eta", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
				}
			} else{
				std::cout << " asphericity_eta global parameter found in .model file but the chosen model is incompatible with this (aj_switch > 0). " << std::endl;
				std::cout << " If you use a model that fits a2 coefficients, please use the 'a2' parameter and set to Auto or Fix.  Or for models a1l_a2a3 or a1a2a3 use " << std::endl;
				std::cout << " the 'a2_0', 'a2_1' and 'a2_2' parameters... a2_0:constant term, a2_1: Linear term, a2_3: Quadratic term" << std::endl;
				exit(EXIT_FAILURE);
			}
			*/
		}
		if(inputs_MS_global.common_names[i] == "a2" && aj_switch >=2){ // a2 is a valid keyword only for a1n_a2a3, a1l_a2a3, a1nl_a2a3 models
			std::cout << "io_ms_global.cpp : NEED TO BE UPDATED TO HANDLES MODELS WITH n or l dependence of a2" << std::endl;
			std::cout << "The program will exit now" << std::endl;
			exit(EXIT_SUCCESS);
		}		
		if(inputs_MS_global.common_names[i] == "a2" && aj_switch <2){ // a2 is a valid keyword only for a1n_a2a3, a1l_a2a3, a1nl_a2a3 models
			std::cout << "   a2 keyword is not valid for fitting a2 coefficients in case of a1n_a2a3, a1l_a2a3, a1nl_a2a3 models " << std::endl;
			std::cout << "   models a1l_a2a3 or a1a2a3 use the 'a2_0', 'a2_1' and 'a2_2' parameters... a2_0:constant term, a2_1: Linear term, a2_3: Quadratic term" << std::endl;
			std::cout << "   Please change the .model accordingly" << std::endl;
			std::cout << "   The program will exit now" << std::endl;
			exit(EXIT_SUCCESS);
		}
		if(inputs_MS_global.common_names[i] == "a2_0" && aj_switch == 1){ // a2 is a valid keyword only for a1n_a2a3, a1l_a2a3, a1nl_a2a3 models
			a2_param_count=a2_param_count+1;
			p0=6; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "a2_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a2_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a2_1"  && aj_switch == 1){ // a2 is a valid keyword only for a1n_a2a3, a1l_a2a3, a1nl_a2a3 models
			a2_param_count=a2_param_count+1;
			p0=6+1; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "a2_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a2_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a2_2"  && aj_switch == 1){ // a2 is a valid keyword only for a1n_a2a3, a1l_a2a3, a1nl_a2a3 models
			a2_param_count=a2_param_count+1;
			p0=6+2; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "a2_2", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a2_2 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}		
		}
		aj_param_count=aj_param_count + settings_aj_splittings(i,inputs_MS_global, &Snlm_in, aj_switch);

//  Dealing with models using epsilon_nl see Gizon 2002, AN
		if(inputs_MS_global.common_names[i] == "epsilon_0" && aj_switch == 5){ // aj_switch=5 is for model_MS_Global_a1etaAlma3_HarveyLike
			a2_param_count=a2_param_count+1;
			p0=6; 
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "epsilon_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for epsilon_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "epsilon_1"  && aj_switch == 5){ // aj_switch=5 is for model_MS_Global_a1etaAlma3_HarveyLike
			a2_param_count=a2_param_count+1;
			p0=6+1; 
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "epsilon_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for epsilon_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "epsilon_2"  && aj_switch == 5){ // aj_switch=5 is for model_MS_Global_a1etaAlma3_HarveyLike
			a2_param_count=a2_param_count+1;
			p0=6+2; 
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "espilon_2", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for epsilon_2 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}		
		}
		if(inputs_MS_global.common_names[i] == "theta0"  && aj_switch == 5){ //theta0: central latitude for the Active Region. aj_switch=5 is for model_MS_Global_a1etaAlma3_HarveyLike
			a2_param_count=a2_param_count+1;
			p0=6+3; 
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "theta0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for epsilon_2 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}		
		}
		if(inputs_MS_global.common_names[i] == "Dtheta"  && aj_switch == 5){ //theta_min and theta_max = theta0 +/- Dtheta/2: Max angle for the Active Region. aj_switch=5 is for model_MS_Global_a1etaAlma3_HarveyLike
			a2_param_count=a2_param_count+1;
			p0=6+4; 
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(&Snlm_in, "Dtheta", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for epsilon_2 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}		
		}
// ---
		if(inputs_MS_global.common_names[i] == "splitting_a3" || inputs_MS_global.common_names[i] == "Splitting_a3"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("splitting_a3", "Fix_Auto", "", "" );
			}
			p0=2;
			io_calls.fill_param(&Snlm_in, "Splitting_a3", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
		}

		if(inputs_MS_global.common_names[i] == "asymetry" || inputs_MS_global.common_names[i] == "Asymetry"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("asymetry", "Fix_Auto", "", "" );
			}
			if (aj_switch != 6){
				p0=5;
			} else{
				p0=Snlm_in.inputs.size()-1; // The new standard for aj fitting is that the Asymetry is at the end of the Snlm 
			}
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
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("sqrt(splitting_a1).cosi", "Fix_Auto", "", "" );
			}
			p0=3;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).cosi", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);

            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0  and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1cosi=1;
		}
		
		if(inputs_MS_global.common_names[i] == "sqrt(splitting_a1).sini"){
			if(inputs_MS_global.common_names_priors[i] == "Fix_Auto"){
				fatalerror_msg_io_MS_Global("sqrt(splitting_a1).sini", "Fix_Auto", "", "" );
			}
			p0=4;
			io_calls.fill_param(&Snlm_in, "sqrt(splitting_a1).sini", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);
			
            if(do_a11_eq_a12 == 0 || do_avg_a1n == 0){ // In that case, we put the same prior for a1(2) than a1(1)
				std::cout << "Warning: do_a11_eq_a12=0 and/or do_avg_a1n=0 (models *_a1n* or *_a1l*) was requested but is not available when fitting sqrt(a1).cosi and sqrt(a1).sini" << std::endl;
				std::cout << "         You must modify the code accordingly if you want to implement a11 and a12 in that scenario" << std::endl;
				std::cout << "         The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
            }
            bool_a1sini=1;
		}
	}

if (aj_switch == 1 && a2_param_count !=3){
	std::cout << " Invalid number of constraints for a2: Please set a2_0, a2_1 and a2_2" <<std::endl;
	exit(EXIT_SUCCESS); 
}
if (aj_switch == 5 && a2_param_count !=5){
	std::cout << " Invalid number of constraints for epsilon: Please set epsilon_0, epsilon_1, epsilon_2 and theta_min, Dtheta for the considered model" <<std::endl;
	exit(EXIT_SUCCESS); 
}
if (aj_switch == 6 && aj_param_count.sum() !=12){
	std::cout << "aj_param_count = " << aj_param_count << std::endl;
	std::cout << " Invalid number of constraints for epsilon: Please set a1_0, a1_1, a2_0, a2_1, a3_0, a3_1, a4_0, a4_1, a5_0, a5_1, a6_0, a6_1 (12 parameters) for the considered model" <<std::endl;
	exit(EXIT_SUCCESS); 	
}
if(bool_a1cosi != bool_a1sini){ // Case when one of the projected splitting quantities is missing ==> Problem
	std::cout << "Warning: Both 'sqrt(splitting_a1).sini' and 'sqrt(splitting_a1).cosi' keywords must appear" << std::endl;
	std::cout << "         It is forbidden to use only one of them" << std::endl;
	std::cout << "         Edit the .MODEL file accordingly" << std::endl;
	std::cout << "The program will exit now" << std::endl;
	exit(EXIT_FAILURE);
}
if((bool_a1cosi == 0) && (bool_a1sini == 0)){ // Case where Inclination and Splitting_a1 are supposed to be used (CLASSIC and CLASSIC_vX models)
	if( (all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike") || (all_in.model_fullname == "model_MS_Global_a1etaa3_Harvey1985") ||
		(all_in.model_fullname == "model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1") || (all_in.model_fullname=="model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2")  ||
		(all_in.model_fullname == "model_MS_Global_a1n_a2a3_HarveyLike") || (all_in.model_fullname == " model_MS_Global_a1nl_a2a3_HarveyLike")  || (all_in.model_fullname == "model_MS_Global_a1a2a3_HarveyLike") ||
		(all_in.model_fullname == "model_MS_Global_a1etaAlma3_HarveyLike") ){
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
	// ----------------                                                                                      ----------------
	// ----------------                                                                                      ----------------
	// The models in this section are not fitting the inclination. Instead they consider Hlm or the Ratios as free parameters
	// ----------------                                                                                      ----------------
	// ----------------    
	//----------------

	if(all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic_v2"){
		tmp=Inc_in.inputs[0]; // Recover the initial stellar inclination as defined by the user in the .model file
		io_calls.initialise_param(&Inc_in, 9, Nmax_prior_params, -1, -1); // 2 param for l=1, 3 for l=2 and 4 for l=3
		ind=0;
		tmpXd.resize(4);
		tmpXd << 0, 1, -9999., -9999.;
		for(int el=1; el<=lmax; el++){
			ratios_l=amplitude_ratio(el, tmp); 
			for(int em=0; em<=el; em++){
				io_calls.fill_param(&Inc_in, "Inc:H" + int_to_str(el) + "," + int_to_str(em), "Uniform", ratios_l[el+em], tmpXd, ind, 0);
				ind=ind+1;
			}
		}
		extra_priors[8]=1; // Impose Sum(H(nlm))_{l=-m,l=+m} =1 for that model (case == 1)
	}
	if(all_in.model_fullname == "model_MS_Global_a1etaa3_HarveyLike_Classic_v3"){
		for(int i=0; i<8;i++){std::cout << "  ------------------------------------------------------------------------------" <<std::endl;}
		std::cout << "WARNING      WARNING     WARNING     WARNING     WARNING     WARNING     WARNING     " << std::endl;
		std::cout << "                   THIS FUNCTION MAY NOT BE SUITABLE FOR A GLOBAL FIT! " <<std::endl;
		std::cout << "            YOUR ARE LIKELY TO HAVE TOO MANY PARAMETERS FOR MODES HEIGHTS " << std::endl;
		std::cout << "      WE RECOMMEND TO USE model_MS_Global_a1etaa3_HarveyLike_Classic_v2 FOR A GLOBAL FIT" << std::endl;
		std::cout << "               OR TO SWITCH TO A LOCAL FIT (see config_default.cfg and model model_MS_local_basic_v2)" <<std::endl;
		std::cout << "WARNING      WARNING     WARNING     WARNING     WARNING     WARNING     WARNING     " << std::endl;
		for(int i=0; i<8;i++){std::cout << "  ------------------------------------------------------------------------------" <<std::endl;}
	
		tmpXd.resize(4);
		tmpXd << -9999, -9999, -9999., -9999.;
		tmp=Inc_in.inputs[0]; // Recover the initial stellar inclination as defined by the user in the .model file
		// Visibilities are not used in the is model ==> deactivated
		const VectorXd Vistmp=Vis_in.inputs; 
		std::cout << Vistmp << std::endl;
		for (int el=1; el<lmax;el++){
			io_calls.fill_param(&Vis_in, "Empty", "Fix", 0, tmpXd, el-1, 0); 
		}
		//
		io_calls.initialise_param(&Inc_in, Nf_el[1]*2 + Nf_el[2]*3 + Nf_el[3]*4, Nmax_prior_params, -1, -1); // 2 param for each l=1, 3 for l=2 and 4 for l=3
		ind=0;
		tmpXd << Hmin, Hmax, -9999., -9999.;
		for(int el=1; el<=lmax; el++){
			ratios_l=amplitude_ratio(el, tmp); 
			for(int en=0; en<Nf_el[el]; en++){
				for(int em=0; em<=el; em++){
					io_calls.fill_param(&Inc_in, "Inc: H" + int_to_str(en) + "," + int_to_str(el) + "," + int_to_str(em), "Jeffreys", height_in.inputs[en]*Vistmp[el-1]*ratios_l[el+em], tmpXd, ind, 0);
					ind=ind+1;
				}
			}
		}
		extra_priors[8]=2; // Cannot impose Sum(H(nlm))_{l=-m,l=+m} =1 for that model (case == 2) because we do not rely on visbilities
	}
	// ----------------                                                                                      ----------------
	// ----------------                                                                                      ----------------
	// ----------------                                                                                      ----------------
	// ----------------                                                                                      ----------------

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
	plength[6]=Snlm_in.inputs.size(); plength[7]=w_inputs.size() ; plength[8]=Noise_in.inputs.size(); 
	plength[9]=Inc_in.inputs.size();
	plength[10]=2; // This is trunc_c and do_amp;
	
	io_calls.initialise_param(&all_in, plength.sum(), Nmax_prior_params, plength, extra_priors);
	
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
	if(do_width_Appourchaux == 0){ // Most models
		io_calls.add_param(&all_in, &width_in, p0);
	} 
	if(do_width_Appourchaux == 1){// Case of: model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1
		width_in=set_width_App2016_params_v1(numax, width_in);
		io_calls.add_param(&all_in, &width_in, p0);
	}
	if(do_width_Appourchaux == 2){// Case of: model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1
		width_in=set_width_App2016_params_v2(numax, width_in);
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
		std::cout << "Warning: trunc_c = " << all_in.inputs[p0] << " <= 0. This is forbidden. Setting default to 10000. (No truncation)" << std::endl;
		all_in.inputs[p0]=10000.; // In case of a non-sense value for c, we use Full-Lorentzian as default
	}
	// -- Add the Amplitude switch --
	p0=all_in.plength[0] + all_in.plength[1] + all_in.plength[2] + all_in.plength[3] + all_in.plength[4] + all_in.plength[5] + all_in.plength[6] + all_in.plength[7] + all_in.plength[8] + all_in.plength[9] + 1;
	io_calls.fill_param(&all_in, "Switch for fit of Amplitudes or Heights", "Fix", do_amp, inputs_MS_global.modes_common.row(0), p0,1);
			
	if(verbose == 1){
		std::cout << " ----------------- Configuration summary -------------------" << std::endl;
		std::cout << "Model Name = " << all_in.model_fullname << std::endl;
		if(all_in.extra_priors[0] == 0){ 
			std::cout << "   freq_smoothness is set to 0 ==> NO smoothness condition on frequencies" << std::endl;
		} else{
			std::cout << "   freq_smoothness is set to 1 ==> APPLIES a smoothness condition on frequencies" << std::endl;
			std::cout << "   smoothness coeficient as specified by the user (of defined by default): " << all_in.extra_priors[1] << " microHz" << std::endl;
		}

		for(int eli=0; eli<6; eli++){ //   aj, with 1<j<6
			std::cout << "    Maximum ratio between a"<<eli+1<<" and a1: " << all_in.extra_priors[2+eli] << std::endl;
		}
			
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
	
	if (aj_switch >=2 && aj_switch != 5 && aj_switch !=6){ // aj_switch = 5 is a model that use eta and Alm to get a2_CF and a2_AR and is already set
		std::cout << " io_ms_global.cpp : aj_switch >=2 not yet properly configured yet" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_SUCCESS);
	}
	//std::cout << "Exiting test " << std::endl;
	//exit(EXIT_SUCCESS);

return all_in;
}


short int set_noise_params(Input_Data *Noise_in, const MatrixXd& noise_s2, const VectorXd& noise_params){
/*
 *
 * A function that prepares the Noise_in Data structure using 
 * inputs of the .model file regrouped inside the noise_s2 matrixXd
 *
*/

	(*Noise_in).inputs_names[0]="Harvey-Noise_H"; (*Noise_in).inputs_names[3]="Harvey-Noise_H"; (*Noise_in).inputs_names[6]="Harvey-Noise_H";
	(*Noise_in).inputs_names[1]="Harvey-Noise_tc"; (*Noise_in).inputs_names[4]="Harvey-Noise_tc"; (*Noise_in).inputs_names[7]="Harvey-Noise_tc";
	(*Noise_in).inputs_names[2]="Harvey-Noise_p"; (*Noise_in).inputs_names[5]="Harvey-Noise_p"; (*Noise_in).inputs_names[8]="Harvey-Noise_p";
	(*Noise_in).inputs_names[9]="White_Noise_N0";

	(*Noise_in).priors_names[0]="Fix"; (*Noise_in).priors_names[3]="Fix"; (*Noise_in).priors_names[6]="Gaussian";
	(*Noise_in).priors_names[1]="Fix"; (*Noise_in).priors_names[4]="Fix"; (*Noise_in).priors_names[7]="Gaussian";
	(*Noise_in).priors_names[2]="Fix"; (*Noise_in).priors_names[5]="Fix"; (*Noise_in).priors_names[8]="Gaussian";
	(*Noise_in).priors_names[9]="Gaussian";
	//(*Noise_in).priors_names[9]="Fix";

	(*Noise_in).relax[0]=0; (*Noise_in).relax[3]=0; (*Noise_in).relax[6]=1;
	(*Noise_in).relax[1]=0; (*Noise_in).relax[4]=0; (*Noise_in).relax[7]=1;
	(*Noise_in).relax[2]=0; (*Noise_in).relax[5]=0; (*Noise_in).relax[8]=1;
	(*Noise_in).relax[9]=1; // 1 WHITE NOISE IS FIXED !
	(*Noise_in).inputs=noise_params;

	// Handle cases with negative H or tc ==> Harvey is Fix to 0 (no Harvey) <==> case of simulations with white noise
	if(((*Noise_in).inputs[0] <= 0) || ((*Noise_in).inputs[1] <= 0) || ((*Noise_in).inputs[2] <= 0)){
		(*Noise_in).priors_names[0]="Fix"; (*Noise_in).priors_names[1]="Fix"; (*Noise_in).priors_names[2]="Fix";
		(*Noise_in).relax[0]=0;            (*Noise_in).relax[1]=0;            (*Noise_in).relax[2]=0;
		(*Noise_in).inputs[0]=0;		    (*Noise_in).inputs[1]=0;		    (*Noise_in).inputs[2]=1;
	}
	if(((*Noise_in).inputs[3] <= 0) || ((*Noise_in).inputs[4] <= 0) || ((*Noise_in).inputs[5] <= 0)){
		(*Noise_in).priors_names[3]="Fix"; (*Noise_in).priors_names[4]="Fix"; (*Noise_in).priors_names[5]="Fix";
		(*Noise_in).relax[3]=0;            (*Noise_in).relax[4]=0;            (*Noise_in).relax[5]=0;
		(*Noise_in).inputs[3]=0;		    (*Noise_in).inputs[4]=0;		    (*Noise_in).inputs[5]=1;
	}
	if(((*Noise_in).inputs[6] <= 0) || ((*Noise_in).inputs[7] <= 0) || ((*Noise_in).inputs[8] <= 0)){
		(*Noise_in).priors_names[6]="Fix"; (*Noise_in).priors_names[7]="Fix"; (*Noise_in).priors_names[8]="Fix";
		(*Noise_in).relax[6]=0;            (*Noise_in).relax[7]=0;            (*Noise_in).relax[8]=0;
		(*Noise_in).inputs[6]=0;		    (*Noise_in).inputs[7]=0;		    (*Noise_in).inputs[8]=1;
	}
	// --- Center of the Gaussian -----
	(*Noise_in).priors(0,6)=noise_s2(6,0); 
	(*Noise_in).priors(0,7)=noise_s2(7,0);
	(*Noise_in).priors(0,8)=noise_s2(8,0);
	(*Noise_in).priors(0,9)=noise_s2(9,0);
	// --- sigma of the Gaussian
	(*Noise_in).priors(1,6)=(noise_s2(6,1) + noise_s2(6,2))*3./2;
	(*Noise_in).priors(1,7)=(noise_s2(7,1) + noise_s2(7,2))*3./2;
	if(noise_s2(8,1) !=0){ // If p has given errors then set uncertainty to 3*error
		(*Noise_in).priors(1,8)=(noise_s2(8,1) + noise_s2(8,2))*3./2;
	} else{ // If p is given with null-errors then set uncertainty to 0.1*p
		(*Noise_in).priors(1,8)=(*Noise_in).priors(0,8)*0.1;
	}
	//(*Noise_in).priors(1,9)=(noise_s2(9,1) + noise_s2(9,2))*10./2;
    if((*Noise_in).priors_names[9] == "Uniform"){
    	(*Noise_in).priors(1,9)=(noise_s2(9,1) + noise_s2(9,2)); // This is valid only if the noise prior is Uniform
	} else{
		(*Noise_in).priors(1,9)=(*Noise_in).priors(0,9)*0.1; // This is valid only i the noise prior is Gaussian
	}
	if((((*Noise_in).priors(1,6)/(*Noise_in).priors(0,6)) <= 0.05) && (*Noise_in).priors_names[6] != "Fix"){ // If the given relative uncertainty on H3 is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the Height of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 5% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 5%" << std::endl;
		(*Noise_in).priors(1,6)=(*Noise_in).priors(0,6)*0.05;
		std::cout << "	       Resuming..." << std::endl;
	}
	if((((*Noise_in).priors(1,7)/(*Noise_in).priors(0,7)) <= 0.005) && (*Noise_in).priors_names[7] != "Fix"){ // If the given relative uncertainty on tc3 is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the timescale of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 0.5% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 0.5%" << std::endl;
		(*Noise_in).priors(1,7)=(*Noise_in).priors(0,7)*0.005;
		std::cout << "	       Resuming..." << std::endl;
	}
	if((((*Noise_in).priors(1,8)/(*Noise_in).priors(0,8)) <= 0.05) && (*Noise_in).priors_names[8] != "Fix"){ // If the given relative uncertainty on p is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the power of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 5% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 5%" << std::endl;
		(*Noise_in).priors(1,8)=(*Noise_in).priors(0,8)*0.05;
		std::cout << "	       Resuming..." << std::endl;
	}
	if((((*Noise_in).priors(1,9)/(*Noise_in).priors(0,9)) <= 0.0005) && (*Noise_in).priors_names[9] != "Fix"){ // If the given relative uncertainty on N0 is smaller than 5%
		std::cout << "Warning: The relative uncertainty on the Height of the high-frequency Harvey profile" << std::endl;
		std::cout << "         is smaller than 0.05% in the .MCMC file. This is too small" << std::endl;
		std::cout << "         ==> Relative uncertainty forced to be of 0.05%" << std::endl;
		(*Noise_in).priors(1,9)=(*Noise_in).priors(0,9)*0.0005;
		std::cout << "	       Resuming..." << std::endl;
	}
	return 0;
}

short int fatalerror_msg_io_MS_Global(const std::string varname, const std::string param_type, const std::string syntax_vals, const std::string example_vals){
/*
* Function that handle error messages and warnings for io_MS_Global
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

double getnumax(const VectorXd& fl, const VectorXd& Hl){
/*
* Function that uses Heights and frequencies of modes in order to calculate numax
* The vector of inputs must be flat
*/
	double numax=0;
	double Htot=0;
	
	for(long i=0; i<fl.size();i++){
	 	numax=numax + fl[i]*Hl[i];
	}
	Htot=Hl.sum();
	numax=numax/Htot;
	
	return numax;
}

Input_Data set_width_App2016_params_v1(const double numax, Input_Data width_in){
/* 
 * Function that calculates the initial guesses for the widths using numax and the linear fit reported in Appourchaux+2016
 * Note that these values are taken by hand from the graphs
 *
 * Note It would be better to have a function that takes the input widths, fit them using Appourchaux relation so that we 
 * get good initial guess, on a case-by-case basis. However, this would require to implement more dependencies in the code
 * (gradient descent minimisation library/algorithms). This is not plan at the moment 
*/
 
 	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure

 	VectorXd out(5); // numax, nudip, alpha, Gamma_alfa, Wdip, DeltaGammadip
 	MatrixXd priors(5,4);
 
 	priors.setConstant(-9999); // Set the default value for priors
 	
 	// Input values
 	out[0]=numax; // nudip
 	out[1]=4./2150.*numax + (1. - 1000.*4./2150.); // alpha
 	out[2]=0.8/2150.*numax + (4.5 - 1000.*0.8/2150.); // Gamma_alpha. Linear for a1.nu + a0... using graphical reading of App2016
 	out[3]=3400./2150.*numax + (1000. - 1000.*3400./2150.); // Wdip
 	out[4]=2.8/2200.*numax + (1. - 2.8/2200. * 1.); //DeltaGammadip
 		
 	// Priors on the parameters... most of those are put completely wildely: Would need to plot the graphs from App2016 to put proper gaussians
 	priors(0,0)=out[0];
 	priors(0,1)=out[0]*0.1; // 10% of numax on nudip
 	priors(1,0)=out[1];
 	priors(1,1)=out[1]*0.2; // 20% of alpha
 	priors(2,0)=out[2];
 	priors(2,1)=out[2]*0.2; // 20% of Gamma_alpha
 	priors(3,0)=out[3];
 	priors(3,1)=out[3]*0.2; // 20% of Wdip
 	priors(4,0)=out[4];
 	priors(4,1)=out[4]*0.4; // 20% of DeltaGammadip

	//io_calls.show_param(width_in, 0);
	
	// Filling the structure of width parameters
	io_calls.fill_param(&width_in, "width:Appourchaux_v1:nudip", "Gaussian", out[0],priors.row(0), 0, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v1:alpha", "Gaussian", out[1],priors.row(1), 1, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v1:Gamma_alpha", "Gaussian", out[2],priors.row(2), 2, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v1:Wdip", "Gaussian", out[3],priors.row(3), 3, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v1:DeltaGammadip", "Gaussian", out[4],priors.row(4), 4, 0);
	
	//io_calls.show_param(width_in, 0);
	return width_in;
}

Input_Data set_width_App2016_params_v2(const double numax, Input_Data width_in){
/* 
 * Function that calculates the initial guesses for the widths using numax and the linear fit reported in Appourchaux+2016
 * Note that these values are taken by hand from the graphs
 *
 * Note It would be better to have a function that takes the input widths, fit them using Appourchaux relation so that we 
 * get good initial guess, on a case-by-case basis. However, this would require to implement more dependencies in the code
 * (gradient descent minimisation library/algorithms). This is not plan at the moment 
*/
 
 	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure

 	VectorXd out(6); // numax, nudip, alpha, Gamma_alfa, Wdip, DeltaGammadip
 	MatrixXd priors(6,4);
 
 	priors.setConstant(-9999); // Set the default value for priors
 	
 	// Input values
 	out[0]=std::abs(numax); // numax
 	out[1]=std::abs(numax); // nudip
 	out[2]=std::abs(4./2150.*numax + (1. - 1000.*4./2150.)); // alpha
 	out[3]=std::abs(0.8/2150.*numax + (4.5 - 1000.*0.8/2150.)); // Gamma_alpha. Linear for a1.nu + a0... using graphical reading of App2016
 	out[4]=std::abs(3400./2150.*numax + (1000. - 1000.*3400./2150.)); // Wdip
 	out[5]=std::abs(2.8/2200.*numax + (1. - 2.8/2200. * 1.)); //DeltaGammadip

	if(numax < 800){
		out[3]=out[3]/5;
 	}	
 	// Priors on the parameters... most of those are put completely wildely: Would need to plot the graphs from App2016 to put proper gaussians
 	/*  OLD PRIORS (CHANGED ON 01/12/2020)
 	priors(0,0)=out[0];
 	priors(0,1)=out[0]*0.1; // 10% of numax on numax
 	priors(1,0)=out[1];
 	priors(1,1)=out[1]*0.1; // 10% of numax on nudip
 	priors(2,0)=out[2];
 	priors(2,1)=out[2]*0.2; // 20% of alpha
 	priors(3,0)=out[3];
 	priors(3,1)=out[3]*0.2; // 20% of Gamma_alpha
 	priors(4,0)=out[4];
 	priors(4,1)=out[4]*0.2; // 20% of Wdip
 	priors(5,0)=out[5];
 	priors(5,1)=out[5]*0.4; // 20% of DeltaGammadip
	
	//io_calls.show_param(width_in, 0);
	
	// Filling the structure of width parameters
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:numax", "Gaussian", out[0],priors.row(0), 0, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:nudip", "Gaussian", out[1],priors.row(1), 1, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:alpha", "Gaussian", out[2],priors.row(2), 2, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:Gamma_alpha", "Uniform", out[3],priors.row(3), 3, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:Wdip", "Gaussian", out[4],priors.row(4), 4, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:DeltaGammadip", "Gaussian", out[5],priors.row(5), 5, 0);
	*/

	// Priors on the parameters... most of those are put completely wildely: Would need to plot the graphs from App2016 to put proper gaussians
	// POSITION PARAMETERS WITH GAUSSIAN PRIORS AND INTENSIVE PARAMETERS AS UNIFORM
 	
    priors(0,0)=out[0];
 	priors(0,1)=out[0]*0.1; // 10% of numax on numax
 	priors(1,0)=out[1];
 	priors(1,1)=out[1]*0.1; // 10% of numax on nudip
 	priors(2,0)=0;
 	priors(2,1)=6; // min/max of the plot in Appourchaux 2016
 	priors(3,0)=0; //
 	priors(3,1)=10; // min/max of the plot in Appourchaux 2016
 	priors(4,0)=out[4];
 	priors(4,1)=out[4]*0.25; // 20% of Wdip
 	priors(5,0)=0.;
 	priors(5,1)=15; // min/max of the plot in Appourchaux 2016
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:numax", "Gaussian", out[0],priors.row(0), 0, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:nudip", "Gaussian", out[1],priors.row(1), 1, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:alpha", "Uniform", out[2],priors.row(2), 2, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:Gamma_alpha", "Uniform", out[3],priors.row(3), 3, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:Wdip", "Gaussian", out[4],priors.row(4), 4, 0);
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:DeltaGammadip", "Uniform", out[5],priors.row(5), 5, 0);
	

	// --- TEST BY FIXING ---
        /*
	io_calls.fill_param(&width_in, "width:Appourchaux_v2:numax", "Fix", out[0],priors.row(0), 0, 0);
        io_calls.fill_param(&width_in, "width:Appourchaux_v2:nudip", "Fix", out[1],priors.row(1), 1, 0);
        io_calls.fill_param(&width_in, "width:Appourchaux_v2:alpha", "Fix", out[2],priors.row(2), 2, 0);
        io_calls.fill_param(&width_in, "width:Appourchaux_v2:Gamma_alpha", "Fix", out[3]/10,priors.row(3), 3, 0);
        io_calls.fill_param(&width_in, "width:Appourchaux_v2:Wdip", "Fix", out[4],priors.row(4), 4, 0);
        io_calls.fill_param(&width_in, "width:Appourchaux_v2:DeltaGammadip", "Fix", out[5],priors.row(5), 5, 0);
	*/
	//io_calls.show_param(width_in, 0);
	return width_in;
}


// NEED TO WORK ON THE INPUT/OUTPUT OF THAT FUNCTION
VectorXd settings_aj_splittings(const int i, const MCMC_files inputs_MS_global, Input_Data* Snlm_in, const int aj_switch){
	int p0;
	VectorXd aj_param_count(6); // counter for knowing where we passed
	IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure
	
	aj_param_count.setZero();
	//std::cout << "inputs_MS_global.common_names[" << i << "] =" << inputs_MS_global.common_names[i] << "   aj_switch=" << aj_switch  << std::endl;
		if(inputs_MS_global.common_names[i] == "a1_0" && aj_switch == 6){  // This is a valid keyword only for aj models
			aj_param_count[0]=aj_param_count[0]+1;
			p0=0; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a1_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a1_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a1_1"  && aj_switch == 6){ // This is a valid keyword only for aj models
			aj_param_count[0]=aj_param_count[0]+1;
			p0=1; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a1_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a1_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a2_0" && aj_switch == 6){  // This is a valid keyword only for aj models
			aj_param_count[1]=aj_param_count[1]+1;
			p0=2; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a2_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a2_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a2_1"  && aj_switch == 6){ // This is a valid keyword only for aj models
			aj_param_count[1]=aj_param_count[1]+1;
			p0=3; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a2_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a2_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a3_0" && aj_switch == 6){  // This is a valid keyword only for aj models
			aj_param_count[2]=aj_param_count[2]+1;
			p0=4; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a3_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a3_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a3_1"  && aj_switch == 6){ // This is a valid keyword only for aj models
			aj_param_count[2]=aj_param_count[2]+1;
			p0=5; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a3_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a3_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a4_0" && aj_switch == 6){  // This is a valid keyword only for aj models
			aj_param_count[3]=aj_param_count[3]+1;
			p0=6; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a4_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a4_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a4_1"  && aj_switch == 6){ // This is a valid keyword only for aj models
			aj_param_count[3]=aj_param_count[3]+1;
			p0=7; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a4_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a3_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a5_0" && aj_switch == 6){  // This is a valid keyword only for aj models
			aj_param_count[4]=aj_param_count[4]+1;
			p0=8; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a5_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a4_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a5_1"  && aj_switch == 6){ // This is a valid keyword only for aj models
			aj_param_count[4]=aj_param_count[4]+1;
			p0=9; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a5_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a5_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a6_0" && aj_switch == 6){  // This is a valid keyword only for aj models
			aj_param_count[5]=aj_param_count[5]+1;
			p0=10; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a6_0", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a6_0 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(inputs_MS_global.common_names[i] == "a6_1"  && aj_switch == 6){ // This is a valid keyword only for aj models
			aj_param_count[5]=aj_param_count[5]+1;
			p0=11; // ONLY VALUD IF WE CONSIDER do_a11_eq_a12 == 0 && do_avg_a1n == 0 THIS SHOULD ALWAYS BE TRUE
			if(inputs_MS_global.common_names_priors[i] != "Fix_Auto"){ 
				io_calls.fill_param(Snlm_in, "a6_1", inputs_MS_global.common_names_priors[i], inputs_MS_global.modes_common(i,0), inputs_MS_global.modes_common.row(i), p0, 1);	
			} else{
				std::cout << "    Fix_Auto requested for a6_1 ... This is not allowed. Please use an explicit prior " << std::endl;
				exit(EXIT_SUCCESS);
			}
		}
	return aj_param_count;
}
