#include <Eigen/Dense>
#include <vector>
#include <string>

#include "data.h" // contains the structure Data
#include "string_handler.h"
#include "io_models.h" // Contains the function that create the final vector

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

aj_files read_ajfile(const std::string cfg_model_file){
	/*
	This function is designed to read .model files for fitting aj parameters
	using a Gaussian Likelihood. 
	The syntax is strict to stay simple. 
	*/
	aj_files i_ajfit;
	std::string line, line0, subline0; //token
	std::vector<std::string> header, words;
	std::ifstream cfg_session;
	int cpt=0;
    cfg_session.open(cfg_model_file.c_str());
    if (cfg_session.is_open()) {
		// [1] Get the header
		std::getline(cfg_session, line0);
		line0=strtrim(line0); // remove any white space at the begining/end of the string
		subline0=strtrim(line0.substr(0, 1)); // pick the first character
		subline0=subline0.c_str();
		if (subline0 == "#"){
			while(subline0 == "#"){		
				header.push_back(strtrim(line0.substr(1, std::string::npos))); // add all characters except the first one (and any space at the begining)
				std::getline(cfg_session, line0);
				line0=strtrim(line0); // remove any white space at the begining/end of the string
				subline0=strtrim(line0.substr(0, 1)); // pick the first character
				subline0=subline0.c_str();
				cpt=cpt+1;
			}
			std::cout << "   [1] " << cpt << " header lines found..." << std::endl;
		} else{
			header.push_back("");
			std::cout << "   [1] Header not found. Header vector set to a blank vector<string> of size 1. Pursuing operations..." << std::endl;
		}
		std::cout << "     [2] Reading Dnu..." << std::endl;
		std::getline(cfg_session, line0);
		i_ajfit.Dnu=str_to_dbl(line0);
		std::cout << "          - Dnu = " << i_ajfit.Dnu << std::endl;
		std::cout << "     [3] els..." << std::endl;
		std::getline(cfg_session, line0); // First skip the '! l' comment
		std::getline(cfg_session, line0);
		i_ajfit.els=str_to_Xdarr(line0, " ");
		std::cout << "          - els = " << i_ajfit.els.transpose() << std::endl;
		std::cout << "     [4] nu_nl..." << std::endl;
		std::getline(cfg_session, line0); // First skip the '! nu_nl' comment
		std::getline(cfg_session, line0);
		i_ajfit.nu_nl=str_to_Xdarr(line0, " ");
		std::cout << "          - nu_nl = " << i_ajfit.els.transpose() << std::endl;
		std::cout << "     [5] setups..." << std::endl;
		std::getline(cfg_session, line0); // First skip the '! setups' comment
		std::getline(cfg_session, line0); // filter_type
		words=strsplit(line0, " \t");
		i_ajfit.filter_type=strtrim(words[1]); // gate or gauss or triangle
		std::cout << "          - filter_type = " << i_ajfit.filter_type << std::endl;
		std::getline(cfg_session, line0); // do_a2
		words=strsplit(line0, " \t");
		i_ajfit.do_a2=str_to_bool(words[1]);
		std::getline(cfg_session, line0); // do_a4
		words=strsplit(line0, " \t");
		i_ajfit.do_a4=str_to_bool(words[1]);
		std::getline(cfg_session, line0); // do_a6
		words=strsplit(line0, " \t");
		i_ajfit.do_a6=str_to_bool(words[1]);
		std::cout << "          - do_a2 = " << i_ajfit.do_a2 << std::endl;
		std::cout << "          - do_a4 = " << i_ajfit.do_a4 << std::endl;
		std::cout << "          - do_a6 = " << i_ajfit.do_a6 << std::endl;
		std::getline(cfg_session, line0); // do_CFonly
		words=strsplit(line0, " \t");
		i_ajfit.do_CFonly=str_to_bool(words[1]);
		std::cout << "          - do_CFonly = " << i_ajfit.do_CFonly << std::endl;
   		cfg_session.close();
   } else {
   		std::cout << "Unable to open the file: " << cfg_model_file << std::endl;
   		std::cout << "Check that the file exist and that the path is correct" << std::endl;
   		std::cout << "Cannot proceed" << std::endl;
   		std::cout << "The program will exit now" << std::endl;
   		exit(EXIT_FAILURE);
   }
	return i_ajfit;
}

/*
Data set_observables_ajfit(aj_files i_ajfit, Data data){
	//
	//	Function that removes irrelevant aj terms from that are specified in the data file
	//	in order to only keep a2, a4, a6 and removing some of these in function of do_a2, do_a4, do_a6
	//
	const short int Nobs=i_ajfit.do_a2 + i_ajfit.do_a4 + i_ajfit.do_a6; // A sum of 0 and 1 gives the number of observables
	int j;
	Data data_out;

	data_out.x.resize(Nobs);
	data_out.y.resize(Nobs);
	data_out.sigma_y.resize(Nobs);
	data_out.xlabel=data.xlabel;
	data_out.ylabel=data.ylabel;
	data_out.xunit=data.xunit;
	data_out.yunit=data.yunit;
	data_out.header=data.header;

	j=0;
	for (int i = 0; i<data.x.size(); i++){
		if ((data.x[i] == 2 && i_ajfit.do_a2 == true) ||  (data.x[i] == 4 && i_ajfit.do_a4 == true) || (data.x[i] == 6 && i_ajfit.do_a6 == true) ){
			data_out.x[j]=data.x[i];
			data_out.y[j]=data.y[i];
			data_out.sigma_y[j]=data.sigma_y[i];
			j=j+1;
		}
	}
	std::cout << data.x.transpose() << std::endl;
	std::cout << data.y.transpose() << std::endl;
	std::cout << data.sigma_y.transpose() << std::endl;
	
	std::cout << "data_out.x = " << data_out.x.transpose() << std::endl;
	std::cout << "data_out.y = " << data_out.y.transpose() << std::endl;
	std::cout << "data_out.sigma_y = " << data_out.sigma_y.transpose() << std::endl;
	std::cout << " Test in set_observables_ajfit " << std::endl;
	exit(EXIT_SUCCESS);
	return data_out;
}
*/

Data_Nd set_observables_ajfit(const aj_files i_ajfit, const Data_Nd data, const int x_col, const int y_col, const int ysig_col){
	/*
		Function that removes irrelevant aj terms from that are specified in the data file
		in order to only keep a2, a4, a6 and removing some of these in function of do_a2, do_a4, do_a6
	*/
	const short int Nobs=i_ajfit.do_a2 + i_ajfit.do_a4 + i_ajfit.do_a6; // A sum of 0 and 1 gives the number of observables
	const int Nrows=data.data.rows(), Ncols=data.data.cols();
	int j;
	Data_Nd data_out;

	data_out.labels=data.labels;
	data_out.units=data.units;
	data_out.header=data.header;
	data_out.data.resize(Nobs, Ncols);

	j=0;
	for (int i = 0; i<Nrows; i++){
		if ((data.data(i,x_col) == 2 && i_ajfit.do_a2 == true) ||  (data.data(i,x_col) == 4 && i_ajfit.do_a4 == true) || (data.data(i,x_col) == 6 && i_ajfit.do_a6 == true) ){
			data_out.data(j, x_col)=data.data(i, x_col);
			data_out.data(j, y_col)=data.data(i, y_col);
			data_out.data(j, ysig_col)=data.data(i, ysig_col);
			j=j+1;
		}
	}	
	return data_out;
}

Input_Data build_init_ajfit(const aj_files i_ajfit, const double a1_obs){
	/*
		The main function that structure the inputs in a way that can be handled by the MCMC code and 
		also by the model "model_ajfit"
	*/
		const int Nmax_prior_params=4; // The maximum number of parameters for the priors. Should be 4 in all my code
		// Hard-coded initial guess... We have so few parameters that it is irrelevant to let the user set this... 
		const double epsilon_nl0=5e-4;
		const double theta0=70;  // in degrees
		const double delta=10; // in degrees
		std::vector<double> els, nu_nl;
		VectorXi plength(5);
		VectorXd priors_epsilon_nl0(4), priors_theta0(4), priors_delta(4), extra_priors(1);
		VectorXd tmpXd(4);
		std::vector<bool> relax;
		Input_Data all_in;
		IO_models io_calls; // function dictionary that is used to initialise, create and add parameters to the Input_Data structure
		int p0;
		tmpXd  << -9999, -9999, -9999, -9999;
		extra_priors(0)=-9999;
		
		all_in.model_fullname = "model_ajfit";
		priors_epsilon_nl0 << 1e-3, 1e-2, -9999, -9999; // FOR JEFFREYS
		priors_theta0 << 0, 90, -9999, -9999; // FOR UNIFORM
		priors_delta  << 10, 45, -9999, -9999; 
		if (i_ajfit.filter_type == "gate"){ // gate is coded by 0
			priors_delta  << 10, 45, -9999, -9999;// FOR JEFFREYS DEFAULT (gate and triangle case). gauss requires smaller range (3 times smaller)
		} else{
			priors_delta  << 3, 45, -9999, -9999; // The sigma of a gauss function is approx 3-5 times smaller than the pure slope of a triangle function
		}
		plength[0]=3; // First epsilon_nl0 / theta0 / delta
		plength[1]=2; // Second Dnu and a1 for the centrifugal term
		plength[2]=i_ajfit.els.size();
		plength[3]=i_ajfit.nu_nl.size();
		plength[4]=5;
		if(plength[2] != plength[3]){
			std::cout << " ERROR: The number of els and the number of nu_nl do not match inside the model file!" << std::endl;
			std::cout << "        Please check your model file " << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "plength =" << plength.transpose() << std::endl;
		// --- conversion --
		for (int i=0; i<plength[2];i++){
			els.push_back(i_ajfit.els[i]);
			nu_nl.push_back(i_ajfit.nu_nl[i]);
			relax.push_back(0);
		}
		std::cout << "Init all_in..." << std::endl;	
		io_calls.initialise_param(&all_in, plength.sum(), Nmax_prior_params, plength, extra_priors);	
		std::cout << "Setting initial guesses and priors for the parameters epsilon_nl0, theta0, delta" << std::endl;
		p0=0;
		io_calls.fill_param(&all_in, "epsilon_0", "Jeffreys", epsilon_nl0, priors_epsilon_nl0, p0, 0);
		p0=1;
		io_calls.fill_param(&all_in, "theta0", "Uniform", theta0, priors_theta0, p0, 0);
		p0=2;
		io_calls.fill_param(&all_in, "delta", "Jeffreys", delta, priors_delta, p0, 0);
		p0=3;
		io_calls.fill_param(&all_in, "Dnu", "Fix", i_ajfit.Dnu, tmpXd, p0, 0);
		p0=4;
		io_calls.fill_param(&all_in, "a1_obs", "Fix", a1_obs, tmpXd, p0, 0);
		p0=5;
		io_calls.fill_param_vect(&all_in, els, relax, "els", "Fix", tmpXd, p0, 0, 0);
		p0=	5 + i_ajfit.els.size();
		io_calls.fill_param_vect(&all_in, nu_nl, relax, "nu_nl", "Fix", tmpXd, p0, 0, 0);
		p0= 5 + 2*i_ajfit.els.size();
		io_calls.fill_param(&all_in, "do_a2", "Fix", i_ajfit.do_a2, tmpXd, p0, 0);
		p0= 5 + 2*i_ajfit.els.size() + 1;
		io_calls.fill_param(&all_in, "do_a4", "Fix", i_ajfit.do_a4, tmpXd, p0, 0);
		p0= 5 + 2*i_ajfit.els.size() + 2;
		io_calls.fill_param(&all_in, "do_a6", "Fix", i_ajfit.do_a6, tmpXd, p0, 0);
		p0= 5 + 2*i_ajfit.els.size() + 3;
		io_calls.fill_param(&all_in, "do_CFonly", "Fix", i_ajfit.do_CFonly, tmpXd, p0, 0);
	 	// Handling the filter type
		p0= 5 + 2*i_ajfit.els.size() + 4;
		if (i_ajfit.filter_type == "gate"){ // gate is coded by 0
			io_calls.fill_param(&all_in, "filter_type", "Fix", 0, tmpXd, p0, 0);
		}
		if (i_ajfit.filter_type == "gauss"){ // gauss is coded by 1
			io_calls.fill_param(&all_in, "filter_type", "Fix", 1, tmpXd, p0, 0);
		}
		if (i_ajfit.filter_type == "triangle"){ // triangle is coded by 2
			io_calls.fill_param(&all_in, "filter_type", "Fix", 2, tmpXd, p0, 0);
		}
		if (i_ajfit.filter_type != "gate" && i_ajfit.filter_type != "gauss" && i_ajfit.filter_type != "triangle"){
			std::cout << " Error: Unrecognized filter type. You must have filter_type = gate or gauss or triangle" << std::endl;
			exit(EXIT_SUCCESS);
		}
		std::cout << " ----------------- Configuration summary -------------------" << std::endl;
		std::cout << " Model Name = " << all_in.model_fullname << std::endl;		
		std::cout << " -----------------------------------------------------------" << std::endl;
		std::cout << " ---------- Configuration of the input vectors -------------" << std::endl;
		std::cout << "    Table of inputs " << std::endl;
		io_calls.show_param(all_in, 1);
		std::cout << " -----------------------------------------------------------" << std::endl;
	return all_in;
}