/*
Date: 16 Nov 2021
Testing of the acoefs.cpp function present in [root]/tamcmc/sources and its header 
Adapted from acoefs.py from the acoefs_check project (see github)
This code verifies creates an executable that allow to verify that when
you provide a central frequency nu_nl, along with aj coefficients j={1,...,6}
(1) the decomposition in acoeficient is properly done and along the way
(2) Tnlm is properly calculated
(3) Snlm is properly calculated
*/
#include <math.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
//#include "../../tamcmc/headers/acoefs.h"
#include "../../tamcmc/headers/build_lorentzian.h"
//#include "../../tamcmc/headers/string_handler.h"
#include "../../tamcmc/headers/function_rot.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

struct Func_info{
	bool error;
	double xmin;
	double xmax;
	std::vector<std::string> func_name;
	std::vector<int> Nparams;
	std::vector<std::string> params_names;
	VectorXd params;
	VectorXd x;
	VectorXd y;
};

void usage(int argc, char* argv[]);
int  options(int argc, char* argv[]);
Func_info test_lorentzian(const Func_info& model_info, const double trunc);
Func_info get_model_info(std::string func_name, const bool ignore_error_msg=false);
Func_info call_function(const VectorXd x, const Func_info& model_info, const double c);


Func_info call_function(const VectorXd x,  const Func_info& model_info, const double c){
	bool OK;
	int l;
	double step=x[10] - x[9];
	VectorXd y(x.size()), V;
	y.setZero();
	Func_info model_out;
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1etaa3"){
		// params: H_l, fc_l, f_s,  eta,  a3,  asym, gamma_l, inclination, l
		l=(int)model_info.params[8];		
		V=amplitude_ratio(l, model_info.params[7]);
		y=optimum_lorentzian_calc_a1etaa3(x, y, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6], l, V, step, c); 
//VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1acta3"){
		// params: H_l,  fc_l,  f_s,  eta,  a3, b,  alpha,  asym,  gamma_l, inclination, l
		l=(int)model_info.params[10];
		V=amplitude_ratio(l, model_info.params[9]);
		y=optimum_lorentzian_calc_a1acta3(x, y, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6], model_info.params[7], model_info.params[8], l, V, step, c); 
 //VectorXd optimum_lorentzian_calc_a1acta3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s,  const double eta,  const double a3, 
		 //const double b,  const double alpha,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1a2a3"){
		// params: H_l,  fc_l,  f_s,  a2,   a3, asym, gamma_l, inclination, l
		l=(int)model_info.params[8];
		V=amplitude_ratio(l, model_info.params[7]);
		y=optimum_lorentzian_calc_a1a2a3(x, y, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6], l, V, step, c); 
 //VectorXd optimum_lorentzian_calc_a1a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1etaa3_v2"){
		// params: H_l, fc_l, f_s,  eta,  a3,  asym, gamma_l, l 
		l=(int)model_info.params[2*l+7];
		//V=amplitude_ratio(l, model_info.params[7]); 
		y=optimum_lorentzian_calc_a1etaa3_v2(x, y, model_info.params.segment(0,2*l+1), model_info.params[2*l+1], model_info.params[2*l+2],  model_info.params[2*l+3],  model_info.params[2*l+4],  model_info.params[2*l+5], model_info.params[2*l+6], l, step, c); 
 //VectorXd   optimum_lorentzian_calc_a1etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1l_etaa3_v2"){
		// params: H_lm,  fc_l,  f_s,  eta,   a3, asym, gamma_l, inclination, l
		l=(int)model_info.params[8];
		//V=amplitude_ratio(l, model_info.params[8]);
		y=optimum_lorentzian_calc_a1l_etaa3_v2(x, y, model_info.params.segment(0,2*l+1), model_info.params[2*l+1], model_info.params[2*l+2],  model_info.params[2*l+3],  model_info.params[2*l+4],  model_info.params[2*l+5], model_info.params[2*l+6], model_info.params[2*l+7], l, step, c); 
//VectorXd  optimum_lorentzian_calc_a1l_etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l, double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1l_etaa3"){
		// params: H_lm,  fc_l,  f_s1,  f_s2,  eta,   a3, asym, gamma_l, inclination, l
		l=(int)model_info.params[9];
		V=amplitude_ratio(l, model_info.params[8]);
		y=optimum_lorentzian_calc_a1l_etaa3(x, y, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6], model_info.params[7], l, V, step, c); 
//VectorXd  optimum_lorentzian_calc_a1l_etaa3(const VectorXd& x, const VectorXd& y, const  double H_l, const  double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1l_a2a3"){
		// params: H_l,  fc_l,  f_s1,  f_s2,  a2, a3, asym, gamma_l, inclination, l
		l=(int)model_info.params[9];
		V=amplitude_ratio(l, model_info.params[8]);
		y=optimum_lorentzian_calc_a1l_a2a3(x, y, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6], model_info.params[7], l, V, step, c); 
 //VectorXd optimum_lorentzian_calc_a1l_a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_a1etaAlma3"){
		// params: H_l,  fc_l,  f_s, eta0, epsilon_nl, [theta0, delta], a3, asym, gamma_l, inclination, l
		l=(int)model_info.params[12];
		V=amplitude_ratio(l, model_info.params[11]);
		y=optimum_lorentzian_calc_a1etaAlma3(x, y, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params.segment(6,2), model_info.params[8], model_info.params[9], model_info.params[10], l, V, step, c); 
// VectorXd optimum_lorentzian_calc_a1etaAlma3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s, 
	//	const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_aj"){
		// params: H_l,  fc_l, a1,  a2, a3, a4, a5, a6, eta0, asym, gamma_l, inclination, l
		l=(int)model_info.params[12];
		V=amplitude_ratio(l, model_info.params[11]);
		y=optimum_lorentzian_calc_aj(x, y, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6], model_info.params[7], model_info.params[8], model_info.params[9], model_info.params[10], l, V, step, c); 
// VectorXd optimum_lorentzian_calc_aj(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, 
//       		const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
//      		const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);
	}
	model_out.func_name=model_info.func_name;
	model_out.error=model_info.error;
	model_out.params=model_info.params;
	model_out.params_names=model_info.params_names;
	model_out.Nparams=model_info.Nparams;
	OK=true;
	if (OK == true){
		//model_out.x.resize(x.size());
		//model_out.y.resize(x.size())
		model_out.x=x;
		model_out.y=y;	
		model_out.error=false;
		model_out.xmin=x.minCoeff();
		model_out.xmax=x.maxCoeff();
	} else{
		model_out.error=true;
	}
	return model_out;
}

Func_info get_model_info(std::string func_name, const bool ignore_error_msg){
	/* 
		A function that allows to retrieve the specifics of a single model if func_name is an existing model name
		Or will retrieve the full list of available models if func_name is set to "get_all"
	*/
	const int Nmodels=9;
	int i=0;
	Func_info struc;
	
	struc.error=false;

	if (func_name != "get_all"){
		struc.func_name.resize(1);
		struc.params_names.resize(1);
		struc.Nparams.resize(1);
		struc.params.resize(1);
	} else{
		struc.func_name.resize(Nmodels);
		struc.params_names.resize(Nmodels);
		struc.Nparams.resize(Nmodels);
	}
	if (func_name == "optimum_lorentzian_calc_a1etaa3" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1etaa3";
		}
		struc.params_names[i]="H_l, fc_l, f_s,  eta,  a3,  asym, gamma_l, inclination, l";
		struc.Nparams[i]=9;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1acta3" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1acta3";
		}		
		struc.params_names[i]="H_l,  fc_l,  f_s,  eta,  a3, b,  alpha,  asym,  gamma_l, inclination, l";
		struc.Nparams[i]=11;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1a2a3" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1a2a3";
		}		
		struc.params_names[i]="H_l,  fc_l,  f_s,  a2,   a3, asym, gamma_l, inclination, l";
		struc.Nparams[i]=9;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1etaa3_v2" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1etaa3_v2";
		}		
		struc.params_names[i]="H_lm,  fc_l,  f_s,  eta,   a3, asym, gamma_l, l";
		struc.Nparams[i]=8;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1l_etaa3_v2" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1l_etaa3_v2";
		}		
		struc.params_names[i]="H_lm,  fc_l,  f_s1,  f_s2,  eta,   a3, asym, gamma_l, l";
		struc.Nparams[i]=9;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1l_etaa3" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1l_etaa3";
		}		
		struc.params_names[i]="H_l,  fc_l,  f_s1,  f_s2,  eta,   a3, asym, gamma_l, inclination, l";
		struc.Nparams[i]=10;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1l_a2a3" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1l_a2a3";
		}		
		struc.params_names[i]="H_l,  fc_l,  f_s1,  f_s2,  a2, a3, asym, gamma_l, inclination, l";
		struc.Nparams[i]=10;
		struc.error=false;
		i=i+1;
	}	
	if (func_name == "optimum_lorentzian_calc_a1etaAlma3" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
			std::cout << " WARNING: THE FUNCTION " << func_name << " IS STILL EVOLVING " << std::endl;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_a1etaAlma3";
		}		
		struc.params_names[i]="H_l,  fc_l,  f_s, eta0, epsilon_nl, theta0, delta, a3, asym, gamma_l, inclination, l";
		struc.Nparams[i]=12;
		struc.error=true;
		i=i+1;
	}	
	if (func_name == "optimum_lorentzian_calc_aj" || func_name == "get_all"){
		if (func_name != "get_all"){
			struc.func_name[i]=func_name;
		} else{
			struc.func_name[i]="optimum_lorentzian_calc_aj";
		}		
		struc.params_names[i]="H_l,  fc_l, a1,  a2, a3, a4, a5, a6, eta0, asym, gamma_l, inclination, l";
		struc.Nparams[i]=13;
		struc.error=false;
		i=i+1;
	}	
	if (struc.error == true && ignore_error_msg == false){
		std::cout << "Error: The name provided in get_params_names() is Unrecognized: " << func_name << std::endl;
		std::cout << "       Please be sure to enter a name of an actual model in build_lorentzian.cpp" << std::endl;
		std::cout << "       Note that only models starting by 'optimum_[...]' are recognized here" << std::endl;
		std::cout << "       The program will exit now" << std::endl;
		//exit(EXIT_FAILURE);
	}
	if (func_name == "get_all"){ // Bypassing the error indicator related to each function 
		struc.error=false;
	} else{
		struc.params.resize(struc.Nparams[0]);
	}
	return struc;
}	

Func_info test_lorentzian(const Func_info model_info, const double trunc, const double resol){
	/*
		With this function you can call any model defined into the build_lorentzian function and get its evaluation over a set x of frequencies and
		provided a set of inputs. The inputs are given as a vector and you need to provide a consistent number of them, in the order these appear into the
		function definition that is called
	*/
	const int Ndata=(model_info.xmax-model_info.xmin)/resol;
	VectorXd x(Ndata);
	Func_info results;

	x.setLinSpaced(Ndata, model_info.xmin, model_info.xmax);
	results=call_function(x, model_info, trunc);
	return results;
}

void usage(int argc, char* argv[]){
			std::cout << "Help" << std::endl;
			std::cout << "     - To execute: " << argv[0] << "  <function> <xmin> <xmax> <params> "<< std::endl;
			std::cout << "                   Note that if <function> is set to 'get_all', all other inputs are ignored and the program returns a list of available functions to test" << std::endl;
			exit(EXIT_FAILURE);
}

int options(int argc, char* argv[]){
	// Check that we have the correct number of arguments
	int val=-1;
	std::string func_name="";
	Func_info model_info;

	if (argc == 2){
		func_name = argv[1];
		if (func_name == "get_all"){
			val=2;
			return val;
		}
		if (func_name == "help"){
			usage(argc, argv);
			return 0;
		}
	}

	if (argc > 2){
		func_name=argv[1];
		model_info=get_model_info(func_name); // Get the information on the requested model in order to do a first consistency check on the number of parameters to run it	
		if (model_info.error == true){ // Stop execution if the model is tagged as having errors / incomplete / under construction / or whatever reason from the dev
			std::cout << " Error: The requested model is tagged as having errors. " << std::endl;
			std::cout << "        Cannot proceed " << std::endl;
			exit(EXIT_FAILURE);
		}
		//std::cout << "model_info.Nparams[0] + 4 = " << model_info.Nparams[0] + 4 << std::endl;
		if(argc == model_info.Nparams[0] + 4){ // We expect to have argv with the (0) name of the program, (1) function name, (2) xmin, (3) xmax, and (4) All parameters of the function
			val=1; 
		} else{
			std::cout << " Error: Unexpected number of parameters provided to run the function " << func_name << " in call_function()" << std::endl;
			std::cout << "        Count                           : " << argc - 4 << std::endl;
			std::cout << "        Expected                        : " << model_info.Nparams[0] << std::endl;
			std::cout << "        Name of the expected parameters : " << model_info.params_names[0] << std::endl;
			std::cout << "        Please check (1) that you did not do a mistake in the number of arguments for the parameters" << std::endl;
			std::cout << "                     (2) that the name of the function you wish to execute is correct" << std::endl;
			std::cout << "                     (3) Check that the order of the parameters is also the same as specified above " << std::endl;
		}
	}
	if (val == -1){ // Error code
		usage(argc, argv);
	} 
	if (val > 0 ){
		return val; // Execution code val
	} else{
		return -1; // Default value is to return an error code
	}
}

Func_info read_inputs(int argc, char* argv[]){
	/*
		A function that takes the argv and organise them in a way readable by test_lorentzian()
	*/
	std::string func_name;
	Func_info model_info;

	model_info=get_model_info(argv[1]);

	std::istringstream(argv[2]) >> model_info.xmin;
	std::istringstream(argv[3]) >> model_info.xmax;
	for (int i=0; i<model_info.params.size();i++){
		std::istringstream(argv[i+4]) >> model_info.params[i];
	}
	return model_info;
}
int main(int argc, char* argv[]){
	const double trunc=1000;
	const double resol=0.05; //resolution in microHz is fixed
	int msg_code;
	int l;
	Func_info model_info, results;

	msg_code=options(argc, argv);
	
	if(msg_code == -1){
		std::cout << "Error detected in options. Cannot proceed. Debug required." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		if (msg_code == 2){ // Case where the function name was "get_all" >> return all the models only
			model_info=get_model_info("get_all");
			std::cout << std::setw(40) << "# available functions   " << std::setw(10) << "              required parameters" << std::endl;
			for (int i=0; i<model_info.func_name.size(); i++){
				std::cout << std::setw(40) << model_info.func_name[i] << "      " << std::setw(10) << model_info.params_names[i] << std::endl;
			}
		}
		if (msg_code == 1){ // Case where a function was explicitly requested and that initial checks of parameters length are passed
			model_info=read_inputs(argc, argv);

			results=test_lorentzian(model_info, trunc, resol);
			std::cout << "# xmin= " << results.xmin << std::endl;
			std::cout << "# xmax= " << results.xmax << std::endl;
			std::cout << "# func_name = " << results.func_name[0] << std::endl;
			std::cout << "# params_names = " << results.params_names[0] << std::endl;
			std::cout << "! params = " << results.params.transpose() << std::endl;
			std::cout << "# Freq     Power " << std::endl;
			for (int i = 0; i<results.x.size(); i++){
				std::cout <<  std::setw(20) << std::setprecision(10) << results.x[i] <<  std::setw(20) << std::setprecision(10) << results.y[i] << std::endl;
			}
		}
	}
}

