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
#include <boost/program_options.hpp>
#include "../../tamcmc/headers/gnuplot-iostream.h"
#include "../../tamcmc/headers/build_lorentzian.h"
#include "../../tamcmc/headers/function_rot.h"
#include "../../external/Alm/Alm_cpp/Alm_interpol.h"
#include "../../external/Alm/Alm_cpp/bilinear_interpol.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

namespace po = boost::program_options;

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
gsl_funcs Alm_grid_loader(const std::string grid_dir="../../external/Alm/data/Alm_grids_CPP/1deg_grids/");
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
	if (model_info.func_name[0] == "optimum_lorentzian_calc_ajAlm"){
		// params: H_l,  fc_l,  a1,  a3, a5, eta0, epsilon_nl, [theta0, delta], asym, gamma_l,  inclination, l
		l=(int)model_info.params[12];
		V=amplitude_ratio(l, model_info.params[11]);
		std::string ftype;
		if (model_info.params[13] == 0){
			//std::cout << " filter_type = 0 ==> Gate" << std::endl;
			ftype="gate";
		} else{
			//std::cout << " filter_type != 0 ==> Triangle" << std::endl;
			ftype="triangle";
		}
		gsl_funcs interp_funcs=Alm_grid_loader();
		//optimum_lorentzian_calc_ajAlm(x, y,        H_l          ,          fc_l       ,           a1        ,           a3         ,           a5         ,         eta0         ,         epsilon_nl  ,            thetas                ,        asym         ,       gamma_l        ,         l            ,       V  , step, c, ftype, gslfuncs);
		Optim_L model_tmp=optimum_lorentzian_calc_ajAlm(x, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6],    model_info.params.segment(7,2), model_info.params[9], model_info.params[10],         l            ,       V  , step, c, ftype, interp_funcs);
		VectorXd y(x.size());
		y.segment(model_tmp.i0, model_tmp.N)=y.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
	}
	if (model_info.func_name[0] == "optimum_lorentzian_calc_aj"){
		// params: H_l,  fc_l, a1,  a2, a3, a4, a5, a6, eta0, asym, gamma_l, inclination, l
		l=(int)model_info.params[12];
		V=amplitude_ratio(l, model_info.params[11]);
		Optim_L model_tmp=optimum_lorentzian_calc_aj(x, model_info.params[0], model_info.params[1], model_info.params[2],  model_info.params[3],  model_info.params[4],  model_info.params[5], model_info.params[6], model_info.params[7], model_info.params[8], model_info.params[9], model_info.params[10], l, V, step, c); 
		VectorXd y(x.size());
		y.segment(model_tmp.i0, model_tmp.N)=y.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
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
	const int Nmodels=7; // Total number of models
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
		struc.func_name[i]="optimum_lorentzian_calc_a1etaa3";
		struc.params_names[i]="H_l, fc_l, f_s,  eta,  a3,  asym, gamma_l, inclination, l";
		struc.Nparams[i]=9;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1etaa3_v2" || func_name == "get_all"){
		struc.func_name[i]="optimum_lorentzian_calc_a1etaa3_v2";
		struc.params_names[i]="H_lm,  fc_l,  f_s,  eta,   a3, asym, gamma_l, l";
		struc.Nparams[i]=8;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1l_etaa3_v2" || func_name == "get_all"){
		struc.func_name[i]="optimum_lorentzian_calc_a1l_etaa3_v2";
		struc.params_names[i]="H_lm,  fc_l,  f_s1,  f_s2,  eta,   a3, asym, gamma_l, l";
		struc.Nparams[i]=9;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1l_etaa3" || func_name == "get_all"){
		struc.func_name[i]="optimum_lorentzian_calc_a1l_etaa3";
		struc.params_names[i]="H_l,  fc_l,  f_s1,  f_s2,  eta,   a3, asym, gamma_l, inclination, l";
		struc.Nparams[i]=10;
		struc.error=false;
		i=i+1;
	}
	if (func_name == "optimum_lorentzian_calc_a1l_a2a3" || func_name == "get_all"){
		struc.func_name[i]="optimum_lorentzian_calc_a1l_a2a3";
		struc.params_names[i]="H_l,  fc_l,  f_s1,  f_s2,  a2, a3, asym, gamma_l, inclination, l";
		struc.Nparams[i]=10;
		struc.error=false;
		i=i+1;
	}	
	if (func_name == "optimum_lorentzian_calc_ajAlm" || func_name == "get_all"){
		struc.func_name[i]="optimum_lorentzian_calc_ajAlm";
		// filter_type : triangle=1  gate=0
		struc.params_names[i]="H_l,  fc_l,  a1,  a3, a5, eta0, epsilon_nl, theta0, delta, asym, gamma_l,  inclination, l, filter_type";
		struc.Nparams[i]=14;
		struc.error=false;
		i=i+1;
	}	
	if (func_name == "optimum_lorentzian_calc_aj" || func_name == "get_all"){
		struc.func_name[i]="optimum_lorentzian_calc_aj";
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
		struc.params.resize(struc.Nparams[0]); // Remove unecessary elements as Nparams defines them
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

int main(int argc, char* argv[]){
	int l;
	Func_info model_info, results;

	// Options 
    po::options_description desc("Program options");
    desc.add_options()
        ("get_all", "Return names of all the available models")
        ("help,H", "Show this help message")
        ("model,M", po::value<std::string>(), "Requested model")
		("xmin", po::value<double>(), "Minimum value for the frequency")
		("xmax", po::value<double>(), "Maximum value for the frequency")
        ("params,P", po::value< std::vector<double> >()->multitoken(), "All parameters for the function. Varies for each model. The separator is an empty space")
		("trunc,T", po::value<int>()->default_value(1000), "Truncation threshold for the Lorentzian")
		("resolution,R", po::value<double>()->default_value(0.05), "Spectrum resolution in microHz");
    
    po::variables_map vm;
	try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        po::notify(vm);  
    } catch(exception& e) {
        cout << "Error: " << e.what() << "\n";
        return -1;
    } 
	const double trunc=vm["trunc"].as<int>();
    const double resol=vm["resolution"].as<double>();
	 
	if (!vm.count("get_all") && (!vm.count("model") || !vm.count("xmin") || !vm.count("xmax") || !vm.count("params")))
    {
        std::cout << "The Number of arguments is incorrect.\n";
        std::cout << desc << std::endl;
        return -1;
    }

    if(vm.count("get_all")) {
		model_info=get_model_info("get_all");
		std::cout << std::setw(40) << "# available functions   " << std::setw(10) << "              required parameters" << std::endl;
		for (int i=0; i<model_info.func_name.size(); i++){
			std::cout << std::setw(40) << model_info.func_name[i] << "      " << std::setw(10) << model_info.params_names[i] << std::endl;
		}
    }

    if(vm.count("model")) { // Case where a function was explicitly requested and that initial checks of parameters length are passed
		Func_info model_info=get_model_info(vm["model"].as<std::string>());
		if (model_info.error == true){ // Stop execution if the model is tagged as having errors / incomplete / under construction / or whatever reason from the dev
			std::cout << " Error: The requested model is tagged as having errors. " << std::endl;
			std::cout << "        Cannot proceed " << std::endl;
			exit(EXIT_FAILURE);
		}
		model_info.xmin=vm["xmin"].as<double>();
		model_info.xmax=vm["xmax"].as<double>();
		try{
			for(int i=0; i<model_info.params.size();i++){
				model_info.params[i]=vm["params"].as<std::vector<double>>()[i];
			}
		} catch(...){
			std::cerr << "Unknown error while attempting to read parameters for the function: " << vm["model"].as<std::string>() << std::endl;
			std::cerr << "Debug required" << std::endl;
		}
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

gsl_funcs Alm_grid_loader(const std::string grid_dir){	
	// Handling the GSL interpolator pointer
	//const std::string grid_dir="external/Alm/data/Alm_grids_CPP/1deg_grids/"; // Using the provided 1 degree grid
	gsl_funcs funcs_data;
	GridData_Alm_fast grids;
	try{
		grids=loadAllData(grid_dir, "gate");
		// Pre-initialisation of the grid into gsl : Flattening + gsl init
		funcs_data.flat_grid_A10=flatten_grid(grids.A10);
		funcs_data.flat_grid_A11=flatten_grid(grids.A11);
		funcs_data.flat_grid_A20=flatten_grid(grids.A20);
		funcs_data.flat_grid_A21=flatten_grid(grids.A21);
		funcs_data.flat_grid_A22=flatten_grid(grids.A22);
		funcs_data.flat_grid_A30=flatten_grid(grids.A30);
		funcs_data.flat_grid_A31=flatten_grid(grids.A31);
		funcs_data.flat_grid_A32=flatten_grid(grids.A32);
		funcs_data.flat_grid_A33=flatten_grid(grids.A33);
		funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
		funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
		funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
		funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
		funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
		funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
		funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
		funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
		funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
		//
		/*
		grids=loadAllData(grid_dir, "gauss");
		// Pre-initialisation of the grid into gsl : Flattening + gsl init
		funcs_data.flat_grid_A10=flatten_grid(grids.A10);
		funcs_data.flat_grid_A11=flatten_grid(grids.A11);
		funcs_data.flat_grid_A20=flatten_grid(grids.A20);
		funcs_data.flat_grid_A21=flatten_grid(grids.A21);
		funcs_data.flat_grid_A22=flatten_grid(grids.A22);
		funcs_data.flat_grid_A30=flatten_grid(grids.A30);
		funcs_data.flat_grid_A31=flatten_grid(grids.A31);
		funcs_data.flat_grid_A32=flatten_grid(grids.A32);
		funcs_data.flat_grid_A33=flatten_grid(grids.A33);
		funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
		funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
		funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
		funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
		funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
		funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
		funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
		funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
		funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
		*/
		grids=loadAllData(grid_dir, "triangle");
		// Pre-initialisation of the grid into gsl : Flattening + gsl init
		funcs_data.flat_grid_A10=flatten_grid(grids.A10);
		funcs_data.flat_grid_A11=flatten_grid(grids.A11);
		funcs_data.flat_grid_A20=flatten_grid(grids.A20);
		funcs_data.flat_grid_A21=flatten_grid(grids.A21);
		funcs_data.flat_grid_A22=flatten_grid(grids.A22);
		funcs_data.flat_grid_A30=flatten_grid(grids.A30);
		funcs_data.flat_grid_A31=flatten_grid(grids.A31);
		funcs_data.flat_grid_A32=flatten_grid(grids.A32);
		funcs_data.flat_grid_A33=flatten_grid(grids.A33);
		funcs_data.interp_A10=init_2dgrid(funcs_data.flat_grid_A10);
		funcs_data.interp_A11=init_2dgrid(funcs_data.flat_grid_A11);
		funcs_data.interp_A20=init_2dgrid(funcs_data.flat_grid_A20);
		funcs_data.interp_A21=init_2dgrid(funcs_data.flat_grid_A21);
		funcs_data.interp_A22=init_2dgrid(funcs_data.flat_grid_A22);
		funcs_data.interp_A30=init_2dgrid(funcs_data.flat_grid_A30);
		funcs_data.interp_A31=init_2dgrid(funcs_data.flat_grid_A31);
		funcs_data.interp_A32=init_2dgrid(funcs_data.flat_grid_A32);
		funcs_data.interp_A33=init_2dgrid(funcs_data.flat_grid_A33);
	}
	catch (exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		std::cerr << "       It is possible that the grid files do not exist in the pointed directory" << std::endl;
		std::cerr << "       grid_dir is set to (relative path):" << grid_dir << std::endl;
		exit(EXIT_FAILURE);
	}
	return funcs_data;
}