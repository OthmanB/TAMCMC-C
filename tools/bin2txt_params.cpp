/*
   A small function that extract one parameter from the binary output
   That parameter is saved into a output file
*/

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include "diagnostics.h"
#include "data.h"
#include "string_handler.h"
#include "version.h"
#include "quick_samples_stats.h"

namespace po = boost::program_options;

void showversion();
bool file_exists (const std::string& name);
void usage(const po::options_description& desc);

int main(int argc, char* argv[]){

	const bool replicate_cte=false;
	bool istar=false; // since 1.85.0, this option is deactivated as old IDL code is obselete
	int ind0, indmax, ind_param, ind_var, ind_cons, ind_chain, Nsamples, Samples_period, Newsize, val, first_param, last_param;
	long cpt, lcpt;
	Eigen::MatrixXd data_array, data_out;
	std::string rootname, filename_params, filename_params_hdr, dir_out, file;
	std::ofstream fileout_stream, fileout_stream_synthese;
	std::string cpath=getcwd(NULL, 0);
	Eigen::VectorXd median, mean, stddev;
	Diagnostics diags;
	Params_hdr hdr;
	po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Print help message")
        ("version,v", "Print version information")
        ("rootname", po::value<std::string>()->required(), "The root name (along with the full path) of the binary/header file containing the parameters ('[out_root_name]_chain-*' file)")
        ("chain-index", po::value<int>()->required(), "The chain index (e.g. 0 for the coldest chain)")
        ("output-dir", po::value<std::string>()->required(), "The output directory (must already exist)")
        ("first-kept-element", po::value<int>()->default_value(0), "[Optional] Index of the first element which we keep. All index below that will be discarded")
        ("last-kept-element", po::value<int>()->default_value(-1), "[Optional] Index of the last element which we keep. All index below that will be discarded")
        ("periodicity", po::value<int>()->default_value(1), "[Optional] The Periodicity at which we keep samples. If <1 then all samples are returned")
        ("single-param-index", po::value<int>()->default_value(-1), "[Optional] If specified by a positive number, it will only extract the parameter with the provided index\n\tThis will not show create the Summary file");
	po::variables_map vm;
	try {
		po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
		if (vm.count("help")) {
			std::cerr << desc << std::endl;
			return EXIT_SUCCESS;
		}
		if (vm.count("version")) {
			showversion();
			return EXIT_SUCCESS;
		}
		po::notify(vm);
	}
	catch (const po::error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		usage(desc);
	}
	rootname = vm["rootname"].as<std::string>();
	ind_chain = vm["chain-index"].as<int>();
	dir_out = vm["output-dir"].as<std::string>();
	ind0 = vm["first-kept-element"].as<int>();
	indmax = vm["last-kept-element"].as<int>();
	Samples_period = vm["periodicity"].as<int>();
	
	filename_params=rootname + "_chain-" + int_to_str(ind_chain) + ".bin";
	if(!file_exists(filename_params.c_str())){
		filename_params=rootname + "_chain-" + int_to_str(ind_chain) + ".tar.gz";
		if(!file_exists(filename_params.c_str())){
			std::cerr << "Error: no binary or tar.gz file of this name found" << std::endl;
			std::cerr << "       check the usage or debug the program" << std::endl;
			std::cerr << "---" << std::endl;
			usage(desc);
		} else{
			std::cout << "  tar.gz file detected... proceeding in reading it..." << std::endl;
		}
		istar=true;
	} else{
		std::cout << "  bin file detected... proceeding in reading it..." << std::endl;
	}
	filename_params_hdr=rootname + ".hdr";

	std::cout << "  1. Reading the binary output files..." << std::endl;
	hdr=diags.read_params_header(filename_params_hdr); // Get the metadata from the header file 
	if(istar == false){
		data_array=diags.read_bin_matrix_dbl(filename_params, hdr.Nvars, hdr.Nsamples, "dbl");
	} else{
		data_array=diags.read_tar_gz_bin_matrix_dbl(filename_params, hdr.Nvars, hdr.Nsamples, "dbl");
	}
	median.resize(hdr.Nvars);
	mean.resize(hdr.Nvars);
	stddev.resize(hdr.Nvars);
	//
	if (indmax < 0){
		std::cout << "    Warning: Negative indmax provided ==> the last sample is going to be the last sample of the chain" << std::endl;
		indmax=data_array.rows();
	}
	if (indmax < ind0+Samples_period){
		std::cerr << "    Error: The last index must be greater than the first sample + Periodicity " << std::endl;
		std::cerr << "           The program will exit" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (Samples_period <=0){
		std::cout << "    Warning: Negative sampling period provided ==> sampling period will be fixed to 1 " << std::endl; 
		Samples_period=1;
	}
	// ---- Selecting only 1 sample out every Sample_period samples -----
	std::cout << "  2. Applying the selection rule..." << std::endl;
	if(Samples_period <= 1){ // Case where we return all samples
		data_out = data_array.block(ind0, 0, indmax-ind0, data_array.cols()) ;
	} else{ // Case where we keep one row out of Sample_period;
		Newsize=(indmax -ind0)/Samples_period;
		data_out.resize(Newsize, data_array.cols());
		cpt=ind0; lcpt=0;
		while(lcpt<data_out.rows()){
			data_out.row(lcpt)=data_array.row(cpt);
			cpt=cpt+Samples_period;
			lcpt=lcpt+1;
			if(lcpt == 1){
				std::cout << " Operation is of the type: ";
				std::cout << "Ouputs[" << lcpt << "]  --> New_Outputs[" << cpt << "]  "  << std::endl;
				std::cout << "..." << std::endl;
			}
		}
	}
	std::cout << "  3. Writing outputs into text files... " << std::endl;
	ind_var=0;
	ind_cons=0;
	fileout_stream_synthese.open((dir_out + "SUMMARY.STATS").c_str());
	if(fileout_stream_synthese.is_open()){
		fileout_stream_synthese << "# Statistical Summary of the MCMC analysis" << std::endl;
		fileout_stream_synthese << "#"  << std::setw(18) << "Variable_Name" << std::setw(18) << "Mean" << std::setw(18) << "Median" << std::setw(18) << "Stddev" << std::endl;
	} else{
		std::cerr << " Unable to open the binary data file " << dir_out.c_str() << "SUMMARY.STATS" << std::endl;	
		std::cerr << " Check that the full path exists" << std::endl;
		std::cerr << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);    			
	}
	if(vm["single-param-index"].as<int>() <= -1){ // We do all of the parameters
		first_param=0;
		last_param=hdr.Nvars+hdr.Ncons;
	} else{ 									 // Or we do the single index provided 
		if (vm["single-param-index"].as<int>() >= hdr.Nvars+hdr.Ncons){
			std::cerr << " Error: the 'single-param-index' parameter was specified, but with a value exceeding the maximum number of parameters" << std::endl;
			std::cerr << "        single-param-index = " <<  vm["single-param-index"].as<int>() << std::endl;
			std::cerr << "        Nvars + Ncons      = " <<  hdr.Nvars+hdr.Ncons << std::endl;
			exit(EXIT_FAILURE);
		}
		first_param=vm["single-param-index"].as<int>();
		last_param=vm["single-param-index"].as<int>() + 1;
	}
	for(int ind_param=0; ind_param < hdr.Nvars+hdr.Ncons; ind_param++){
		if (ind_param >= first_param && ind_param <last_param){ //
			file=diags.formated_int_to_str(ind_param);
			fileout_stream.open((dir_out + file + ".ASCII").c_str());
			if(fileout_stream.is_open()){
				if(hdr.relax[ind_param] == 1){
					fileout_stream << "! variable_name= " << hdr.variable_names[ind_var] << std::endl;
					fileout_stream << std::setprecision(12) << data_out.col(ind_var) << std::endl;
					mean[ind_var]=mean_fct(data_out.col(ind_var));
					median[ind_var]=median_fct(data_out.col(ind_var));
					stddev[ind_var]=stddev_fct(data_out.col(ind_var));
					fileout_stream_synthese << " " << std::setw(18) <<  hdr.variable_names[ind_var] << std::setw(18) << std::setprecision(12) << mean[ind_var] <<  std::setw(18) << std::setprecision(12) << median[ind_var] << std::setw(18) << std::setprecision(12) << stddev[ind_var] <<  std::endl;
					std::cout << "    - File: " << file + ".ASCII" << "    variable: " << hdr.variable_names[ind_var] << "   median: " << median[ind_var] << "   stddev: " << stddev[ind_var] << std::endl;
					ind_var=ind_var+1;
				} else{
					fileout_stream << "! constant_name= " << hdr.constant_names[ind_cons] << std::endl;
					fileout_stream_synthese << " " << std::setw(18) << hdr.constant_names[ind_cons] << std::setw(18) << std::setprecision(12) << hdr.constant_values[ind_cons] << std::setw(18) << std::setprecision(12) << hdr.constant_values[ind_cons] << std::setw(18) << "0"<<  std::endl;
					std::cout << "    - File: " << file + ".ASCII" << "  (Constant value)" << "    variable: " << hdr.constant_names[ind_cons] << "   value: " << hdr.constant_values[ind_cons] << std::endl;
					if(replicate_cte == 1){
						for(long repeat=0; repeat<data_array.rows(); repeat++){
							fileout_stream << hdr.constant_values[ind_cons] << std::endl;    
						}
					} else{
							fileout_stream << hdr.constant_values[ind_cons] << std::endl;
					}
					ind_cons=ind_cons + 1;				
				}
			} else{
				std::cerr << " Unable to open the binary data file " << dir_out.c_str() + file + ".ASCII" << std::endl;	
				std::cerr << " Check that the full path exists" << std::endl;
				std::cerr << " The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
			}
			fileout_stream.close();
		} else{
			if(hdr.relax[ind_param] == 1){
				ind_var=ind_var+1;
			} else{
				ind_cons=ind_cons + 1;	
			}
		}
	}
	
	// --- Write a single line with the plength vector on a file ---
	file=dir_out + "plength.txt";
	fileout_stream.open(file.c_str());
	if(fileout_stream.is_open()){
			fileout_stream << hdr.plength << std::endl;
	} else{
		std::cerr << " Unable to open the binary data file " << dir_out.c_str() + file  << std::endl;	
		std::cerr << " Check that the full path exists" << std::endl;
		std::cerr << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
	fileout_stream.close();

	if(vm["single-param-index"].as<int>() == -1){ // We do all of the parameters
		// --- Show a summary of all the values in a list format that can be easily taken to python ---
		std::cout << "  4. Statistics of outputs in simple format (variable names, mean, median, standard deviation) " << std::endl;
		std::cout << "      ";
		for (int i=0; i<hdr.Nvars;i++){
			std::cout << hdr.variable_names[i] <<" , " ;
		}
		std::cout << std::endl;
		std::cout << "      ";
		for (int i=0; i<hdr.Nvars;i++){
			std::cout << mean[i];
			if (i != hdr.Nvars-1){
				std::cout <<" , " ;
			} 
		}
		std::cout << std::endl;
		std::cout << "      ";
		for (int i=0; i<hdr.Nvars;i++){
			std::cout << median[i] ;
			if (i != hdr.Nvars-1){
				std::cout <<" , " ;
			}     	
		}
		std::cout << std::endl;
		std::cout << "      ";
		for (int i=0; i<hdr.Nvars;i++){
			std::cout << stddev[i] ;
			if (i != hdr.Nvars-1){
				std::cout <<" , " ;
			} 
		}
		std::cout << std::endl;
	} else{
		std::cout << "  4. Single parameter extraction requested... not showing a statistics table" << std::endl;
	}
}


void showversion()
{
    std::cout << APP_NAME " - bin2txt tool - " APP_VERSION "\n built on " __DATE__ << std::endl;

#   if defined(__clang__)
    	printf(" with clang " __clang_version__);
#   elif defined(__GNUC__)
    	printf(" with GCC");
    	printf(" %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#   elif defined(_MSC_VER)
    	printf(" with MSVC");
    	printf(" %d", MSVC_VERSION);
#   else
    printf(" unknown compiler");
#   endif

    std::cout << "\n features:";
#   if defined(__i386__) || defined(_M_IX86)
    std::cout << " i386" << std::endl;
#   elif defined(__x86_64__) || defined(_M_AMD64)
    std::cout << " x86_64" << std::endl;
#   endif
    std::cout << " Author: " << APP_COPYRIGHT << std::endl;

}

void usage(const po::options_description& desc){

    std::cout << " You need to provide at least 3 arguments to that function. The available arguments are: " << std::endl;
    std::cout << desc << std::endl;
	std::cout << "      WARNING: Since 1.85.0, Calls to arguments was drastically changed.\n\t Any code that calls bin2txt need to be upgraded accordingly" << std::endl;

    exit(EXIT_FAILURE);
}

bool file_exists(const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}