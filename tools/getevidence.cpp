/*
   A small function that extract the likelihoods / priors / posteriors
   of a given parallel chain, given a .bin and .hdr file
*/
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <iomanip>
#include <archive.h>
#include <archive_entry.h>
#include "string_handler.h"
#include "version.h"
#include <boost/program_options.hpp>
#include "data.h"
//#include "interpol.h"
#include "diagnostics.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
namespace po = boost::program_options;

struct Proba_out{
	/*
	 * This is a structure that contains most of values that are saved
	 * for the likelihood/prior/posterior information for each chain
	*/
	MatrixXd logL;
	MatrixXd logPr;
	MatrixXd logP;
};

struct Proba_hdr{
	std::vector<std::string> header; // Any comment
	std::vector<std::string> labels; // The name of the variables
	long Nsamples_done, Nchains;
};

struct Ptemp_hdr{
	std::vector<std::string> header; // Any comment
	std::vector<std::string> labels; // The name of the variables
	VectorXd Tcoefs;
	long Nsamples_done;	
};

void showversion();
void usage(const po::options_description& desc);
bool file_exists (const std::string& name);
VectorXd reduce_data(VectorXd& samples, const int Samples_period, const int ind0, const int indmax);
Proba_hdr read_proba_header(const std::string file);
Ptemp_hdr read_ptemp_header(const std::string file);
Proba_out read_bin_proba_params(const std::string binfile, const long Nrows, const long Ncols);
Proba_out read_tar_gz_bin_proba_params(const std::string tarGzFile, const long Nrows, const long Ncols);

int main(int argc, char* argv[]){
/*
 * This is a procedure that read the binary files that contain information
 * about the parallel chains. See the structure Ptempering_out for further
 * information about the retrieved variables.
*/
	bool istar=false;
    const int Nchars=20;
    const int precision=10;
	int interpolation_factor;
	int cpt, lcpt, testval; // for the options
	int ind0, indmax, Samples_period; // The Samples_period defines out of all samples, how many we keep.
	long Nrows, Nchains;
	std::string file_in, file_proba,  file_proba_bin,  file_out, bootstrap_method;
	std::ifstream file;
	std::ostringstream strg;
	Proba_hdr proba_header;
	Ptemp_hdr ptemp_header;
	Proba_out proba;
	Diagnostics diags;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Print help message")
        ("version,v", "Print version information")
		("file-in,i", po::value<std::string>(&file_in)->required(), "input file name, without extension and qualifying names. Example: '4671239_A'. It is used to read [file-in]_stats.bin, [file-in]_stats.hdr and [file-in]_stats.bin [file-in]_parallel_tempering.hdr")
		("file-out,o", po::value<std::string>(&file_out)->required(), "output file name")
		("interpol-factor,f", po::value<int>(&interpolation_factor)->default_value(1000), "[Optional] Define the interpolation factor for f(beta) = P<P(D|M,I)>(beta). Recommended value if the default (1000)")
		("first-kept-element,S", po::value<int>(&ind0)->default_value(0), "[Optional] index of the first kept sample. All index below that will be discarded. By default, it will start at the element 0")
		("last-kept-element,L", po::value<int>(&indmax)->default_value(-1), "[Optional] index of the last kept sample. All index above that will be discarded. By default, it will take up to the last element (-1)")
		("periodicity,p", po::value<int>(&Samples_period)->default_value(7), "[Optional] sampling periodicity to perform the bootstrap with independent samples. The default value assumes a MH sampling rejection rate of 70%")
		("bootstrap-method,b", po::value<std::string>(&bootstrap_method)->default_value("none"), "Bootstrap method to compute uncertainties: 'none', 'standard' or 'block' (default: 'none')");
	po::variables_map vm;
	try {
		po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
		if (vm.count("help")) {
			std::cout << " This program allows you to calculate the evidence from a MCMC computation from TAMCMC." << std::endl;
			std::cout << " It also incorporate an uncertainty computation based on the bootstrap method. " << std::endl;
			usage(desc);
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
	file_proba=file_in + "_stat_criteria";
	const std::string file_proba_hdr=file_in + "_stat_criteria.hdr";
	const std::string file_ptempering_hdr=file_in + "_parallel_tempering.hdr"; // Use to get the Tcoefs vector

	std::cout << "  0. Configuration: " << std::endl;
	std::cout << "      - Input Statistical file: " << file_proba << " (.bin to .tar.gz automatic switch supported)"  << std::endl;
	std::cout << "      - Output ASCII file: " << file_out << std::endl;
	
	std::cout << "      Optional Arguments: " << std::endl;
	if (vm.count("interpol-factor")) {
		std::cout << "          - Interpolation factor for P<P(D|M,I)>(beta): " << interpolation_factor << std::endl;
	}
	if (vm.count("first-kept-element")) {
		ind0 = vm["first-kept-element"].as<int>();
		std::cout << "          - Index of the first kept sample: " << ind0 << std::endl;
	}
	if (vm.count("last-kept-element")) {
		indmax = vm["last-kept-element"].as<int>();
		std::cout << "          - Index of the last kept sample: " << indmax;
		if(indmax == -1){
			std::cout << " (until last sample)" << std::endl;
		} else{
			std::cout << std::endl;
		}
	}
	if (vm.count("periodicity")) {
		Samples_period = vm["periodicity"].as<int>();
		std::cout << "          - Samples_period: " << Samples_period << " ==> ";
		std::cout << "            Keep 1 sample every " << Samples_period << " samples" << std::endl;
	}
	if(Samples_period < 1){
		std::cout << "            Warning: The given periodicity value is smaller than 1... The program will use the default value instead (Period =1)" << std::endl;
		std::cout << "                     ===> All samples will be returned" << std::endl;
		Samples_period=1;
	}
	if(bootstrap_method != "standard" and bootstrap_method != "block" and bootstrap_method != "none"){
		std::cerr << "Error while setting the bootstrap argument" << std::endl;
		usage(desc);
	}
		
	
	std::cout << "  1. Reading the data files..." << std::endl;

	std::cout << "      # Reading Header files:" << std::endl;
	std::cout << "          - " <<  file_proba_hdr << std::endl;
	proba_header=read_proba_header(file_proba_hdr);
	std::cout << "          - " <<  file_ptempering_hdr << std::endl;
	ptemp_header=read_ptemp_header(file_ptempering_hdr);
	std::cout << "      - Temperature Law : Tcoefs=" << ptemp_header.Tcoefs[1] << "^(chain-1),  with Nchains=" <<ptemp_header.Tcoefs.size() << std::endl; 
	std::cout << "                    Gives Tcoefs = " << ptemp_header.Tcoefs.transpose() << std::endl; 


	std::cout << "      # Reading stats file..." << std::endl;// << file_proba << std::endl;
	file_proba_bin=file_proba + ".bin";
	if(!file_exists(file_proba_bin.c_str())){
		file_proba_bin=file_proba + ".tar.gz";
		if(!file_exists(file_proba_bin.c_str())){
			std::cerr << "Error: no binary or tar.gz file of this name found" << std::endl;
			std::cerr << "       check the usage or debug the program" << std::endl;
			std::cerr << "---" << std::endl;
		} else{
			std::cout << "         tar.gz  >> " << file_proba_bin << " file detected... proceeding in reading it..." << std::endl;
		}
		istar=true;
	} else{
		std::cout << "        bin   >> " << file_proba_bin <<  " file detected... proceeding in reading it..." << std::endl;
	}
	if(istar == false){
		proba=read_bin_proba_params(file_proba_bin,  proba_header.Nsamples_done, proba_header.Nchains);
	}
	else{
		proba=read_tar_gz_bin_proba_params(file_proba_bin, proba_header.Nsamples_done, proba_header.Nchains);
	}
	// Consistency checks on optional arguements
	if (indmax < 0){
		std::cout << std::endl;
		std::cout << "    Warning: Negative indmax provided ==> the last sample is going to be the last sample of the chain" << std::endl;
		std::cout << std::endl;
		indmax=proba_header.Nsamples_done;
	}
	if (indmax < ind0+Samples_period){
		std::cout << "    Error: The last index must be greater than the first sample + Periodicity " << std::endl;
		std::cout << "           The program will exit" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (Samples_period <=0){
		std::cout << "    Warning: Negative sampling period provided ==> sampling period will be fixed to 1 " << std::endl; 
		Samples_period=1;
	}
	// Apply filtering rules
	std::cout << "      # Reducing data by applying selection rules..." << std::endl;
	MatrixXd logL_chains;
	for (int n=0;n<ptemp_header.Tcoefs.size();n++){ 
		std::cout << "       >> chain: " << n << " ... " << std::endl;
		VectorXd tmp0=proba.logL.col(n);
		VectorXd tmp=reduce_data(tmp0, Samples_period, ind0, indmax);
		if (n==0){ // Declare a matrix of size determined by the size of the reduced data
			logL_chains.resize(tmp.size(), ptemp_header.Tcoefs.size());
		}
		logL_chains.col(n)=tmp;
	}
	/////// Write the parameters ////////	
	// Compute evidence
	Evidence_out evidence=diags.evidence_calc(ptemp_header.Tcoefs, logL_chains, interpolation_factor);
	const bool firstpass=1;
    diags.write_evidence(file_out, evidence, proba_header.Nsamples_done, firstpass);
	std::cout << "Evidence : " << evidence.evidence << std::endl;
	std::cout << "All done successfully" << std::endl;

return 0;
}

Proba_out read_bin_proba_params(const std::string binfile, const long Nrows, const long Ncols){
/*
 * Function that read the outputs file that contains inputs as writen in the file that contains
 * the information about the logLikelihood, logPrior and logPosterior. Return a structure of type Buffer_parallel_tempering
 * It requires as input:
 * 	- The name of the file with the data
*/

	double val_dbl=0;
	size_t size_dbl=sizeof(val_dbl);
	std::ifstream file;
	Proba_out proba_read;
	MatrixXd vals(Nrows, Ncols);

	proba_read.logL.resize(Nrows, Ncols); 
	proba_read.logPr.resize(Nrows, Ncols); 
	proba_read.logP.resize(Nrows, Ncols); 

	// Reading the file until the end;
	file.open(binfile.c_str(), std::ios::binary); // Open the binary file in read only
	if (file.is_open()){
		for (int i=0; i<Nrows; i++){ // we read the buffer
			for(int j=0; j<Ncols;j++){ // LogLikelihoods for all chains
				file.read(reinterpret_cast<char*>(&proba_read.logL(i,j)), size_dbl);
			}
			for(int j=0; j<Ncols;j++){ // LogPriors for all chains
				file.read(reinterpret_cast<char*>(&proba_read.logPr(i,j)), size_dbl);
			}
			for(int j=0; j<Ncols;j++){ // LogPosteriors for all chains
				file.read(reinterpret_cast<char*>(&proba_read.logP(i,j)), size_dbl);
			}	
		}
		file.close();
  	}
  	else {
		std::cout << " Unable to open file " << binfile  << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
return proba_read;
}

Proba_out read_tar_gz_bin_proba_params(const std::string tarGzFile, const long Nrows, const long Ncols) {
    double val_dbl = 0;
    size_t size_dbl = sizeof(val_dbl);
    Proba_out proba_read;
    proba_read.logL.resize(Nrows, Ncols);
    proba_read.logPr.resize(Nrows, Ncols);
    proba_read.logP.resize(Nrows, Ncols);

    struct archive_entry *entry;
    struct archive* ar = archive_read_new();
    archive_read_support_filter_gzip(ar);
    archive_read_support_format_tar(ar);
    int r = archive_read_open_filename(ar, tarGzFile.c_str(), 10240);
    if (r != ARCHIVE_OK) {
        std::cout << "Unable to open tar.gz file " << tarGzFile << std::endl;
        std::cout << "Check that the file exists" << std::endl;
        std::cout << "The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    }
    entry = archive_entry_new();
	r = archive_read_next_header(ar, &entry);
    if (r != ARCHIVE_OK) {
        std::cout << "Unable to read header from tar.gz file" << std::endl;
        std::cout << "The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::vector<char> buffer(archive_entry_size(entry));
    archive_read_data(ar, buffer.data(), buffer.size());
    std::istringstream iss(std::string(buffer.begin(), buffer.end()));
    for (int i = 0; i < Nrows; i++) {
        for (int j = 0; j < Ncols; j++) {
            iss.read(reinterpret_cast<char*>(&proba_read.logL(i, j)), size_dbl);
        }
        for (int j = 0; j < Ncols; j++) {
            iss.read(reinterpret_cast<char*>(&proba_read.logPr(i, j)), size_dbl);
        }
        for (int j = 0; j < Ncols; j++) {
            iss.read(reinterpret_cast<char*>(&proba_read.logP(i, j)), size_dbl);
        }
    }
    archive_read_close(ar);
    archive_read_free(ar);
    return proba_read;
}

Proba_hdr read_proba_header(const std::string file){

  int cpt=0;
  size_t pos=0;

  long varcase;
  std::string line;
  std::vector<std::string> keyword, varnames;
  std::ifstream myfile;

  Proba_hdr data;

  myfile.open(file.c_str());
  if (myfile.is_open()){
    while(getline(myfile,line)){
	keyword=strsplit(line, "=");

	varcase=3; // Defaut case is header case
	if(keyword[0]== "! labels"){
		varcase=0;
	}
	if(keyword[0]== "! Nchains"){
		varcase=1;
	}
	if(keyword[0]== "! Nsamples_done"){
		varcase=2;
	}
	switch(varcase){
		case 0: 
		  varnames=strsplit(keyword[1], " ");
		  for(int j=0; j<varnames.size(); j++) {data.labels.push_back(strtrim(varnames[j]));}
		  break;
		case 1: 
		  data.Nchains=0;
		  varnames=strsplit(keyword[1], " ");
		  std::stringstream(strtrim(varnames[0])) >> data.Nchains;
		  break;
		case 2: 
		  data.Nsamples_done=0;
		  varnames=strsplit(keyword[1], " ");
		  std::stringstream(strtrim(varnames[0])) >> data.Nsamples_done;
		  break;
		default:
		  data.header.push_back(strtrim(line)); // any comment or extra parameters are put in the header
	}
    }
    myfile.close();
  } else {
	std::cout << "    Unable to open the following parameter file: " << file << std::endl; 
	std::cout << "    The program will stop now" << std::endl;
	exit(EXIT_FAILURE);
  }

return data;
}


Ptemp_hdr read_ptemp_header(const std::string file){

  int cpt=0;
  size_t pos=0;

  long varcase;
  std::string line;
  std::vector<std::string> keyword, varnames;
  std::ifstream myfile;

  Ptemp_hdr data;

  myfile.open(file.c_str());
  if (myfile.is_open()){
    while(getline(myfile,line)){
		keyword=strsplit(line, "=");

		varcase=3; // Defaut case is header case
		if(keyword[0]== "! labels"){
			varcase=0;
		}
		if(keyword[0]== "! Tcoefs"){
			varcase=1;
		}
		if(keyword[0]== "! Nsamples_done"){
			varcase=2;
		}
		switch(varcase){
			case 0: 
				varnames=strsplit(keyword[1], " ");
				for(int j=0; j<varnames.size(); j++) {data.labels.push_back(strtrim(varnames[j]));}
				break;
			case 1: 
				data.Tcoefs=str_to_Xdarr(keyword[1], " ");
				break;
			case 2: 
				data.Nsamples_done=0;
				varnames=strsplit(keyword[1], " ");
				std::stringstream(strtrim(varnames[0])) >> data.Nsamples_done;
				break;
			default:
				data.header.push_back(strtrim(line)); // any comment or extra parameters are put in the header
		}
    }
    myfile.close();
  } else {
	std::cout << "    Unable to open the following parameter file: " << file << std::endl; 
	std::cout << "    The program will stop now" << std::endl;
	exit(EXIT_FAILURE);
  }

return data;
}

void showversion()
{
    std::cout << APP_NAME " - getstats tool - " APP_VERSION "\n built on " __DATE__ << std::endl;

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

    std::cerr << " You need to provide at least 3 arguments to that function. The available arguments are: " << std::endl;
	std::cerr << desc << std::endl;
	exit(EXIT_FAILURE);
}


bool file_exists(const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

VectorXd reduce_data(VectorXd& samples, const int Samples_period, const int ind0, const int indmax){
		int cpt, lcpt;
		VectorXd samples_out;
		if(Samples_period < 1){ // Case where we return all samples
			std::cout << "           Warning : Sample_period = " << Samples_period << " <= 1   ... Imposing Sample_period=1" << std::endl;
		} 
		samples_out.resize((indmax -ind0)/Samples_period + 1);
		// ---- Selecting only 1 sample out every Sample_period samples -----
		if (Samples_period != 1){ 
			std::cout << "           Operation is of the type: ";
			std::cout << "           Outputs[" << ind0 << "]  --> New_Outputs[" << ind0 + Samples_period-1 << "]  "  << std::endl;		
			cpt=ind0; lcpt=0;
			while(cpt<samples.size()){
				samples_out[lcpt]=samples[cpt];
				cpt=cpt+Samples_period;
				lcpt=lcpt+1;
			}
			samples_out.conservativeResize(lcpt-1);
			std::cout << "                       New size : " << (indmax -ind0)/Samples_period << std::endl;
		} else{
			samples_out=samples;
			std::cout << "                       size kept at: " << samples.size() << std::endl;
		}
	return samples_out;
}