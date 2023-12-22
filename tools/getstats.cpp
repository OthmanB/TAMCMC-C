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

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
namespace po = boost::program_options;

void showversion();
void usage(const po::options_description& desc);
bool file_exists (const std::string& name);

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

Proba_hdr read_proba_header(const std::string file);
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
	int readchain; // index of the chain to be read
	int cpt, lcpt, testval; // for the options
	int ind0, indmax, Samples_period; // The Samples_period defines out of all samples, how many we keep.
	long Nrows, Nchains;
	std::string file_proba, file_proba_bin, file_proba_hdr, file_out;
	std::ifstream file;
	std::ostringstream strg;
	std::ofstream outfile_proba;
	Proba_hdr hdr;
	Proba_out proba;
	
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Print help message")
        ("version,v", "Print version information")
		("file-stats,i", po::value<std::string>(&file_proba)->required(), "input file name, without extension. Example: '4671239_A_stat_criteria'")
		("chain-index,c", po::value<int>(&readchain)->required(), "index of the chain to be read")
		("file-out,o", po::value<std::string>(&file_out)->required(), "output file name")
		("first-kept-element,S", po::value<int>(&ind0)->default_value(0), "[Optional] index of the first kept sample. All index below that will be discarded")
		("last-kept-element,L", po::value<int>(&indmax)->default_value(-1), "[Optional] index of the last kept sample. All index above that will be discarded")
		("periodicity,p", po::value<int>(&Samples_period)->default_value(1), "[Optional] sampling periodicity");

	po::variables_map vm;
	try {
		po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
		if (vm.count("help")) {
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
	file_proba = vm["file-stats"].as<std::string>();
	file_out = vm["file-out"].as<std::string>();
	readchain = vm["chain-index"].as<int>();

	std::cout << "  0. Configuration: " << std::endl;
	std::cout << "      - Header file: " << file_proba_hdr << std::endl;        
	std::cout << "      - Binary file: " << file_proba << " (.bin or .tar.gz supported)"  << std::endl;
	std::cout << "      - Reading chain: " << readchain << std::endl;
	std::cout << "      - Output ASCII file: " << file_out << std::endl;

	if (vm.count("first-kept-element")) {
		ind0 = vm["first-kept-element"].as<int>();
		std::cout << "      Optional Arguments: " << std::endl;
		std::cout << "          - Index of the first kept sample: " << ind0 << std::endl;
	}
	if (vm.count("last-kept-element")) {
		indmax = vm["last-kept-element"].as<int>();
		std::cout << "          - Index of the last kept sample: " << indmax << std::endl;
	}
	if (vm.count("periodicity")) {
		Samples_period = vm["periodicity"].as<int>();
		std::cout << "          - Samples_period: " << Samples_period << " ==> ";
		std::cout << " Keep 1 sample every " << Samples_period << " samples" << std::endl;
	}
	if(Samples_period < 1){
		std::cout << "Warning: The given periodicity value is smaller than 1... The program will use the default value instead (Period =1)" << std::endl;
		std::cout << "          ===> All samples will be returned" << std::endl;
		Samples_period=1;
	}
		
	file_proba_hdr=file_proba + ".hdr";
	file_proba_bin=file_proba + ".bin";
	std::cout << "  1. Reading the data file..." << std::endl;

	std::cout << "      # Reading Header file:" << file_proba_hdr << std::endl;
	hdr=read_proba_header(file_proba_hdr);

	std::cout << "# labels= ";
	for (int i=0; i<hdr.labels.size(); i++){
		std::cout << hdr.labels[i] << "    ";
	}
	std::cout << std::endl;
	std::cout << "      # Nsamples_done=" << hdr.Nsamples_done << std::endl;

	std::cout << "      # Nchains=" << hdr.Nchains << std::endl;

	std::cout << "      # Reading Filename:" << std::endl;// << file_proba << std::endl;
	if(!file_exists(file_proba_bin.c_str())){
		file_proba_bin=file_proba + ".tar.gz";
		if(!file_exists(file_proba_bin.c_str())){
			std::cerr << "Error: no binary or tar.gz file of this name found" << std::endl;
			std::cerr << "       check the usage or debug the program" << std::endl;
			std::cerr << "---" << std::endl;
		} else{
			std::cout << "   tar.gz  >> " << file_proba_bin << " file detected... proceeding in reading it..." << std::endl;
		}
		istar=true;
	} else{
		std::cout << "  bin   >> " << file_proba_bin <<  " file detected... proceeding in reading it..." << std::endl;
	}
	if(istar == false){
		proba=read_bin_proba_params(file_proba_bin, hdr.Nsamples_done, hdr.Nchains);
	}
	else{
		proba=read_tar_gz_bin_proba_params(file_proba_bin, hdr.Nsamples_done, hdr.Nchains);
	}
	// Consistency checks on optional arguements
	if (indmax < 0){
		std::cout << std::endl;
		std::cout << "    Warning: Negative indmax provided ==> the last sample is going to be the last sample of the chain" << std::endl;
		std::cout << std::endl;
		indmax=hdr.Nsamples_done;
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
	/////// Write the parameters ////////	
	outfile_proba.open(file_out.c_str()); 
		 if (outfile_proba.is_open()){
				outfile_proba << "# This File contains the likelihood/priors/posteriors of a MCMC process \n";
				// ----			
				outfile_proba << "# Header file:" << file_proba_hdr << "\n";
				outfile_proba << "# Stats file:" << file_proba_bin << "\n";
				outfile_proba << "# Nsamples= " <<  hdr.Nsamples_done << "\n";
				outfile_proba << "# Nchains= " <<hdr.Nchains << "\n";
				outfile_proba << "# readchain= " << readchain << "\n";
				outfile_proba << "# ";
				outfile_proba << "     logLikelihood  " << "   logPrior  " << "  logPosterior     " << "\n";
				outfile_proba.flush(); // Explicitly specify to flush the header into the disk
				// Writing Data 
				if(Samples_period <= 1){ // Case where we return all samples
					for(int i=0; i<hdr.Nsamples_done; i++){	
						strg.str(std::string());	
						strg << std::setw(Nchars) << std::setprecision(precision) << proba.logL(i,readchain)  << std::setw(Nchars) << std::setprecision(precision) << proba.logPr(i,readchain)  << std::setw(Nchars) << std::setprecision(precision) << proba.logP(i,readchain) << "\n";
    					outfile_proba << strg.str().c_str();
					}
				} else{
					// ---- Selecting only 1 sample out every Sample_period samples -----
					std::cout << "  2. Applying the selection rule..." << std::endl;
					std::cout << "     Operation is of the type: ";
					std::cout << "     Outputs[" << ind0 << "]  --> New_Outputs[" << ind0 + Samples_period-1 << "]  "  << std::endl;		
					std::cout << "     ..." << std::endl;
					cpt=ind0; lcpt=0;
					while(cpt<hdr.Nsamples_done){
						strg.str(std::string());	
						strg << std::setw(Nchars) << std::setprecision(precision) << proba.logL(cpt,readchain) << std::setw(Nchars) << std::setprecision(precision) << proba.logPr(cpt,readchain)  <<std::setw(Nchars) << std::setprecision(precision) << proba.logP(cpt,readchain) << "\n";
    					outfile_proba << strg.str().c_str();
						cpt=cpt+Samples_period;
						lcpt=lcpt+1;
					}
					std::cout << "                       New size : " << (indmax -ind0)/Samples_period << std::endl;
				}
			outfile_proba.flush(); // Explicitly specify to flush the data into the disk
			strg.str(std::string());
			outfile_proba.close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << file_out << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

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

	std::cout << "Ncols = " << Ncols << std::endl;
	std::cout << "Nrows = " << Nrows << std::endl;
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
		//std::cout << "-----------" << std::endl;
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