/*
   A small function that extract one parameter from the binary output
   That parameter is saved into a output file
*/

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>
#include "diagnostics.h"
#include "data.h"
#include "string_handler.h"
#include "version.h"

void showversion();
int options(int argc, char* argv[]);
void usage(int argc, char* argv[]);

int main(int argc, char* argv[]){

		bool replicate_cte;
		int ind0, ind_param, ind_var, ind_cons, ind_chain, Nsamples, Samples_period, Newsize, val; // The Samples_period defines out of all samples, how many we keep.
		long cpt, lcpt;
		Eigen::MatrixXd data_array, data_out;
		std::string rootname, filename_params, filename_params_hdr, dir_out, file;
		std::ofstream fileout_stream;
		std::string cpath=getcwd(NULL, 0);
		
		Diagnostics diags;
		Params_hdr hdr;

		val=options(argc, argv);
		rootname=argv[1]; 
		std::istringstream(argv[2]) >> ind_chain; 
		dir_out=argv[3];
		std::istringstream(argv[4]) >> ind0;
		std::istringstream(argv[5]) >> Samples_period;
		if(Samples_period < 1){
			std::cout << "Warning: The given periodicity value is smaller than 1... The program will use the default value instead (Period =1)" << std::endl;
			std::cout << "          ===> All samples will be returned" << std::endl;
			Samples_period=1;
		}
		if (val == 11){
			std::istringstream(argv[6]) >> replicate_cte;
		} else{
			std::cout << " Binary variable set to default: 0 ==> Constant values are written once" << std::endl;
			replicate_cte=0;
		}

		filename_params=rootname + "_chain-" + int_to_str(ind_chain) + ".bin";
		filename_params_hdr=rootname + ".hdr";
		
		std::cout << "  0. Configuration: " << std::endl;
		std::cout << "      - Binary file: " << filename_params << std::endl;
		std::cout << "      - Header file: " << filename_params_hdr << std::endl;
		std::cout << "      - Output directory: " << dir_out << std::endl;
		std::cout << "      - Index of the first kept sample: " << ind0 << std::endl;
		std::cout << "      - Samples_period: " << Samples_period << " ==> ";
		if(Samples_period >1){
			std::cout << " Keep 1 sample every " << Samples_period << " samples" << std::endl;
		} else{
			std::cout << " Keep all samples" << std::endl;
		}
		std::cout << "      - IDL compatibility: " << replicate_cte << std::endl;

		std::cout << "  1. Reading the binary output files..." << std::endl;
		hdr=diags.read_params_header(filename_params_hdr); // Get the metadata from the header file 
		data_array=diags.read_bin_matrix_dbl(filename_params, hdr.Nvars, hdr.Nsamples, "dbl");
		
		// ---- Selecting only 1 sample out every Sample_period samples -----
		std::cout << "  2. Applying the selection rule..." << std::endl;
		if(Samples_period <= 1){ // Case where we return all samples
			data_out = data_array.block(ind0, 0, data_array.rows()-ind0, data_array.cols()) ;
		} else{ // Case where we keep one row out of Sample_period;
			Newsize=(data_array.rows() -ind0)/Samples_period;
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
    	for(int ind_param=0; ind_param < hdr.Nvars+hdr.Ncons; ind_param++){
    		file=diags.formated_int_to_str(ind_param);
    		fileout_stream.open((dir_out + file + ".ASCII").c_str());
    		if(fileout_stream.is_open()){
				if(hdr.relax[ind_param] == 1){
					fileout_stream << "! variable_name= " << hdr.variable_names[ind_var] << std::endl;
					std::cout << "    - File: " << file + ".ASCII" << "   " << std::endl;
					fileout_stream << std::setprecision(12) << data_out.col(ind_var) << std::endl;
					ind_var=ind_var+1;
    			} else{
					fileout_stream << "! constant_name= " << hdr.constant_names[ind_cons] << std::endl;
					std::cout << "    - File: " << file + ".ASCII" << "  (Constant value)" << std::endl;
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
				std::cout << " Unable to open the binary data file " << dir_out.c_str() + file + ".ASCII" << std::endl;	
				std::cout << " Check that the full path exists" << std::endl;
				std::cout << " The program will exit now" << std::endl;
				exit(EXIT_FAILURE);
    		}
    		fileout_stream.close();
		}
		
		// --- Write a single line with the plength vector on a file ---
		file=dir_out + "plength.txt";
		fileout_stream.open(file.c_str());
		if(fileout_stream.is_open()){
				fileout_stream << hdr.plength << std::endl;
   		} else{
			std::cout << " Unable to open the binary data file " << dir_out.c_str() + file + ".ASCII" << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
    	}
    	fileout_stream.close();
		
	std::cout << "All done" << std::endl;
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

int options(int argc, char* argv[]){

	std::string arg1, arg2;
	int val;
	
	val=-1;
	arg1="";
	
	if(argc == 2){
		arg1=argv[1];
		if(arg1 == "version"){
			 val=0;
		} 
	}
	if(argc == 6){
		val=10;
	}
	if(argc == 7){
		val=11;
	}

	if (val == -1){ // Error code
		usage(argc, argv);
	} 
	if (val == 0){ // Version code
		showversion(); 
		exit(EXIT_SUCCESS);
	}
	if (val > 0 ){
		return val; // Execution code val
	} else{
		return -1; // Default value is to return an error code
	}
}

void usage(int argc, char* argv[]){

			std::cout << " You need to provide at least 5 arguments to that function. The available arguments are: " << std::endl;
			std::cout << "     [1] The root name (along with the full path) of the binary/header file containing the parameters ('[out_root_name]_chain-*' file)" << std::endl;
			std::cout << "     [2] The chain index (e.g. 0 for the coldest chain)" << std::endl;
			std::cout << "     [3] The output directory (must already exist)" << std::endl;
			std::cout << "     [4] Index of the first element which we keep. All index below that will be discarded" << std::endl;
			std::cout << "     [5] The Periodicity at which we keep samples. If <1 then all samples are returned" << std::endl;	
			std::cout << "     [6] [Optional] Binary variable. If 1, then constant values will be written Nsamples/Nperiod times in the ASCII file (For IDL code compatibility)" << std::endl;
			std::cout << "                                     If 0, then constant values will be written once in the ASCII file" << std::endl;
			std::cout << "                                     Default is 0" << std::endl;
			std::cout << " Call sequence: " << std::endl;
			std::cout << "     " << argv[0] << " [rootname] [chain index] [output directory] [First kept element] [Periodicity] [IDL compatibility]" << std::endl;
			exit(EXIT_FAILURE);

}
