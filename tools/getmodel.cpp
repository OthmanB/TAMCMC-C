/*
   A small function that extract one parameter from the binary output
   That parameter is save into a output file
*/

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>
//#include <filesystem> // ONLY FOR C++17
#include <cstdio>
#include <ios>
#include <fstream>
#include <cstring>
#include <cerrno>

#include "models.h"
#include "model_def.h"
#include "data.h"
#include "string_handler.h"
#include "linspace.h"
#include "config.h"
#include "version.h"

void showversion();
int options(int argc, char* argv[], int Nmaxlines);
void usage(int argc, char* argv[], int Nmaxlines);
int check_retrocompatibility(VectorXi plength, std::string modelname);
VectorXd adapt2new_MSGlobal(const VectorXi plength, const VectorXd params, const double c0);
void mvfile(std::string file_in, std::string file_out);
//inline bool exists_test (const std::string& name);

Data Data_Nd2Data(Data_Nd dat);

int main(int argc, char* argv[]){

		bool verbose=0;
		int i, Nmaxlines, Nmodels, modelname_switch, testval, status_compat, pos_option_file, pos_option_type;
		VectorXi plength;
		VectorXd arr0, model0, data_range;
		MatrixXd data_out;
		MatrixXd models;
		Data_Nd data;
		Data data_model;
		Config cfg;
		std::string filename_data, filename_params, dirout, fileout, modelname;
		std::string char0, line0;
		std::ifstream cfg_session;
		std::ofstream fileout_stream;
		std::vector<std::string> tmp;
		Model_def model_list;
		Data_Basic listoutputs;

		// Set the current path variables
		std::string cpath=getcwd(NULL, 0);
		Nmaxlines=5; // Maximum of models that can be provided (and therefore generated)
		testval=options(argc, argv, Nmaxlines);
		std::cout << testval << std::endl;
		filename_data=argv[1];
		filename_params=argv[2];
		modelname=argv[3]; 
		if(testval == 10 || testval == 13){
			std::cout << " Output directory / Output file not provided " << std::endl;
			std::cout << "       - Generating the output files within the current directory" << std::endl;
			dirout="";
			fileout="output_model.ascii";
			std::cout <<"       - Default filename used: " << fileout << std::endl;
			if(testval == 13){
				pos_option_type=4;
			}
		} else{ //Then we expect the first option, followed by the second one
			if(testval == 11 || testval == 12){
				pos_option_file=4;
				pos_option_type=5;
			}
			if(testval == 15){ // Case where the data_type argument is before the output dir
				pos_option_file=5;
				pos_option_type=4;
			}
			//exit(EXIT_SUCCESS);
			// Identifies if the user provided a directory (end by a '/') or a file (does not end by a '/')
			tmp=strsplit(argv[pos_option_file], "/");
			if(tmp[tmp.size()-1] == ""){ // In that case, we have a directory as input. Not a file. The default file for outputs will be provided
				std::cout <<" Directory was provided as output... using it for outputs of the ascii and the param file" << std::endl;
				dirout=argv[pos_option_file];
				fileout="output_model.ascii";
				std::cout <<"       Default filename used for the ASCII file: " << fileout << std::endl;
			} else{
				dirout="/";
				if(tmp[0] != "."){
					for(int i=0; i<tmp.size()-1;i++){
						dirout=dirout + tmp[i] + "/";
					}
				} // Otherwise, the user specifically asked the current dir
				fileout=argv[pos_option_file];
				std::cout << " Directory and filename provided for outputs" << std::endl;
				if (dirout !="/"){
					std::cout << "       - Identified directory: " << dirout << std::endl;
				} else{
					std::cout << "       - Identified directory: " << "Current Directory" << std::endl;
				}
				std::cout << "       - Provided filename used for the ASCII file :" << tmp[tmp.size()-1] << std::endl;
			}
		}
		std::cout << "  0. Configuration: " << std::endl;
		std::cout << "      - Model name: " << modelname << std::endl;		
		if(testval == 10 || testval == 11){
			std::cout << "      - Data file: " << filename_data << std::endl;
		} else{
			tmp=strsplit(argv[pos_option_type], "="); // Dealing with case where data_type option is requested, but IS GIVEN BEFORE output dir provided
			if(tmp[1] == "range"){
				data_range=str_to_Xdarr(filename_data, ",");// filename_data contains 3 values the string in that scenario
				if (data_range.size() !=3){
					std::cout << " Error: A range was provided but it seems that it does not contain 3 elements" << std::endl;
					std::cout << "        The expected syntax is: data_range=x_min,x_max,resolution (no spaces)" << std::endl;
					std::cout << "        Inputs for the argument: "<< filename_data << std::endl;
					exit(EXIT_FAILURE);
				}
				std::cout << "      - Data range: " << data_range.transpose() << std::endl;
				filename_data="";
			} else{
				std::cout << "      - Data file: " << filename_data << std::endl;
			}
		}
		std::cout << "      - Parameters file: " << filename_params << std::endl;
		std::cout << "      - Output dir: " << dirout << std::endl;
		
		if(filename_data != ""){
			std::cout << "  1. Reading the data file..." << std::endl;
			verbose=0;
			data=cfg.read_data_ascii_Ncols(filename_data, " \t", verbose); // the data in a matricial form				
			data_model=Data_Nd2Data(data);
		} else{
			std::cout << "  1. Generating x-axis data according to the provided data_range values..." << std::endl;
			data_model.x=linspace(data_range[0], data_range[1], std::ceil(1 + (data_range[1]-data_range[0])/data_range[2])); // convert the resolution into a number of data points
			data_model.xrange.resize(2); data_model.y.resize(1); data_model.sigma_y.resize(1);
			data_model.xrange[0]=data_model.x.minCoeff(); // Contains the min and max of x... in practice, it is used to limit the data range if requested by the cfg file (e.g. freq_range option in the .MCMC)
			data_model.xrange[1]=data_model.x.maxCoeff(); // Contains the min and max of x... in practice, it is used to limit the data range if requested by the cfg file (e.g. freq_range option in the .MCMC)
			data_model.y[0]=-1; data_model.sigma_y[0]=-1;
			data_model.Nx=data_model.x.size(); // Ny is not checked but should be as long as x
			data_model.xlabel="x-axis"; // label for x-axis. In ASCII file, the axis labels are identified by the symbol ! at the begining of the line
			data_model.ylabel="NOT USED";
			data_model.xunit="Unknown";
			data_model.yunit="NOT USED";
			// To ensure that we write the x-axis in the output file, we also set a minimalistic data_Nd raw structure (without y-axis, header, units and labels)
			data.data.resize(data_model.Nx,1);
			data.data.col(0)=data_model.x;
		}
		// Getting the list of models that are currently implemented
		std::string file_list=cpath + "/models_ctrl.list";
		listoutputs=cfg.read_listfiles(file_list, 1);	 // The file must be in the same folder as the main getmodel compiled file
		modelname_switch=cfg.convert_model_fct_name_to_switch(modelname, listoutputs); // look for the case number that is going to be used in call_model

		std::cout << "  2. Reading the file with the parameters of the model and computing model(s)..." << std::endl;
		cfg_session.open(filename_params.c_str());
		if (cfg_session.is_open()) {
			char0="#";
			std::getline(cfg_session, line0);	
			while(!cfg_session.eof() && char0 == "#"){ // Jump comments lines in the header
				line0=strtrim(line0);
				char0=strtrim(line0.substr(0, 1));
				if (char0 == "#"){
					std::getline(cfg_session, line0);
				}
			}
			// After all the comments, the first line must contain plength
			plength=str_to_Xiarr(line0, " \t"); // The separator is either a space of a tabulation
			status_compat=check_retrocompatibility(plength, modelname); // Verify the consistency of the vector with earlier version of the code
			// All the following lines must correspond to a vector of inputs. 1 line <==> 1 model. Maximum of Nmaxlines lines permitted
			i=0;
			std::getline(cfg_session, line0);
			while(!cfg_session.eof() && i < Nmaxlines){ 
				//std::getline(cfg_session, line0);
				line0=strtrim(line0);
				char0=strtrim(line0.substr(0, 1));
				
				arr0=str_to_Xdarr(line0, " \t"); // Read the input line
				if(status_compat == 1){ // Attempting adjustement of the size of the vectors
					std::cout << " >>> Incompatibility in vector sizes is attempted to be solved by adding:" << std::endl;
					std::cout << "		(1) c0=20 at the end of the parameter vector" << std::endl;
					std::cout << "		(2) conversion of inclination / splitting into nus.sin(i), nus.cos(i)" << std::endl;
					std::cout << "     This because since version 1.3.0 the Lorentzian truncation parameter is an adjustable variable and the splitting projecions are fitted" << std::endl;
					plength.conservativeResize(plength.size() +1);
					plength[plength.size()-1]=1; // adding c0
					arr0=adapt2new_MSGlobal(plength, arr0, 20.); // last value is c0
				}	
				model0=model_list.call_model_explicit(&data_model, plength, arr0, modelname_switch, true);
				if(i==0){// Initialisation of the Matrix of parameters
					models.resize(Nmaxlines, model0.size());
				} 
				models.row(i) = model0; 
				std::getline(cfg_session, line0);
				// Taking care of the params.model file that is created since version 1.61
				mvfile(cpath + "/params.model", dirout + "/params_" + int_to_str(i) + ".model");
				std::cout << "    " << cpath + "/params.model" << " --> " << dirout + "/params_" + int_to_str(i) + ".model" << std::endl;
				i=i+1;
			}	
			Nmodels=i;	
			cfg_session.close();
 	 	} else {
   			std::cout << "Unable to open the file: " << filename_params << std::endl;
   			std::cout << "Check that the file exist and that the path is correct" << std::endl;
   			std::cout << "Cannot proceed" << std::endl;
   			std::cout << "The program will exit now" << std::endl;
   			exit(EXIT_FAILURE);
		}
		std::cout << "  3. Writing outputs into the output file in a simple matricial format... " << std::endl;
		data_out.resize(data.data.rows(), data.data.cols() + Nmodels) ; //We take the initial data and the models as additional columns
		data_out.block(0, 0, data.data.rows(), data.data.cols())=data.data; // Put the data first
		for( i=0; i<Nmodels; i++){
			data_out.col(data.data.cols() + i)=models.row(i);
		}
    	fileout_stream.open(fileout.c_str());
    	if(fileout_stream.is_open()){
			fileout_stream << std::setprecision(12) << data_out << std::endl;
    	} else{
			std::cout << " Unable to open the binary data file " << fileout.c_str() << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
    	}
    	fileout_stream.close();
	
	std::cout << "Output model file successfully written" << std::endl;
}

Data Data_Nd2Data(Data_Nd dat){

		Data dat_out;
		dat_out.x=dat.data.col(0);
		dat_out.y=dat.data.col(1);
		if (dat.data.cols() == 3){
			dat_out.sigma_y=dat.data.col(2);
		}
		dat_out.Nx=dat.data.rows();
	return dat_out;
}

void showversion()
{
    std::cout << APP_NAME " - getmodel tool - " APP_VERSION "\n built on " __DATE__ << std::endl;

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

int options(int argc, char* argv[], int Nmaxlines){

	std::string arg1;
	std::vector<std::string> tmp_str;
	int val;
	
	val=-1;
	arg1="";
	 // Deal with the case that any of the argument happens to be requesting for parameters only (no spectrum)
	 // This case is useful to get for example only the l=1 mixed modes, which are a solution of fitting parameters
	if(argc == 2){
		arg1=argv[1];
		if(arg1 == "version"){
			 val=0;
		} 
	}
	if(argc == 4){
		val=10;
	}
	if(argc == 5){
		val=11; // Case where we had ONE optional parameter
		tmp_str=strsplit(argv[4], "="); // Dealing with case where data_type option is requested, but NO output dir provided (no first option)
		if (tmp_str[0] == "data_type"){
			val=13;
			//if(tmp_str[1] == "range"){
			//	val=13;
			//}
		}
	}
	if(argc == 6){ // Case where we had TWO optional parameters 
		val=12; // The expected order for the argument (directory first, then data_type)
		tmp_str=strsplit(argv[4], "="); // Dealing with case where data_type option is requested, but IS GIVEN BEFORE output dir provided
		if (tmp_str[0] == "data_type"){
			val=15;
			//if(tmp_str[1] == "range"){
			//	val=15;
			//}
		}

	}

	if (val == -1){ // Error code
		usage(argc, argv, Nmaxlines);
	} 
	if (val == 0){ // Version code
		showversion(); 
		exit(EXIT_SUCCESS);
	}
	if (val > 0 ){
		return val; // Execution code val
	} else{
		usage(argc, argv, Nmaxlines);
		return -1; // Default value is to return an error code
	}
}

void usage(int argc, char* argv[], int Nmaxlines){

		if( argc < 4){
			std::cout << " You need to provide at least three argument to that function. The available arguments are: " << std::endl;
			std::cout << "     [1] The data filename or range. By default a file is expected and it should be in the same format as those provided to TAMCMC (*.data file)" << std::endl;
			std::cout << "     [2] The filename for the parameters to read. After comments ('#'), this files contains on row(1) plength and on row(2:2+Nmaxlines) the model parameters. Maximum number of models is Nmaxlines=" << Nmaxlines << std::endl;
			std::cout << "         THIS FILE IS NOT A DIRECT OUTPUT OF THE TAMCMC program. You must make it yourself, using e.g. the mean of the parameters (that can be extracted the bin2txt program)" <<std::endl;
			std::cout << "     [3] The model name among the family of MS_Global models (e.g. model_MS_Global_a1etaa3_HarveyLike)" << std::endl;
			std::cout << "     [4] [Optional] the output directory for files. If not given, then the output directory is the current one" << std::endl;
			std::cout << "     [5] [Optional argument] data_type=[string]. [string] can be:" << std::endl;
			std::cout << "               - 'range' in that case arg[1] must be 3 values (xmin,xmax,resolution) without spaces" << std::endl; 
			std::cout << "               - 'file' in that case arg[1] must be the full path to a data file " << std::endl;
			std::cout << "         Remark concerning data_type: If data_type is 'range', the output file will not have a copy of the observational data" << std::endl;
			std::cout << "         Remark concerning the options: There is no order requirement. However, there is no consistency check between the data_type and the arg1 content" << std::endl;
			std::cout << " Call sequence: " << std::endl;
			std::cout << "     " << argv[0] << " [data filename] [parameter filename]  [model name] [output filename] [Options]" << std::endl;
			exit(EXIT_FAILURE);
		}
	
}

// Small function that verify that plength is of the correct length knowing the modelname. See the full list of referenced modelnames (as listed in config.cpp)
// This function needs a thorough check when problems are detected... 
int check_retrocompatibility(VectorXi plength, std::string modelname){

	int Nplength_expected=-1; // By defaut, consider that plength is not compatible with the current program
	int status=2; // DEFAULT IS EVERYTHING OK (EXPERIMENTAL... IF FAILS NEED TO BE TO -1)
	if (modelname == "model_MS_Global_a1etaa3_HarveyLike" || modelname == "model_MS_Global_a1etaa3_Harvey1985" || modelname == "model_MS_Global_a1l_etaa3_HarveyLike" ||
	   modelname == "model_MS_Global_a1n_etaa3_HarveyLike" || modelname == "model_MS_Global_a1nl_etaa3_HarveyLike"){
		Nplength_expected = 11;
		if(Nplength_expected != plength.size()){
			std::cout << "   >> Structure of the parameters incorrect. Incompatibility detected" << std::endl;
			std::cout << "   >> Detected size for plength: "<< plength.size() << " ...Expected size:"<< Nplength_expected << " for model '" << modelname << "'" << std::endl;
			if(plength.size() == 10){ // Case where we might know how to get back to our feet
				status=1;
                //std::cout << "Status set to 1" << std::endl;
			} else{
                //std::cout << "Status set to -1" << std::endl;
				status=-1; // There is problem that requires a carefull check/debuging
			}
        } else{
            //std::cout << "Status set to 2" << std::endl;
            status=2;
        }
	}
    if (modelname == "model_MS_Global_a1etaa3_HarveyLike_Classic" ){ // This new function appears in version 1.3.2 so it is ok
        //std::cout << "Status set to 2 (Classic model)" << std::endl;
        status=2;
    }
    if (modelname == "model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v2" ){ // This new function appears in version 1.5 so it is ok
        //std::cout << "Status set to 2 (Classic model)" << std::endl;
        status=2;
    }
    /*if (modelname =="model_MS_local_basic"){
    	std::cout << "   >> The model 'model_MS_local_basic' is not handled by getmodel yet" << std::endl;
    	std::cout << "      Currently working on its support... please wait that it gets released" << std::endl;
    	exit(EXIT_FAILURE);
    }
    */
	if(status >1){
		std::cout << "   >> Compatibility/Consistency test with earlier version passed..." << std::endl;
	}
	if(status == -1){
		std::cout << "   >> Check \n (1) The compatibility / consistency of your model names with the main program. \n (2) The version of this code matches the version of the main TAMCMC code" << std::endl;
		exit(EXIT_FAILURE);
	}

return status;
}

// Perform the vector adjustements required in case the MS_Global model use the standard predating version 1.3.0
VectorXd adapt2new_MSGlobal(const VectorXi plength, const VectorXd params, const double c0){

	const long double pi = 3.141592653589793238L;
	const int posa1=plength[0]+plength[1]+plength[2]+plength[3]+plength[4]+plength[5]; // old variable a1
	const int posinc=plength.size()-1; // old variable inclination
	const int posa1cosi=plength[0]+plength[1]+plength[2]+plength[3]+plength[4]+plength[5] + 3; // Position of the new variable that are used instead of a1-inclination
	const int posa1sini=plength[0]+plength[1]+plength[2]+plength[3]+plength[4]+plength[5] + 4; // Position of the new variable that are used instead of a1-inclination

	double inclination, a1;	
	VectorXd newparams=params;

	inclination=params[posinc]; // Last elements of the vector was the inclination in version earlier than 1.3.0
	a1=params[posa1]; ; //splitting a1

	newparams[posa1cosi]=sqrt(a1)*cos(inclination*pi/180.);
	newparams[posa1sini]=sqrt(a1)*sin(inclination*pi/180.);

	newparams.conservativeResize(params.size()+1);
	newparams[newparams.size()-1]=c0; //Put a default c0 value at the end of the new vector of parameters

	return newparams;
}

/*
   ONLY FOR C++17
*/
/*void mvfile(std::string file_in, std::string file_out) {
  try {
    std::filesystem::rename(file_in, file_out);
  } catch (std::filesystem::filesystem_error& e) {
    std::cout << e.what() << '\n';
  }
}
*/
// Another version that should work always
void mvfile(std::string file_in, std::string file_out) {
  if(std::rename(file_in.c_str(), file_out.c_str()) < 0) {
    std::cout << strerror(errno) << '\n' << std::endl;
  }
}

/*
void mvfile(std::string file_in, std::string file_out) {
  std::ifstream in(file_in.c_str(), std::ios::in | std::ios::binary);
  std::ofstream out(file_out.c_str(), std::ios::out | std::ios::binary);
  out << in.rdbuf();
  std::remove(file_in.c_str());
}
*/

/*
inline bool exists_test1 (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}
*/