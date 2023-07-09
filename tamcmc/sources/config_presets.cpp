/*
 * Config_presets.cpp
 *
 * Header file that contains all kind of class/structures
 * used to tune the configuration using preconfigured 
 * step-by-step analysis
 * 
 *  Created on: 29 May 2016
 *      Author: obenomar
 */

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <sys/stat.h>
#include "config_presets.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

Config_presets::Config_presets(std::string cfg_file_presets_in, Config *cfg0){

	cfg_file_presets=cfg_file_presets_in;

	// ---- Read the preset configuration files and perform the setup ----
	read_cfg_file_presets();
	// ------------------------------------------------------------

	current_id_ind=0; // index pointing to the current object id that we have to process
	current_process_ind=0; // index pointing to the current processing step that we have to execute
	current_slice_ind=0; // The index pointing to the current subdataset range (e.g. frequency range for io_local models)
	Nt_learn0=cfg0->MALA.Nt_learn; // We need this to be sure that this configuration is static when iterating on id_ind and process_ind
	dN_mixing0=cfg0->MALA.dN_mixing; // We need this to be sure that this configuration is static when iterating on id_ind and process_ind

}


void Config_presets::apply_presets(Config *cfg, const bool verbose){
/*
 * Method that is used to superseed parameters of the config_default.cfg
 * by applying the configuration as it is defined into the config_presets.cfg
 * 
*/
	bool passed;
    std::string process_id;

	// ---------- Setting the directories ----------
	cfg->outputs.dir_out=cfg_out_dir + "/" + table_ids[current_id_ind].at(0) + "/outputs/";
	cfg->outputs.restore_dir=cfg_out_dir + "/" + table_ids[current_id_ind].at(0) + "/restore/"; 
	cfg->diags.output_dir=cfg_out_dir + "/" + table_ids[current_id_ind].at(0) + "/diags/";

	// need to perform a check whether the required directories exist: If not, need to create them ---
	generate_default_dirtree();

	if(verbose == 1){
		std::cout << "     current_process_ind=" << current_process_ind << std::endl;
		std::cout << "     current_id_ind =" << current_id_ind << std::endl;
		std::cout << "     current_slice_ind =" << current_slice_ind << std::endl;
	}
	
	if(force_manual_config == 1){
		cfg->cfg_file=manual_config_file; // Override the default configuration and the presets
		cfg->read_cfg_file(1); // read the user-specified configuration file
	} else{
		// We assume that the default configuration is already loaded by the class Config
        passed=0;
		// ----- Specify the rules according to the config_presets.cfg ------
		if (processing[current_process_ind] == "Burn-in"){
			for(int i=0; i< cfg->MALA.Nt_learn.size()-1; i++){
				cfg->MALA.Nt_learn[i]=Nt_learn0[i]; // Ensure that the configuration matches the default configuration of these variables
			}
			cfg->MALA.Nt_learn[cfg->MALA.Nt_learn.size()-1]=Nsamples[current_process_ind]+1; // Means never stop to learn..
			cfg->MALA.dN_mixing=dN_mixing0; // Re-establish the value given by default
            
            passed=1;
		}
		if (processing[current_process_ind] == "Learning"){
			for(int i=0; i< cfg->MALA.Nt_learn.size()-1; i++){
				cfg->MALA.Nt_learn[i]=Nt_learn0[i]; // Ensure that the configuration matches the default configuration of these variables
			}
			cfg->MALA.Nt_learn[cfg->MALA.Nt_learn.size()-1]=Nsamples[current_process_ind]+1; // Means never stop to learn..
			cfg->MALA.dN_mixing=Nsamples[current_process_ind]+1; // MEANS THAT WE NEVER MIX THE PARALLEL CHAINS
            
            passed=1;
		}
		if (processing[current_process_ind] == "Acquire"){
			for(int i=0; i<cfg->MALA.Nt_learn.size();i++){
				cfg->MALA.Nt_learn[i]=Nsamples[current_process_ind]+1+i; // Means never learn
				cfg->MALA.dN_mixing=dN_mixing0; // Re-establish the value given by default
			}
            passed=1;
		}
        if(passed == 0){
            std::cout << "Warning: Unrecognized phase name!" << std::endl;
            std::cout << "         Recognized phase names: " << std::endl;
            std::cout << "                - 'Burn-in'" << std::endl;
            std::cout << "                - 'Learning'" << std::endl;
            std::cout << "                - 'Acquire'" << std::endl;
            std::cout << "The program will exit now" << std::endl;
            exit(EXIT_FAILURE);
        }
		if (cfg->outputs.Nbuffer > cfg->outputs.Nsamples){
			cfg->outputs.Nsamples=cfg->outputs.Nbuffer; // If Nsamples<Nbuffer, no need to initialize a large buffer
		}
		if(verbose == 1){
			std::cout << "  - Used configuration for Nt_learn and periods_Nt_learn:" << std::endl;
			std::cout << "      ";
		}
		for(int i=0; i<cfg->MALA.Nt_learn.size();i++){
			std::cout << cfg->MALA.Nt_learn[i] << "  ";
		}
		std::cout << std::endl;

		// We get the object id which is in a Matrix of strings (vectors of vectors)
		process_id=table_ids[current_id_ind].at(0);

		// We define the root_name for all files
		if (Nslices == 1){ // If there is only one slice of data to analyse, no need to add the slice index
			cfg->outputs.output_root_name=process_id +"_" + core_out[current_process_ind] + "_";
			cfg->diags.output_root_name=process_id +"_" + core_out[current_process_ind] + "_";
		} else{ // If there is more than one slice, we need to add the slice index on output files
			cfg->outputs.output_root_name=process_id + "_" + int_to_str(current_slice_ind+1) + "_" + core_out[current_process_ind] + "_";
			cfg->diags.output_root_name=process_id + "_" + int_to_str(current_slice_ind+1) + "_" + core_out[current_process_ind] +  "_";		
		}
		// We define the name of the model file
		cfg->modeling.cfg_model_file=cfg_models_dir +  table_ids[current_id_ind].at(0) + ".model";

		// We define the name of the data file
		cfg->data.data_file=cfg_models_dir + table_ids[current_id_ind].at(0) + ".data";
	
		// Added on 1st July 2023: We define variables required for reading tabulated .priors files, if these are provided
		cfg->modeling.cfg_model_dir=cfg_models_dir;
		cfg->modeling.input_root_name= table_ids[current_id_ind].at(0);
		// -

		// We set the Number of samples and the damping coeficient for the Learning phases
		cfg->outputs.Nsamples=Nsamples[current_process_ind];
		cfg->MALA.c0=c0[current_process_ind];

		// We define the name for the restoration files
		if (Nslices ==1){
			cfg->outputs.restore_file_in=process_id +"_restore_" + core_in[current_process_ind] + "_";
			cfg->outputs.restore_file_out=process_id +"_restore_" + core_out[current_process_ind] + "_";
		} else{
			cfg->outputs.restore_file_in=process_id + "_" + int_to_str(current_slice_ind+1) + "_restore_" + core_in[current_process_ind] +  "_";
			cfg->outputs.restore_file_out=process_id + "_" + int_to_str(current_slice_ind+1)+ "_restore_" + core_out[current_process_ind] + "_";		
		}
		// We deal with restoration conditions
		if(restore[current_process_ind] <= 3){
		
			if(restore[current_process_ind] == 0){
				if(verbose == 1){
					std::cout << "The user requested to start a new analysis ===> PREVIOUS OUTPUT FILES OF SAME NAME MIGHT BE OVERWRITTEN!" << std::endl;
					std::cout << "No restore will be performed" << std::endl;
				}	
				cfg->outputs.do_restore_proposal=0;
				cfg->outputs.do_restore_variables=0;
				cfg->outputs.do_restore_last_index=0;
				cfg->outputs.erase_old_files=1;
			}
			if(restore[current_process_ind] == 1){
				if(verbose == 1){
					std::cout << "The user requested to start a new analysis ===> PREVIOUS OUTPUT FILES OF SAME NAME MIGHT BE OVERWRITTEN!" << std::endl;
					std::cout << "Only the last position in the parameter space will be restored using the user-specified file:" << std::endl;
					std::cout << "                     [current_directory]/Data/restore/" << cfg->outputs.output_root_name + cfg->outputs.restore_file_in << std::endl;
				}
				cfg->outputs.do_restore_proposal=0;
				cfg->outputs.do_restore_variables=1;
				cfg->outputs.do_restore_last_index=0;
				cfg->outputs.erase_old_files=1;
			}
			if(restore[current_process_ind] == 2){
				if(verbose == 1){
					std::cout << "The user requested to start a new analysis ===> PREVIOUS OUTPUT FILES OF SAME NAME MIGHT BE OVERWRITTEN!" << std::endl;
					std::cout << "The last position in the parameter space and the covariance matrix will be restored using the user-specified file:" << std::endl;
					std::cout << "                     [current_directory]/Data/restore/" << cfg->outputs.output_root_name + cfg->outputs.restore_file_in << std::endl;
				}
				cfg->outputs.do_restore_proposal=1;
				cfg->outputs.do_restore_variables=1;
				cfg->outputs.do_restore_last_index=0;
				cfg->outputs.erase_old_files=1;
			}
			if(restore[current_process_ind] == 3){
				if(verbose == 1){
					std::cout << "The user requested to complete an analysis that was previously made ===> PREVIOUS OUTPUT FILES OF SAME NAME WILL BE APPENDED!" << std::endl;
					std::cout << "The last position in the parameter space and the covariance matrix will be restored using the user-specified file:" << std::endl;
					std::cout << "                     [current_directory]/Data/restore/" << cfg->outputs.output_root_name + cfg->outputs.restore_file_in << std::endl;
				}
				cfg->outputs.do_restore_proposal=1;
				cfg->outputs.do_restore_variables=1;
				cfg->outputs.do_restore_last_index=1;
				cfg->outputs.erase_old_files=0;
			}
		} else {
			msg_handler("", "restore", "Config_presets::apply_presets()", "value", 1);
		}
	}
	std::cout << "--------------------------------------------------------------------" << std::endl;
}


std::string Config_presets::format_line(const std::string str){
	/*
	* This program is intended to prepare the string to be splitted
	* and properly formated.
	* [1] Detect if the line has a group indicator.
	*   If it is the case, return the group indicator with the
	*   group indicator symbol
	*   If not,
	* [2] Detect if the line is just a comment line. Comment
	* lines are indicated by a # symbol. 
	*   If not,
	*  	- Remove the comments at the end of a line (indicated 
	*  	  by the symbol ";"
	*  	- Remove blank at the begining/end of the lines
	*   If it is a comment line, returns an empty string

	* 	
	*/

	std::string str_out, str0, char0, terminator;
	size_t pos;

	terminator=";";

	str0=strtrim(str); // remove blanks
	char0=strtrim(str0.substr(0, 1)); 
	if(char0 != "!"){ // If the line has no group indicator
		char0=strtrim(str0.substr(0, 1)); 
		if(char0 != "#" && str0 != ""){ // If the line is not a comment line or an empty line
			pos = str0.find_first_of(terminator); // detect the end of the line
			if(pos != std::string::npos){	
				str_out=str0.substr(0, pos); // remove comments
				str_out=strtrim(str_out); // remove blanks
			} else {
				msg_handler("", "text_terminator", "Config_presets::format_line()", "", 1);
			}
		} else {
			str_out="";
		}
	} else { // If there is a group indicator
		str_out=str0;  // return the line without comments and any blank and without the ":" which terminate
		pos=str0.find_first_of(":"); 
		if(pos != std::string::npos){	
				str_out=str0.substr(0, pos); // remove comments
				str_out=strtrim(str_out); // remove blanks
		} else {
				msg_handler("", "text_terminator", "Config_presets::format_line()", "", 1);
			}
	}

return str_out;
}

void Config_presets::read_cfg_file_presets(){
/*
 * This function is meant to read the config_presets.cfg file
 * This one allows to execute a typical MCMC process, as defined
 * by the parameters contained into the config_default, but with 
 * presets that simplify the configuration and enables to process
 * several configuration/stars at a time
*/
    //int cpt=0;
    int keyword_found;
    std::string line0, char0; // token
    std::vector<std::string> word, tmp_elements;
    std::ifstream cfg_file_session;
    VectorXd sizes;

    cfg_file_session.open(cfg_file_presets.c_str());
    if (cfg_file_session.is_open()) {
	std::cout << "Main presets file opened... processing lines" << std::endl;

	keyword_found=0;
	std::getline(cfg_file_session, line0);
	while(!cfg_file_session.eof()){

		line0=strtrim(line0);
		char0=strtrim(line0.substr(0, 1)); 
		if (line0 == "" || char0 == "#") { 
			std::getline(cfg_file_session, line0); // skip the blank lines or comment lines
		} else{
			char0=strtrim(line0.substr(0, 1)); 
			while(line0 != "/END" && !cfg_file_session.eof()){ // Extract parameters
				line0=format_line(line0);
				word=strsplit(line0, "=");
					if (strtrim(word[0]) == "force_manual_config"){ 
						keyword_found=keyword_found +1;
						force_manual_config=str_to_bool(word[1]);
						std::cout << "      force_manual_config= " << force_manual_config << std::endl;
					}
					if (strtrim(word[0]) == "manual_config_file"){ 
						keyword_found=keyword_found +1;
						manual_config_file=strtrim(word[1]);
						std::cout << "      manual_config_file= " << manual_config_file << std::endl;
					}
					if (strtrim(word[0]) == "cfg_models_dir"){ 
						keyword_found=keyword_found +1;
						cfg_models_dir=strtrim(word[1]);
						std::cout << "      cfg_models_dir= " << cfg_models_dir << std::endl;
					}
					if (strtrim(word[0]) == "cfg_out_dir"){ 
						keyword_found=keyword_found +1;
						cfg_out_dir=strtrim(word[1]);
						std::cout << "      cfg_out_dir= " << cfg_out_dir << std::endl;
					}		
					if (strtrim(word[0]) == "processing"){ 
						keyword_found=keyword_found +1;
						std::cout << "      processing= ";
						word=strsplit(word[1], ",");
						for(int k=0; k<word.size(); k++){
							processing.push_back(strtrim(word[k]));
							std::cout << processing[k] << "  ";
						}
						std::cout << std::endl;
					}	
					if (strtrim(word[0]) == "Nsamples"){ 
						keyword_found=keyword_found +1;
						Nsamples=str_to_Xiarr(word[1], ",");
						std::cout << "      Nsamples= " << Nsamples.transpose() << std::endl;
					}	
					if (strtrim(word[0]) == "c0"){ 
						keyword_found=keyword_found +1;
						c0=str_to_Xdarr(word[1], ",");
						std::cout << "      c0= " << c0.transpose() << std::endl;
					}
					if (strtrim(word[0]) == "restore"){ 
						keyword_found=keyword_found +1;
						restore=str_to_Xiarr(word[1], ",");
						std::cout << "      restore= " << restore.transpose() << std::endl;
					}
					if (strtrim(word[0]) == "core_out"){ 
						keyword_found=keyword_found +1;
						std::cout << "      core_out= ";
						word=strsplit(word[1], ",");
						for(int k=0; k<word.size(); k++){
							core_out.push_back(strtrim(word[k]));
							std::cout << core_out[k] << "  ";
						}
						std::cout << std::endl;
					}
					if (strtrim(word[0]) == "core_in"){ 
						keyword_found=keyword_found +1;
						std::cout << "      core_in= ";
						word=strsplit(word[1], ",");
						for(int k=0; k<word.size(); k++){
							core_in.push_back(strtrim(word[k]));
							std::cout << core_in[k] << "  ";
						}
						std::cout << std::endl;
					}
					if (strtrim(word[0]) == "start_index_processing"){ 
						keyword_found=keyword_found +1;
						std::cout << "      start_index_processing= ";
						first_process_ind=str_to_dbl(word[1]);
							std::cout << first_process_ind << std::endl;
					}
					if (strtrim(word[0]) == "last_index_processing"){ 
						keyword_found=keyword_found +1;
						std::cout << "      last_index_processing= ";
						last_process_ind=str_to_dbl(word[1]);
							std::cout << last_process_ind << std::endl;
					}
					if (strtrim(word[0]) == "table_ids"){ 
						keyword_found=keyword_found +1;
						sizes=str_to_Xdarr(word[1], ",");
						std::cout << "      table_ids= " << sizes.transpose() << std::endl;
						table_ids.resize(sizes[0]);//, sizes[1]);
						for(int i=0; i<sizes[0]; i++){
							std::getline(cfg_file_session, line0); // Read the matrix
							line0=strtrim(line0);
							tmp_elements=strsplit(line0, " ");							
							std::cout << "             [" << i << "]= ";				
							for(int j=0; j<sizes[1]; j++){
								table_ids[i].push_back(strtrim(tmp_elements[j]));
								std::cout << table_ids[i].at(j) << "   ";
							}
							std::cout << std::endl;
							
						}
					}	
				// If we reach the end of file indicator see whether we detected a known keyword
				if ((keyword_found != 13 && line0 == "/END;") || (keyword_found != 13 && cfg_file_session.eof())){
					msg_handler(cfg_file_presets, "text_incorrect_keywords", "Config_presets::read_cfg_file_presets()", "13", 1);
				}
				std::getline(cfg_file_session, line0); // get a new line
				char0=strtrim(line0.substr(0, 1));
			} // endwhile
		} // endelse

	} // endwhile
    } else {
	msg_handler(cfg_file_presets, "openfile", "Config_presets::read_cfg_file_presets()", "Could not open the master configuration file config_presets.cfg !", 1);
    }

	std::cout << "---------------------------------------------- " << std::endl;
	std::cout << std::endl;
}


bool Config_presets::isdir(const std::string pathname){
// Small function that test whether a directory exists
// Return 0 if not, 1 if exists

	struct stat sb;

	if (stat(pathname.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
	{
	    return 1;
	} else {
	    return 0;
	}
}

bool Config_presets::generate_dir_tree(const std::string rootdir, const std::vector<std::string> subdirs){
// Small function that look into the tree of subdirs whether any directory does not exist
// If it is not, it creates it

	bool v, vout;
	std::string path, cmd, rcmd;

	vout=1;
	for(int i=0; i<subdirs.size(); i++){
	        path=rootdir + subdirs[i];
		v=isdir(path);
		if(v == 0){
			cmd="mkdir " + path;
			rcmd=shell_exec(cmd);
			std::cout << rcmd << std::endl;
		} else{
			std::cout << path << " already exists " << std::endl;
		}
		vout=v*vout; // If all dirs exist, vout should be 1. Otherwise 0.
	}

return vout;
}

void Config_presets::generate_default_dirtree(){
// The main function that generate the directories required for the outputs, as defined by config_presets.cpp

	// Define all the directory that have to be tested (and in the correct order)
	std::string root_object_dir;
	std::vector<std::string> subdirs;
	
	root_object_dir="/" + table_ids[current_id_ind].at(0);

	subdirs.push_back(root_object_dir); // The main directory for the studied object
	subdirs.push_back(root_object_dir + "/diags"); 	subdirs.push_back(root_object_dir + "/diags/pdfs");
	subdirs.push_back(root_object_dir + "/restore");
	subdirs.push_back(root_object_dir + "/outputs"); 
	generate_dir_tree(cfg_out_dir, subdirs);

}

std::string Config_presets::shell_exec(const std::string cmd){
// Small function that execute a script on the shell
// Returns the result of the execution on the form of a string
    char buffer[128];
    std::string result = "";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (!feof(pipe)) {
            if (fgets(buffer, 128, pipe) != NULL)
                result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}


/*
 * Function that show an error message depending on error_type
 * It receives the filename, fonction name specific message values (arguments) and returns an error code
 * A negative error code correspond to a fatal error. A positive
 * code correspond to a warning.
 * Exit from the program is controlled by the boolean fatal
*/
int Config_presets::msg_handler(const std::string file, const std::string error_type, const std::string fct_name, const std::string arguments, const bool fatal){

	bool err_msg=0;
	if (fatal == 1){
		std::cout << "Warning Fatal error :" << std::endl;
	}
	if (file != ""){ std::cout << "          Accessing File: " << file << std::endl;}
	if (fct_name !=""){ std::cout << "          Fonction name: " << fct_name << std::endl;}

	if(error_type == "text_incorrect_keywords"){
		std::cout << "           Incorrect number of keywords. Expected number of keywords: " << arguments << std::endl;
		std::cout << "           Probable reason: mispelt keywords ==> unrecognized. Check the syntax of the configuration file " << std::endl;
		err_msg=1;
	}
	if(error_type == "text_terminator"){
		std::cout << "Terminator not found!" << std::endl;
		std::cout << "   - Each uncommented line that not indicate a group must finish by a ; symbol" << std::endl;
		std:: cout << "  - Each group name must finish by a : symbol" << std::endl;
		std::cout << "Comments may follow the terminator symbol in both cases" << std::endl;
		err_msg=1;
	}
	if(error_type =="restore"){
		std::cout << " restore[" << arguments << "] is not allowed" << std::endl;
		std::cout << " restore[i] must be a number not greater than 3" << std::endl;
		err_msg=1;
	}
        if(error_type == "openfile"){
        	std::cout << arguments << std::endl;
        	std::cout << "Check that the destination path exists" << std::endl;
		err_msg=1;
    	}

    	if(error_type == "" || err_msg == 0){
    	    std::cout << "Unknown Fatal Error in " << fct_name << std::endl;
    	    std::cout << "Could not open the file: " << file << std::endl;
    	    std::cout << "Neeed debuging..." << std::endl;
    	}
	if(fatal == 1){
		std::cout << "The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
		return -1;
	} else{
		return 1;
	}
	
}

