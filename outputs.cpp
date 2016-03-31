/*
 * outputs.cpp
 *
 * Contains the class and all kind of functions
 * used to write the data on files
 * 
 *  Created on: 25 Mar 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "outputs.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

//////////////// Constructors //////////////
Outputs::Outputs(std::string diag_cfg_file, int Nchs, int Ndat, int Nsples, VectorXd params, 
		 std::vector<std::string> params_names, VectorXi relax, VectorXi plgth, VectorXd Tchains){

	Nbuffer=10000.; // How many samples are we keeping before saving
	Nchains=Nchs;
	Ndata_obs=Ndat;
	Nsamples=Nsples;
	
	erase_old_files=1;

	get_statcriteria=1;
	get_proposal_params=1;
	get_params=1; // if 1 save all chains
	get_parallel_tempering_params=1;
	get_models=1;
	file_format="text"; // text, bin or HDF5. For the moment only text possible

/*	params_txtbin_fileout="/home/obenomar//Dropbox/Temporary/Cpp-Playground/data_tests/outputs/params";
	proposal_txtbin_fileout="/home/obenomar//Dropbox/Temporary/Cpp-Playground/data_tests/outputs/proposals";
	parallel_tempering_txtbin_fileout="/home/obenomar//Dropbox/Temporary/Cpp-Playground/data_tests/outputs/parallel_tempering";
	model_txtbin_fileout="/home/obenomar//Dropbox/Temporary/Cpp-Playground/data_tests/outputs/models";
	stat_criteria_txtbin_fileout="/home/obenomar//Dropbox/Temporary/Cpp-Playground/data_tests/outputs/stat_criteria";
*/
	output_dir="/Users/obenomar//Dropbox/Temporary/Cpp-Playground/data_tests/outputs/";
	params_txtbin_fileout=output_dir + "params";
	proposal_txtbin_fileout=output_dir + "proposals";
	parallel_tempering_txtbin_fileout=output_dir + "parallel_tempering";
	model_txtbin_fileout=output_dir + "models";
	stat_criteria_txtbin_fileout=output_dir + "stat_criteria";

	hdf5_fileout="";

	restore_file="restore.dat";

	////////////////////////////
	/// Define cons and vars ///
	///////////////////////////
	int Nparams=params.size();
	VectorXd cons0;
	std::vector<std::string> vars_names, cons_names;

	//vars.resize(Nchains, Nparams); // set to maximum possible size. This must be 2D
	cons0.resize(Nparams);
	vars_names.resize(Nparams); // set to maximum possible size this is 1D
	cons_names.resize(Nparams);
	
	Ncons=0;
	Nvars=0;
	for (int i=0; i<Nparams;i++){
		if(relax[i] == 1){
			vars_names[Nvars]=params_names[i];
			Nvars=Nvars+1;
		} else{
			cons0[Ncons]=params(i);
			cons_names[Ncons]=params_names[i];
			Ncons=Ncons+1;
		}
	}

	//std::cout << "Ncons=" << Ncons << std::endl;
	
	vars_names.resize(Nvars); // This is a vector... we can do a normal resize
	if (Ncons>0){
		cons0.conservativeResize(Ncons); // This is also a vector, not a VectorXd or a MatrixXd
		cons_names.resize(Ncons); // again a vector...
	} else {
		cons0.resize(1);
		cons0[0]=-1;
		cons_names.resize(1);
		cons_names[0]="None";
	}
	//std::cout << "cons0=" << cons0 << std::endl;
	////////////////////////////////////////////


	//std::cout << "Ndata=" << Ndata_obs << std::endl;
	//std::cout << "Nsamples=" << Nsamples << std::endl;
	
	//std::cout << "init_buffer_stat_criteria" << std::endl;
	init_buffer_stat_criteria();
	//std::cout << "init_buffer_models" << std::endl;
	init_buffer_models();
	//std::cout << "init_buffer_params" << std::endl;
	init_buffer_params(cons0, vars_names, cons_names, params_names, relax, plgth);
	//std::cout << "init_buffer_parallel_tempering" << std::endl;
	init_buffer_parallel_tempering(Tchains);
	//std::cout << "init_buffer_proposals" << std::endl;
	init_buffer_proposals(vars_names);
	
	//exit(EXIT_SUCCESS);
}

//--------------
void Outputs::init_buffer_proposals(std::vector<std::string> vrs_nmes){

	//std::cout << "init_buffer_proposals()" << std::endl;
	//std::cout << Nvars << std::endl;

	buf_proposal.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_proposal.Ncopy=0;
	buf_proposal.target_acceptance=-1; // Undefined yet... need it as an input
	buf_proposal.sigmas.resize(Nbuffer, Nchains); // (i, Nchains)
	MatrixXd **mat_mus=initialize_3dMatrix(Nbuffer, Nchains, Nvars);
	MatrixXd ***supermat=initialize_4dMatrix(Nbuffer, Nchains, Nvars, Nvars); // (i, Nchains, Nvars, Nvars)

	buf_proposal.mus=mat_mus;
	buf_proposal.covarmats=supermat;
	buf_proposal.Pmoves.resize(Nbuffer, Nchains);

	bool init_value=0;
	//buf_proposal.moveds.resize( Nchains , std::vector<bool>( Nbuffer , init_value ) ); // A matrix of booleans
	buf_proposal.moveds.resize( Nbuffer , std::vector<bool>( Nchains , init_value ) ); // A matrix of booleans

	buf_proposal.vars_names=vrs_nmes;

}

//--------------
void Outputs::init_buffer_parallel_tempering(VectorXd Tchains){
	
	//std::cout << "init_buffer_parallel_tempering()" << std::endl;
	//std::cout << Nvars << std::endl;

	buf_parallel_temp.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_parallel_temp.Ncopy=0; 
	buf_parallel_temp.chain0s.resize(Nbuffer); // by definition chain1 = chain0 + 1... no need to save it then
	buf_parallel_temp.Pswitchs.resize(Nbuffer);
	bool init_value=0;
	buf_parallel_temp.switcheds.resize(Nbuffer); // A vector of booleans
	buf_parallel_temp.attempt_mixing.resize(Nbuffer); 
	buf_parallel_temp.Tcoefs=Tchains; // Invariant toward counts in the MALA algorithm (because the adaptive scheme does not tune the Temperatures.
}

//--------------
void Outputs::init_buffer_params(VectorXd cons_in, std::vector<std::string> vrs_nmes,
					std::vector<std::string> cns_nmes, std::vector<std::string> prms_nmes, VectorXi relax, VectorXi plength){

	//std::cout << "init_buffer_params(()" << std::endl;
	//std::cout << Nvars << std::endl;

	buf_params.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_params.Ncopy=0;

	MatrixXd **vmat=initialize_3dMatrix(Nbuffer, Nchains, Nvars);
	buf_params.vars=vmat; // (i, Nchains, Nvars)
	buf_params.cons=cons_in;
	buf_params.vars_names=vrs_nmes;
	buf_params.cons_names=cns_nmes;
	buf_params.params_names=prms_nmes;
	buf_params.plength=plength;
	
	//exit(EXIT_SUCCESS);
}

//--------------
void Outputs::init_buffer_models(){

	//std::cout << "init_buffer_models()" << std::endl;
	//std::cout << Nvars << std::endl;

	buf_models.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_models.Ncopy=0;
	//std::cout << "Nbuffer = " << Nbuffer << "   Nchains = " << Nchains << "   Ndata = " << Ndata << std::endl;
	MatrixXd **mat=initialize_3dMatrix(Nbuffer, Nchains, Ndata_obs);
	buf_models.models=mat; // (i, Nchains, Ndata)
		
}

//--------------
void Outputs::init_buffer_stat_criteria(){

	//std::cout << "init_buffer_stat_criteria()" << std::endl;
	//std::cout << Nvars << std::endl;

	buf_stat_crit.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_stat_crit.Ncopy=0;
	buf_stat_crit.logLikelihoods.resize(Nbuffer, Nchains);
	buf_stat_crit.logPriors.resize(Nbuffer, Nchains);
	buf_stat_crit.logPosteriors.resize(Nbuffer, Nchains);	

}

                             /////////////// Methods ////////////////

void Outputs::write_txt_prop_params(long Nrest){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_sigmas, filename_moves; //, str_tmp;
	 std::vector< std::string > filename_covarmats, filename_mus;
	 std::ofstream outfile_sigmas, outfile_moves;
	 std::vector<std::ofstream*> outfile_covarmats(Nchains), outfile_mus(Nchains);

	 filename_sigmas=proposal_txtbin_fileout + "_sigmas.txt" ;
	 filename_moves=proposal_txtbin_fileout + "_moves.txt" ;
	
	
	 for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_covarmats.push_back(proposal_txtbin_fileout + "_covarmats_chain-" + ind_str.str() + ".txt") ;
		filename_mus.push_back(proposal_txtbin_fileout + "_mus_chain-" + ind_str.str() + ".txt") ;
		ind_str.str(std::string());
	 }

	if(file_exists(filename_moves.c_str()) && (erase_old_files == 0)){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the sigmas ////////
	//std::cout << "Write the sigmas" << std::endl;
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_sigmas.open(filename_sigmas.c_str()); // Overwrite any existing file
	 } else {
		outfile_sigmas.open(filename_sigmas.c_str(), std::ios::app); // std::app is for append
	 }
	 if (outfile_sigmas.is_open()){
		if (buf_proposal.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
    			outfile_sigmas << "# This is an output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
			outfile_sigmas << "# This file contains only values for sigma[0:Nchains-1]\n" ;
			outfile_sigmas << "! Nchains= " << Nchains << "\n";
		}
		for (int i=0; i<Ntot; i++){ // we write the buffer
			strg.str(std::string());
			strg << buf_proposal.sigmas.row(i)  << "\n";
    			outfile_sigmas << strg.str().c_str();
		}
		outfile_sigmas.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_sigmas.close();
  	}
  	else {
		std::cout << " Unable to open file " << proposal_txtbin_fileout << ".txt" << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	/////// Write the mus ////////
	//std::cout << "Write the mus" << std::endl;
	for (chain=0; chain<Nchains; chain++){

		 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
			outfile_mus[chain]=new std::ofstream((filename_mus[chain]).c_str());
		} else {
			outfile_mus[chain]=new std::ofstream((filename_mus[chain]).c_str(), std::ios::app);
		 }
		 if (outfile_mus[chain]->is_open()){
			if (buf_proposal.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
	    			*outfile_mus[chain] << "# This is an output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
				*outfile_mus[chain] << "# This file contains only values for mu[0:Nchains-1][ 0:Nvars-1]. Each matrix is in a different file, indexed by the chain number\n" ;
				*outfile_mus[chain] << "! Nchains= " << Nchains << "\n";
				*outfile_mus[chain] << "! Nvars= " << Nvars << "\n";
				*outfile_mus[chain] << "! chain= " << chain << "\n";
				// ---- The variable names
				strg.str(std::string());
				strg << "! variable_names=";
				for (int i=0; i<Nvars; i++){
					strg << buf_params.vars_names[i] << "   ";
				}
				strg << "\n";
				*outfile_mus[chain] << strg.str().c_str();
			}
			for (int i=0; i<Ntot; i++){ // we write the buffer	
				strg.str(std::string());	
				strg << (*buf_proposal.mus[i]).row(chain)  << "\n";
	    			*outfile_mus[chain] << strg.str().c_str();
			}
			outfile_mus[chain]->flush(); // Explicitly specify to flush the data into the disk
			strg.str(std::string());
			
			outfile_mus[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << proposal_txtbin_fileout << ".txt" << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

	/////// Write the Pmoves and moveds ////////
	//std::cout << "Write the Pmoves and moveds" << std::endl;
	strg.str(std::string());
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_moves.open(filename_moves.c_str()); // Overwrite any existing file
	 } else {
		outfile_moves.open(filename_moves.c_str(), std::ios::app); // std::app is for append
	 }
	 if (outfile_moves.is_open()){
		if (buf_proposal.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
    			outfile_moves << "# This is an output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
			outfile_moves << "# This file contains only values for Pmove[0:Nchains-1] (first) and for moved[0:Nchains-1] (second group of Nchain values)\n" ;
			outfile_moves << "! Nchains= " << Nchains << "\n";
		}
		for (int i=0; i<Ntot; i++){ // we write the buffer
			strg.str(std::string());
			strg << buf_proposal.Pmoves.row(i) << "     "; // << "   |";
			//std::cout << "\n buf_proposal.moveds[" << i << "] = " <<	std::endl;
			for (int j=0; j<Nchains; j++){ 
				strg << "   " << buf_proposal.moveds[i][j];
			}
			strg << "\n";
    			outfile_moves << strg.str().c_str();
		}
		outfile_moves.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string());
		
		outfile_moves.close();
  	}
  	else {
		std::cout << " Unable to open file " << proposal_txtbin_fileout << ".txt" << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	/////// Write the covarmats ////////
	//std::cout << "Write the covarmats" << std::endl;
	for (chain=0; chain<Nchains; chain++){
		//std::cout << "  chain = " << chain << std::endl; 

		 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
			outfile_covarmats[chain]= new std::ofstream(filename_covarmats[chain].c_str()); // Overwrite any existing file
		 } else {
			outfile_covarmats[chain]= new std::ofstream(filename_covarmats[chain].c_str(), std::ios::app); // std::app is for append
		 }
		 if (outfile_covarmats[chain]->is_open()){
			if (buf_proposal.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
	    			*outfile_covarmats[chain] << "# This is an output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
				*outfile_covarmats[chain] << "# This file contains only values for covarmat[0:Nchains-1][ 0:Nvars-1][ 0:Nvars-1]. Each matrix is in a different file, indexed by the chain number" ;
				*outfile_covarmats[chain] << "# In each file, a new iteration is indicated by !n, with n the iteration number\n";
				*outfile_covarmats[chain] << "! Nchains= " << Nchains << "\n";
				*outfile_covarmats[chain] << "! Nvars= " << Nvars << "\n";
				*outfile_covarmats[chain] << "! chain= " << chain << "\n";
				// ---- The variable names
				strg.str(std::string());
				strg << "! variable_names=";
				for (int i=0; i<Nvars; i++){
					strg << buf_params.vars_names[i] << "   ";
				}
				strg << "\n";
				*outfile_covarmats[chain] << strg.str().c_str();
			}
			for (int i=0; i<Ntot; i++){ // we write the buffer
				strg.str(std::string());
				index= Nbuffer * (buf_proposal.Ncopy) + i;
				*outfile_covarmats[chain] << "*" << index << "\n"; // put the indicator for a new iteration
				for (ind_row=0; ind_row<Nvars; ind_row++){
					//std::cout << "      (*buf_proposal.covarmats[" << i << "][" << chain << "]).row(" << ind_row << ") = " << (*buf_proposal.covarmats[i][chain]).row(ind_row) << std::endl;
					strg << (*buf_proposal.covarmats[i][chain]).row(ind_row) << "\n"; // copy each row and go to the next line
				}
				*outfile_covarmats[chain] << strg.str().c_str();
			}
			outfile_covarmats[chain]->flush(); // Explicitly specify to flush the data into the disk
			strg.str(std::string());
			
			outfile_covarmats[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << proposal_txtbin_fileout << ".txt" << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

}

//------------------
void Outputs::write_txt_params(long Nrest){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::vector< std::string > filename_vars;
	 std::vector<std::ofstream*> outfile_vars(Nchains);

	for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_vars.push_back(params_txtbin_fileout + "_chain-" + ind_str.str() + ".txt") ;
		ind_str.str(std::string());
	 }

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the parameters ////////
	for (chain=0; chain<Nchains; chain++){

		 if (erase_old_files == 1 && buf_params.Ncopy == 0) {
			outfile_vars[chain]=new std::ofstream((filename_vars[chain]).c_str());
		 } else {
			outfile_vars[chain]=new std::ofstream((filename_vars[chain]).c_str(), std::ios::app);
		 }
		 if (outfile_vars[chain]->is_open()){
			if (buf_params.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
	    			*outfile_vars[chain] << "# This is an output file for the model parameters \n";
				*outfile_vars[chain] << "# This file contains values for vars[0:Nchains-1][ 0:Nvars-1]. Each matrix is in a different file, indexed by the chain number\n" ;
				// ----				
				*outfile_vars[chain] << "! Nchains= " << Nchains << "\n";
				*outfile_vars[chain] << "! Nvars= " << Nvars << "\n";
				*outfile_vars[chain] << "! Ncons= " << Ncons << "\n";
				*outfile_vars[chain] << "! chain= " << chain << "\n";
				*outfile_vars[chain] << "! relax= " << buf_params.relax.transpose() << "\n";
				// ---- The constant names
				strg.str(std::string());
				strg << "! plength= ";
				strg << buf_params.plength.transpose() << "\n";
				*outfile_vars[chain] << strg.str().c_str();
				// ---- The constant names
				strg.str(std::string());
				strg << "! constant_names= ";
				for (int i=0; i<buf_params.cons.size(); i++){
					strg << buf_params.cons_names[i] << "   ";
				}
				strg << "\n";
				*outfile_vars[chain] << strg.str().c_str();
				// ---- The constant values
				strg.str(std::string());
				strg << "! constant_values= ";
				if ( buf_params.cons_names[0] == "None"){
					strg << "-1";
				} else {
					strg << buf_params.cons.transpose() << "\n";
				}
				*outfile_vars[chain] << strg.str().c_str();
				// ---- The variable names
				strg.str(std::string());
				strg << "! variable_names=";
				for (int i=0; i<Nvars; i++){
					strg << buf_params.vars_names[i] << "   ";
				}
				strg << "\n";
				*outfile_vars[chain] << strg.str().c_str();
				// ---- Follows below the variable values
			}
			for (int i=0; i<Ntot; i++){ // we write the buffer	
				strg.str(std::string());	
				strg << (*buf_params.vars[i]).row(chain)  << "\n";
	    			*outfile_vars[chain] << strg.str().c_str();
			}
			outfile_vars[chain]->flush(); // Explicitly specify to flush the data into the disk
			strg.str(std::string());
			
			outfile_vars[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << params_txtbin_fileout << "_[chain].txt" << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

}

//------------------
void Outputs::write_txt_parallel_temp_params(long Nrest){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_parallel_temp;
	 std::ofstream outfile_parallel_temp;

 	filename_parallel_temp=parallel_tempering_txtbin_fileout + ".txt" ;

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the parameters of the parallel tempering ////////
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_parallel_temp.open(filename_parallel_temp.c_str()); // Overwrite any existing file
	 } else {
		outfile_parallel_temp.open(filename_parallel_temp.c_str(), std::ios::app); // std::app is for append
	 }
	 if (outfile_parallel_temp.is_open()){
		if (buf_parallel_temp.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
    			outfile_parallel_temp << "# This is an output file for the parameters of the parallel tempering.\n";
			outfile_parallel_temp << "# This file contains values for \n" ;
			outfile_parallel_temp << "# Correspondance between chain0=[0:Nchains-1] and temperature Tcoefs[chain] \n" ;
			strg.str(std::string());
			strg << "! Tcoefs = ";
			strg << buf_parallel_temp.Tcoefs.transpose() << "\n";
			outfile_parallel_temp << strg.str().c_str();
			outfile_parallel_temp << "! labels= attempt_mixing    chain0    Pswitch    switched \n" ;
		}
		for (int i=0; i<Ntot; i++){ // we write the buffer
			strg.str(std::string());
			strg << buf_parallel_temp.attempt_mixing[i]  << "   " << buf_parallel_temp.chain0s[i] << "   " << buf_parallel_temp.Pswitchs[i] 
			     << "   " << buf_parallel_temp.switcheds[i] << "\n";
    			outfile_parallel_temp << strg.str().c_str();
			
		}
		outfile_parallel_temp.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_parallel_temp.close();
  	}
  	else {
		std::cout << " Unable to open file " << parallel_tempering_txtbin_fileout << ".txt" << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
}

//------------------
void Outputs::write_txt_models(long Nrest){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::vector< std::string > filename_models;
	 std::vector<std::ofstream*> outfile_models(Nchains);

	 for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_models.push_back(model_txtbin_fileout + "_models_chain-" + ind_str.str() + ".txt") ;
		ind_str.str(std::string());
	 }


	if(file_exists(filename_models[0].c_str()) && (erase_old_files == 1)){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}
	/////// Write the models ////////
	for (chain=0; chain<Nchains; chain++){
		 //if(buf_models.Ncopy == 0){ // We open the file only at the first execution of the function
		 	if (erase_old_files == 1 && buf_models.Ncopy == 0) {
				outfile_models[chain]= new std::ofstream(filename_models[chain].c_str()); // Overwrite any existing file
			 } else {
				//if(file_exists(filename_mus[chain].c_str())){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists
				outfile_models[chain]= new std::ofstream(filename_models[chain].c_str(), std::ios::app); // std::app is for append
		 	}
		 //} // End of the if with Ncopy
		 if (outfile_models[chain]->is_open()){
			if (buf_models.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
	    			*outfile_models[chain] << "# This is an output file for the models. \n";
				*outfile_models[chain] << "! Nchains= " << Nchains << "\n";
				*outfile_models[chain] << "! Nvars= " << Nvars << "\n";
				*outfile_models[chain] << "! chain= " << chain << "\n";
			}
			for (int i=0; i<Ntot; i++){ // we write the buffer	
				strg.str(std::string());	
				strg << (*buf_models.models[i]).row(chain)  << "\n";
	    			*outfile_models[chain] << strg.str().c_str();
			}
			outfile_models[chain]->flush(); // Explicitly specify to flush the data into the disk
			strg.str(std::string());
			
			outfile_models[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << model_txtbin_fileout << ".txt" << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
	} // End of the loop on chain

}

//------------------
void Outputs::write_txt_stat_criteria(long Nrest){

	
	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_statcrit;
	 std::ofstream outfile_statcrit;
	 std::vector<std::string> labels(3);

	labels[0]= "logLikelihood";
	labels[1]= "logPrior";
	labels[2]= "logPosteriors";
	filename_statcrit=stat_criteria_txtbin_fileout + ".txt" ;

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}
	/////// Write the statistical criteria (logLikelihood, logPrior, logPosterior) ////////
	 if (erase_old_files == 1 && buf_stat_crit.Ncopy == 0) {
		outfile_statcrit.open(filename_statcrit.c_str()); // Overwrite any existing file
	 } else {
		outfile_statcrit.open(filename_statcrit.c_str(), std::ios::app); // std::app is for append
	 }
	 if (outfile_statcrit.is_open()){
		if (buf_stat_crit.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
    			outfile_statcrit << "# This is an output file for the statistical information.\n";
			outfile_statcrit << "# This file contains values for the logLikelihood (columns 0:Nchains-1), logPrior (columns Nchains:2*Nchains-1) and logPosterior (columns 2*Nchains:3*Nchains-1),  \n" ;

			outfile_statcrit << "! Nchains= " << Nchains << "\n";
			outfile_statcrit << "! labels= ";
			for (int k=0; k<labels.size();k++){
				for (int i=0; i<Nchains; i++){			
					outfile_statcrit << labels[k] << "[" << i << "]   ";
				}
			}
			outfile_statcrit << "\n";
		}
		for (int i=0; i<Ntot; i++){ // we write the buffer
			strg.str(std::string());
			strg << buf_stat_crit.logLikelihoods.row(i)  << "     " << buf_stat_crit.logPriors.row(i) << "     " << buf_stat_crit.logPosteriors.row(i) << "\n";
    			outfile_statcrit << strg.str().c_str();
			//std::cout << buf_stat_crit.logLikelihoods.row(i)  << "     " << buf_stat_crit.logPriors.row(i) << "     " << buf_stat_crit.logPosteriors.row(i) << std::endl;
			
		}
		outfile_statcrit.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_statcrit.close();
  	}
  	else {
		std::cout << " Unable to open file " << stat_criteria_txtbin_fileout << ".txt" << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

}

//------------------
//------------------
//------------------
void Outputs::update_buffer_stat_criteria(VectorXd logLikelihood, VectorXd logPrior, VectorXd logPosterior){

 if(get_statcriteria == 1) {
	long rest=Nsamples - Nbuffer*buf_stat_crit.Ncopy;

	if ((buf_stat_crit.counts < Nbuffer) && (buf_stat_crit.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		buf_stat_crit.logLikelihoods.row(buf_stat_crit.counts)=logLikelihood;
		buf_stat_crit.logPriors.row(buf_stat_crit.counts)=logPrior;
		buf_stat_crit.logPosteriors.row(buf_stat_crit.counts)=logPosterior;
		buf_stat_crit.counts=buf_stat_crit.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text"){
			write_txt_stat_criteria(rest);
		}		
		buf_stat_crit.counts=0;
		buf_stat_crit.Ncopy=buf_stat_crit.Ncopy +1;

		std::cout << "     The output file with Statistical information has been updated" << std::endl;

		// The we write the new entry on the empty buffer
		buf_stat_crit.logLikelihoods.row(buf_stat_crit.counts)=logLikelihood;
		buf_stat_crit.logPriors.row(buf_stat_crit.counts)=logPrior;
		buf_stat_crit.logPosteriors.row(buf_stat_crit.counts)=logPosterior;

		buf_stat_crit.counts=buf_stat_crit.counts +1; 
	}
  }
}

void Outputs::update_buffer_models(MatrixXd models_in){

 if(get_models == 1) {
	long rest=Nsamples - Nbuffer*buf_models.Ncopy;
	if ((buf_models.counts < Nbuffer) && (buf_models.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		set_3dMatrix(buf_models.models, models_in, buf_models.counts); // (Nbuffer, Nchains, Ndata_obs);
		buf_models.counts=buf_models.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text") {
			write_txt_models(rest);
		}
		buf_models.counts=0;
		buf_models.Ncopy=buf_models.Ncopy +1;

		std::cout << "     The output file with model information has been updated" << std::endl;

		set_3dMatrix(buf_models.models, models_in, buf_models.counts); // (Nbuffer, Nchains, Ndata_obs);
		buf_models.counts=buf_models.counts +1; 
	}
  }
}

void Outputs::update_buffer_params(MatrixXd params_in){

 if(get_params == 1) {
	long rest=Nsamples - Nbuffer*buf_params.Ncopy;
	if ((buf_params.counts < Nbuffer) && (buf_params.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		set_3dMatrix(buf_params.vars, params_in, buf_params.counts); // (Nbuffer, Nchains, Nvars);
		buf_params.counts=buf_params.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text") {
			write_txt_params(rest);
		}
		buf_params.counts=0;
		buf_params.Ncopy=buf_params.Ncopy +1;
		std::cout << "     The output file with parameters information has been updated" << std::endl;

		set_3dMatrix(buf_params.vars, params_in, buf_params.counts); // (Nbuffer, Nchains, Nvars);
		buf_params.counts=buf_params.counts +1; 

	}
  }
}

void Outputs::update_buffer_ptempering(bool tempted_mixing, int chain_A, long double Probaswitch, bool bool_switched){

 if(get_parallel_tempering_params == 1) {

	long rest=Nsamples - Nbuffer*buf_parallel_temp.Ncopy;
	if ((buf_parallel_temp.counts < Nbuffer) && (buf_parallel_temp.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		buf_parallel_temp.attempt_mixing[buf_parallel_temp.counts]=tempted_mixing; // Tells us if we actually tried to mix the chains (should be the term i%dN_mixing == 0 in MALA::execute())
		buf_parallel_temp.chain0s[buf_parallel_temp.counts]=chain_A;
		buf_parallel_temp.Pswitchs[buf_parallel_temp.counts]=Probaswitch; // Probability of switching chain_A
		buf_parallel_temp.switcheds[buf_parallel_temp.counts]=bool_switched; // Did we switched the chain_A with the chain_A+1 ?
	
		buf_parallel_temp.counts=buf_parallel_temp.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text") {
			write_txt_parallel_temp_params(rest);
		}
		buf_parallel_temp.counts=0;
		buf_parallel_temp.Ncopy=buf_parallel_temp.Ncopy +1;
		std::cout << "     The output file with information on parallel tempering parameters has been updated" << std::endl;

		buf_parallel_temp.attempt_mixing[buf_parallel_temp.counts]=tempted_mixing; // Tells us if we actually tried to mix the chains (should be the term i%dN_mixing == 0 in MALA::execute())
		buf_parallel_temp.chain0s[buf_parallel_temp.counts]=chain_A;
		buf_parallel_temp.Pswitchs[buf_parallel_temp.counts]=Probaswitch; // Probability of switching chain_A
		buf_parallel_temp.switcheds[buf_parallel_temp.counts]=bool_switched; // Did we switched the chain_A with the chain_A+1 ?
	
		buf_parallel_temp.counts=buf_parallel_temp.counts +1; 

	}
	
  }
}

void Outputs::update_buffer_proposals(VectorXd sigma_chains, MatrixXd mu_chains, VectorXd Pmove_chains, 
			std::vector<bool> moved_chains, MatrixXd **covarmat_chains){

  if(get_proposal_params == 1) {
	long rest=Nsamples - Nbuffer*buf_proposal.Ncopy;
	if ((buf_proposal.counts < Nbuffer) && (buf_proposal.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		buf_proposal.sigmas.row(buf_proposal.counts)=sigma_chains;
		set_3dMatrix(buf_proposal.mus, mu_chains, buf_proposal.counts);
		buf_proposal.Pmoves.row(buf_proposal.counts)=Pmove_chains;
		
		for(int m=0; m<Nchains; m++){
			buf_proposal.moveds.at(buf_proposal.counts)[m]=moved_chains[m]; // 2D vector of booleans. can be called by doing something like this: resize[counts][m]
			set_4dMatrix(buf_proposal.covarmats, *covarmat_chains[m], buf_proposal.counts, m); // can be called by doing something like this: covarmats[counts][chain]
		}
		
		buf_proposal.counts=buf_proposal.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		write_txt_prop_params(rest);
		buf_proposal.counts=0;
		buf_proposal.Ncopy=buf_proposal.Ncopy +1;
		std::cout << "     The output file with information on proposal parameters has been updated" << std::endl;

		buf_proposal.sigmas.row(buf_proposal.counts)=sigma_chains;
		set_3dMatrix(buf_proposal.mus, mu_chains, buf_proposal.counts);
		buf_proposal.Pmoves.row(buf_proposal.counts)=Pmove_chains;		

		for(int m=0; m<Nchains; m++){
			buf_proposal.moveds.at(buf_proposal.counts)[m]=moved_chains[m]; // 2D vector of booleans. can be called by doing something like this: resize[counts][m]
			set_4dMatrix(buf_proposal.covarmats, *covarmat_chains[m], buf_proposal.counts, m); // can be called by doing something like this: covarmats[counts][chain]
		}
		
		buf_proposal.counts=buf_proposal.counts +1; 
	}
  }
}


bool Outputs::file_exists(const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

/*
///////////////////// MAIN PROCEDURE FOR TESTS ONLY ////////////////////
int main(){

	std::string diag_cfg_file=""; // not used yet because the configuration is for the moment hard-coded
	int Nchains=4;
	int Nvars=3;
	int Nsamples=203;
	VectorXd cons_in(4);
	//cons_in << 1, 2, 3, 4;
	cons_in[0]=1; cons_in[1]=2; cons_in[2]=3; cons_in[3]=4;

	std::vector<std::string> vars_names(Nvars), cons_names(4), params_names(Nvars + 4);
	VectorXi relax(Nvars + 4);
	VectorXd Tchains;
	Tchains.setLinSpaced(1., 200., Nchains);
	std::vector< std::vector<bool> > test_2dvec( Nchains , std::vector<bool>( Nsamples , 0 ) ); // initialize a 2D vector of booleans

	std::vector<bool> onerow(Nchains, 1); // initialize to 1 a 1D vector of booleans
	test_2dvec[1]=onerow; // update row 1
	std::cout << test_2dvec[1][0] << std::endl; // read its first element


	for(int i=0; i<vars_names.size(); i++){
		vars_names[i]="";
	}
	for(int i=0; i<cons_names.size(); i++){
		cons_names[i]="";
	}
	for(int i=0; i<params_names.size(); i++){
		params_names[i]="";
	}
	relax.setConstant(1);
	for(int i=0; i<4; i++){
		relax[i]=0;
	}

	std::cout << "Testing Outputs::Outputs (The constructor)..." ;
	Outputs diags(diag_cfg_file, Nchains, Nvars, Nsamples, cons_in, 
			vars_names, cons_names, params_names, relax, Tchains);
	std::cout << "Done" << std::endl;

	VectorXd logLikelihood(Nchains), logPrior(Nchains), logPosterior(Nchains), 
		 sigma_chains(Nchains), Pmove_chains(Nchains);

	MatrixXd models_in(Nchains, Nsamples), params_in(Nchains, Nvars), tmp0(Nvars,Nvars), mu_chains(Nchains, Nvars);
	MatrixXd **covarmat_chains=initialize_3dMatrix(Nchains, Nvars, Nvars);

	bool tempted_mixing=1, bool_switched=0;
	std::vector<bool> moved_chains(Nchains);
	int chain_A=3;
	long double Probaswitch=0.4252;

	logLikelihood.setLinSpaced(10, 100);
	logPrior.setLinSpaced(1, 10);
	logPosterior= logLikelihood - logPrior;
	sigma_chains.setConstant(0.1);
	mu_chains.setConstant(1.5);
	Pmove_chains.setConstant(0.65);
	models_in.setConstant(1.);
	params_in.setConstant(10.);
	tmp0.setRandom();
	for(int i=0; i<Nchains; i++){
		moved_chains[i]=1;
		set_3dMatrix(covarmat_chains, tmp0, i);
	}
	std::cout << "Testing In a loop the following update_buffer functions:" << std::endl;
	std::cout << "       - Outputs::update_buffer_stat_criteria..." << std::endl;
	std::cout << "       - Outputs::update_buffer_models..." << std::endl;
	std::cout << "       - Outputs::update_buffer_params..." << std::endl;
	std::cout << "       - Outputs::update_buffer_ptempering..." << std::endl;	
	std::cout << "       - Outputs::update_buffer_proposals..." << std::endl;
	for(int iter=0; iter<Nsamples; iter++){
		//std::cout << "Testing Outputs::update_buffer_stat_criteria...";
		diags.update_buffer_stat_criteria(logLikelihood, logPrior, logPosterior); 	
		
		diags.update_buffer_models(models_in);
		
		diags.update_buffer_params(params_in);
		
		diags.update_buffer_ptempering(tempted_mixing, chain_A, Probaswitch, bool_switched);
		
		diags.update_buffer_proposals(sigma_chains, mu_chains, Pmove_chains, 
			    moved_chains, covarmat_chains);
		std::cout << "[" << iter << "] counts / Ncopy: " << std::endl;
		std::cout << "  - Proposal: " << diags.buf_proposal.counts << " / " << diags.buf_proposal.Ncopy <<  std::endl;
		std::cout << "  - P. temp.: " << diags.buf_parallel_temp.counts << " / " << diags.buf_parallel_temp.Ncopy << std::endl;
		std::cout << "  - Params  : " << diags.buf_params.counts << " / " << diags.buf_params.Ncopy << std::endl;
		std::cout << "  - Models  : " << diags.buf_models.counts << " / " << diags.buf_models.Ncopy << std::endl;
		std::cout << "  - Stat.   : " << diags.buf_stat_crit.counts << " / " << diags.buf_stat_crit.Ncopy << std::endl;
	}
	//std::cout << "Done" << std::endl;

}
*/

