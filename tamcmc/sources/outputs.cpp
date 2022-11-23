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
Outputs::Outputs(Config *cfg, const VectorXd& Tchains){
	
	Nbuffer=cfg->outputs.Nbuffer;
	Nchains=cfg->MALA.Nchains;
	Ndata_obs=cfg->data.data.Nx;
	Nsamples=cfg->outputs.Nsamples;
	
	erase_old_files=cfg->outputs.erase_old_files;

	get_statcriteria=cfg->outputs.get_statcriteria;
	get_proposal_params=cfg->outputs.get_proposal_params;
	get_params=cfg->outputs.get_params; // if 1 save all chains
	get_parallel_tempering_params=cfg->outputs.get_parallel_tempering;
	get_models=cfg->outputs.get_models;
	file_format=cfg->outputs.file_out_format; // text, bin or HDF5. For the moment only text and binary is possible

	output_dir=cfg->outputs.dir_out;
	params_txtbin_fileout=output_dir + cfg->outputs.output_root_name + cfg->outputs.params_txt_fileout;
	proposal_txtbin_fileout=output_dir + cfg->outputs.output_root_name + cfg->outputs.proposal_txt_fileout;
	parallel_tempering_txtbin_fileout=output_dir + cfg->outputs.output_root_name + cfg->outputs.parallel_tempering_txt_fileout;
	model_txtbin_fileout=output_dir + cfg->outputs.output_root_name + cfg->outputs.model_txt_fileout;
	stat_criteria_txtbin_fileout=output_dir + cfg->outputs.output_root_name + cfg->outputs.stat_txt_fileout;
	acceptance_txtbin_fileout=output_dir + cfg->outputs.output_root_name + cfg->outputs.acceptance_txt_fileout; 
	
	if(cfg->outputs.file_out_format == "text"){
		file_ext="txt";
	}
	if(cfg->outputs.file_out_format == "binary"){
		file_ext="bin";
	}
	if(cfg->outputs.file_out_format == "debug"){
		file_ext="dbg";
	}

	if(cfg->outputs.file_out_format != "binary" && cfg->outputs.file_out_format != "text" && cfg->outputs.file_out_format != "debug"){
		std::cout << "  Problem in the file_out_format option" << std::endl;
		std::cout << "  At the moment, only 'text' or 'binary' or 'debug' is accepted (respect the lower case)" << std::endl;
		std::cout << "  The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	//hdf5_fileout=output_dir + cfg->outputs.hdf5_fileout;

	do_restore_variables=cfg->outputs.do_restore_variables;
	do_restore_proposal=cfg->outputs.do_restore_proposal;
	do_restore_last_index=cfg->outputs.do_restore_last_index;

	restore_dir=cfg->outputs.restore_dir;
	//restore_file_in=restore_dir + cfg->outputs.output_root_name + cfg->outputs.restore_file_in; // This gives the name of the core of the restoring output files
	//restore_file_out=restore_dir + cfg->outputs.output_root_name + cfg->outputs.restore_file_out; // This gives the name of the core of the restoring output files
	restore_file_in=restore_dir + cfg->outputs.restore_file_in;
	restore_file_out=restore_dir + cfg->outputs.restore_file_out;

	//std::cout << "Before cons and vars definition" <<std::endl;
	//exit(EXIT_SUCCESS);
	////////////////////////////
	/// Define cons and vars ///
	///////////////////////////
	int Nparams=cfg->modeling.inputs.inputs.size();
	VectorXd cons0;
	std::vector<std::string> vars_names, cons_names;

	cons0.resize(Nparams);
	vars_names.resize(Nparams); // set to maximum possible size this is 1D
	cons_names.resize(Nparams);
	
	Ncons=0;
	Nvars=0;
	for (int i=0; i<Nparams;i++){
		if(cfg->modeling.inputs.relax[i] == 1){
			vars_names[Nvars]=cfg->modeling.inputs.inputs_names[i];
			Nvars=Nvars+1;
		} else{
			cons0[Ncons]=cfg->modeling.inputs.inputs(i);
			cons_names[Ncons]=cfg->modeling.inputs.inputs_names[i];
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

	//std::cout << "cons0=" << cons0 << std::endl;\
	////////////////////////////////////////////
	Ndata_obs_limit=5000; // Condition that defines when we allow to save models. If Ndata_obs > Ndata_obs_limit, we NEVER save the models
	if (Ndata_obs <= Ndata_obs_limit && get_models == 1){
		std::cout << "       ************ WARNING **********" << std::endl;
		std::cout << "          Too many data points! We    " << std::endl;
		std::cout << "        We override of get_models = 1 "  << std::endl;
		std::cout << "           No models will be saved    " << std::endl;
		std::cout << "       *******************************" << std::endl;
		get_models=0;
	}

	std::cout << "Ndata=" << Ndata_obs << std::endl;
	std::cout << "Nsamples=" << Nsamples << std::endl;
	
	//std::cout << "init_buffer_stat_criteria" << std::endl;
	init_buffer_stat_criteria();
	//std::cout << "init_buffer_models" << std::endl;
	init_buffer_models();
	//std::cout << "init_buffer_params" << std::endl;
	init_buffer_params(cons0, vars_names, cons_names, cfg->modeling.inputs.inputs_names, cfg->modeling.inputs.relax, cfg->modeling.inputs.plength);
	//std::cout << "init_buffer_parallel_tempering" << std::endl;
	init_buffer_parallel_tempering(Tchains);
	//std::cout << "init_buffer_proposals" << std::endl;
	init_buffer_proposals(vars_names);
	init_buffer_acceptance();
	init_buffer_restore(cfg->restored_vals.iteration);

	//std::cout << "End of initialisation of the outputs" <<std::endl;
	//exit(EXIT_SUCCESS);
}

Outputs::~Outputs(){


	//destroy_3dMatrix(buf_restore.covarmats, Nchains);
	destroy_3dMatrix(buf_proposal.mus, Nbuffer);
	destroy_4dMatrix(buf_proposal.covarmats, Nbuffer, Nchains);
	destroy_3dMatrix(buf_params.vars, Nbuffer);
	if(get_models == 1){
		destroy_3dMatrix(buf_models.models, Nbuffer);
	} else{
		destroy_3dMatrix(buf_models.models, 1); // This Tensor is of depth 1 if no model is saved
	}
	//std::cout << "All destroyed" << std::endl;

	buf_restore.covarmats=NULL;
	buf_proposal.mus=NULL;
	buf_proposal.covarmats=NULL;
	buf_params.vars=NULL;
	buf_models.models=NULL;
}


void Outputs::destroy_3dMatrix(MatrixXd** m3d, const int depth){

	for(int d1=0; d1 < depth; d1++){
		delete m3d[d1];
	}
	delete m3d;
}

void Outputs::destroy_4dMatrix(MatrixXd*** m4d, const int depth, const int Nch){
	
	for(int d1=0; d1< depth; d1++){
		for(int d2=0; d2<Nch; d2++){
			delete m4d[d1][d2];
		}
		delete m4d[d1];
	}
	delete m4d;
}

//--------------
void Outputs::init_buffer_acceptance(){

	buf_acceptance.Ncopy=0;
	buf_acceptance.xaxis=0;
	buf_acceptance.acceptance_rate.resize(Nchains);
}

void Outputs::init_buffer_proposals(const std::vector<std::string> vrs_nmes){

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
	buf_proposal.moveds.resize( Nbuffer , std::vector<bool>( Nchains , init_value ) ); // A matrix of booleans

	buf_proposal.vars_names=vrs_nmes;

}

//--------------
void Outputs::init_buffer_parallel_tempering(const VectorXd& Tchains){

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
void Outputs::init_buffer_params(const VectorXd& cons_in, const std::vector<std::string> vrs_nmes,
					const std::vector<std::string> cns_nmes, const std::vector<std::string> prms_nmes, const VectorXi& relax, const VectorXi& plength){

	buf_params.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_params.Ncopy=0;

	MatrixXd **vmat=initialize_3dMatrix(Nbuffer, Nchains, Nvars);
	buf_params.vars=vmat; // (i, Nchains, Nvars)
	buf_params.cons=cons_in;
	buf_params.vars_names=vrs_nmes;
	buf_params.cons_names=cns_nmes;
	buf_params.params_names=prms_nmes;
	buf_params.relax=relax;
	buf_params.plength=plength;
	
}

//--------------
void Outputs::init_buffer_models(){

	buf_models.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_models.Ncopy=0;
	if(get_models == 1){
		MatrixXd **mat=initialize_3dMatrix(Nbuffer, Nchains, Ndata_obs);
		buf_models.models=mat; // (i, Nchains, Ndata)
	} else{
		MatrixXd **mat=initialize_3dMatrix(1, Nchains, Ndata_obs); // Declase the model buffer of size 1 to save memory
		buf_models.models=mat; // (i, Nchains, Ndata)
	}
	
		
}

//--------------
void Outputs::init_buffer_stat_criteria(){

	buf_stat_crit.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_stat_crit.Ncopy=0;
	buf_stat_crit.logLikelihoods.resize(Nbuffer, Nchains);
	buf_stat_crit.logPriors.resize(Nbuffer, Nchains);
	buf_stat_crit.logPosteriors.resize(Nbuffer, Nchains);	

}


void Outputs::init_buffer_restore(const long init_i){

	buf_restore.Nsamples_sofar=init_i;
	buf_restore.counts=0; // How many samples are written so far? Should be iteratively updated
	buf_restore.Ncopy=0;
	
	buf_restore.sigmas.resize(Nchains);
	MatrixXd **cov=initialize_3dMatrix(Nchains, Nvars, Nvars);
	buf_restore.mus.resize(Nchains, Nvars);
	buf_restore.covarmats=cov;

	buf_restore.sigmas_mean.resize(Nchains);
	buf_restore.mus_mean.resize(Nchains, Nvars);
	buf_restore.covarmats_mean=cov;


	buf_restore.sigmas_mean.setZero();
	buf_restore.mus_mean.setZero();
	for (int chain=0; chain<Nchains; chain++){
		buf_restore.covarmats_mean[chain]->setZero();
	}
	
	buf_restore.vars.resize(Nchains, Nvars);
	buf_restore.vars_mean.resize(Nchains, Nvars);
	buf_restore.vars_mean.setZero();

}

                             /////////////// Methods FOR WRITTING IN PLAIN ASCII ////////////////

void Outputs::write_txt_prop_params(const long Nrest, const bool dbg){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot, it;
	 std::ostringstream strg, ind_str;
	 std::string filename_sigmas, filename_moves; //, str_tmp;
	 std::vector< std::string > filename_covarmats, filename_mus;
	 std::ofstream outfile_sigmas, outfile_moves;
	 std::vector<std::ofstream*> outfile_covarmats(Nchains), outfile_mus(Nchains);

	std::string ext_all;
	if(dbg == 1){	
        	 ext_all=file_ext + ".txt";

	} else{
		ext_all=file_ext;
	}
	 filename_sigmas=proposal_txtbin_fileout + "_sigmas." + ext_all; //".txt" ;
	 filename_moves=proposal_txtbin_fileout + "_moves." + ext_all; //"txt" ;


	 for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_covarmats.push_back(proposal_txtbin_fileout + "_covarmats_chain-" + ind_str.str() + "." + ext_all) ;
		filename_mus.push_back(proposal_txtbin_fileout + "_mus_chain-" + ind_str.str() + "." + ext_all) ;

		ind_str.str(std::string());
	 }

	if(file_exists(filename_moves.c_str()) && (erase_old_files == 0)){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists
	
	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the sigmas ////////
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_sigmas.open(filename_sigmas.c_str()); // Overwrite any existing file in PLAIN ASCII
	 } else {
		outfile_sigmas.open(filename_sigmas.c_str(), std::ofstream::app); // std::app is for append
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
		std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << ext_all << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	/////// Write the mus ////////
	for (chain=0; chain<Nchains; chain++){

		 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
			outfile_mus[chain]=new std::ofstream((filename_mus[chain]).c_str()); // Write in plain ascii
		} else {
			outfile_mus[chain]=new std::ofstream((filename_mus[chain]).c_str(), std::ofstream::app);
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
			std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << ext_all << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

	/////// Write the Pmoves and moveds ////////
	 strg.str(std::string());
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_moves.open(filename_moves.c_str()); // Overwrite any existing file in PLAIN ASCII
	 } else {
		outfile_moves.open(filename_moves.c_str(), std::ofstream::app); // std::app is for append
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
		std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << ext_all << std::endl;	
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
			outfile_covarmats[chain]= new std::ofstream(filename_covarmats[chain].c_str(), std::ofstream::app); // std::app is for append
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
			std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << ext_all << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

	for( std::vector<std::ofstream*>::iterator it=outfile_mus.begin() ; it != outfile_mus.end(); ++it){
		delete (*it);
	}
	for( std::vector<std::ofstream*>::iterator it=outfile_covarmats.begin() ; it != outfile_covarmats.end(); ++it){
		delete (*it);
	}
	
}

//------------------
void Outputs::write_txt_params(const long Nrest, const bool dbg){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::vector< std::string > filename_vars;
	 std::vector<std::ofstream*> outfile_vars(Nchains);
	
	std::string ext_all;
	if(dbg == 1){	
        	 ext_all=file_ext + ".txt";

	} else{
		ext_all=file_ext;
	}
	for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_vars.push_back(params_txtbin_fileout + "_chain-" + ind_str.str() + "." + ext_all) ;
		ind_str.str(std::string());
	 }

	if(erase_old_files == 0){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists
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
			//outfile_vars[chain]=new std::ofstream((filename_vars[chain]).c_str(), std::ios::app);
			outfile_vars[chain]=new std::ofstream((filename_vars[chain]).c_str(), std::ofstream::app);
		 }
		 if (outfile_vars[chain]->is_open()){
			if (buf_params.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
	    			*outfile_vars[chain] << "# This is an output file for the model parameters \n";
				*outfile_vars[chain] << "# This file contains values for vars[0:Nchains-1][ 0:Nvars-1]. Each matrix is in a different file, indexed by the chain number\n" ;
				// ----				
				*outfile_vars[chain] << "! Nsamples= " << Nsamples << "\n";
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
					strg << "-1\n";
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
			std::cout << " Unable to open file " << params_txtbin_fileout << "_[chain]." + ext_all << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

	for( std::vector<std::ofstream*>::iterator it=outfile_vars.begin() ; it != outfile_vars.end(); ++it){
		delete (*it);
	}
}

//------------------
void Outputs::write_txt_parallel_temp_params(const long Nrest, const bool dbg){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_parallel_temp;
	 std::ofstream outfile_parallel_temp;
	std::string ext_all;

	if(dbg == 1){	
        	 ext_all=file_ext + ".txt";

	} else{
		ext_all=file_ext;
	}

 	filename_parallel_temp=parallel_tempering_txtbin_fileout + "." + ext_all;

	if(erase_old_files == 0){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists
	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the parameters of the parallel tempering ////////
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_parallel_temp.open(filename_parallel_temp.c_str()); // Overwrite any existing file in PLAIN ASCII
	 } else {
		outfile_parallel_temp.open(filename_parallel_temp.c_str(), std::ofstream::app); // std::app is for append
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
		std::cout << " Unable to open file " << parallel_tempering_txtbin_fileout << "." + ext_all << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
}

//------------------
void Outputs::write_txt_models(const long Nrest, const bool dbg){

	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::vector< std::string > filename_models;
	 std::vector<std::ofstream*> outfile_models(Nchains);

	std::string ext_all;
	if(dbg == 1){	
        	 ext_all=file_ext + ".txt";

	} else{
		ext_all=file_ext;
	}
	 for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_models.push_back(model_txtbin_fileout + "_models_chain-" + ind_str.str() + "." + ext_all) ;
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
				//outfile_models[chain]= new std::ofstream(filename_models[chain].c_str(), std::ios::app); // std::app is for append
				outfile_models[chain]= new std::ofstream(filename_models[chain].c_str(), std::ofstream::app); // std::app is for append
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
			std::cout << " Unable to open file " << model_txtbin_fileout << "." + ext_all << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
	} // End of the loop on chain

	for( std::vector<std::ofstream*>::iterator it=outfile_models.begin() ; it != outfile_models.end(); ++it){
		delete (*it);
	}
}

//------------------
void Outputs::write_txt_acceptance(){

	 bool need_header=1; // By defaut we write a header
	 std::ostringstream strg;
	 std::string filename_acceptance; 
	 std::ofstream outfile_acceptance;

	filename_acceptance=acceptance_txtbin_fileout + ".txt" ;

	if(file_exists(filename_acceptance.c_str()) && (erase_old_files == 0)){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists

	/////// Write the acceptance_rate ////////
	 if (erase_old_files == 1 && buf_acceptance.Ncopy == 0) {
		outfile_acceptance.open(filename_acceptance.c_str()); // Overwrite any existing file in PLAIN ASCII
	 } else {
		outfile_acceptance.open(filename_acceptance.c_str(), std::ofstream::app); // std::app is for append
	 }
	 if (outfile_acceptance.is_open()){
		if (buf_acceptance.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
    			outfile_acceptance << "# This is an output file for the acceptance rate. \n";
			outfile_acceptance << "# This file contains values for the acceptance_rate[0:Nchains-1] in function of the average sample position\n" ;
			outfile_acceptance << "# Averaging is done over Nbuffer \n";
			outfile_acceptance << "! Nchains= " << Nchains << "\n";
		}
		strg.str(std::string());
		strg << buf_acceptance.xaxis;
		strg << " ";
		strg << buf_acceptance.acceptance_rate.transpose() << "\n";

		outfile_acceptance << strg.str().c_str();
		outfile_acceptance.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_acceptance.close();
  	}
  	else {
		std::cout << " Unable to open file " << acceptance_txtbin_fileout << "." + file_ext << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

}

void Outputs::write_txt_stat_criteria(const long Nrest, const bool dbg){

	
	 bool need_header=1; // By defaut we write a header
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_statcrit;
	 std::ofstream outfile_statcrit;
	 std::vector<std::string> labels(3);

	std::string ext_all;

	if(dbg == 1){	
        	 ext_all=file_ext + ".txt";

	} else{
		ext_all=file_ext;
	}

	labels[0]= "logLikelihood";
	labels[1]= "logPrior";
	labels[2]= "logPosteriors";

	filename_statcrit=stat_criteria_txtbin_fileout + "." + ext_all ;

	if(erase_old_files == 0){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}
	/////// Write the statistical criteria (logLikelihood, logPrior, logPosterior) ////////
	 if (erase_old_files == 1 && buf_stat_crit.Ncopy == 0) {
		outfile_statcrit.open(filename_statcrit.c_str()); // Overwrite any existing file
	 } else {
		outfile_statcrit.open(filename_statcrit.c_str(), std::ofstream::app); // std::app is for append
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
			
		}
		outfile_statcrit.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_statcrit.close();
  	}
  	else {
		std::cout << " Unable to open file " << stat_criteria_txtbin_fileout << "." + ext_all << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

}

void Outputs::write_buffer_restore(){

	std::string filename_vars= restore_file_out + "1.dat"; // 1 corresponds to the variables 
	std::string filename_musigmas= restore_file_out + "2.dat"; // 2 corresponds to sigma and mu
	std::string filename_covarmats= restore_file_out + "3.dat"; // 3 corresponds to the covariance matrix

	std::ostringstream strg;
	std::ofstream outfile_vars, outfile_musigmas, outfile_covarmats;

	/////// Write the variables ////////
	outfile_vars.open(filename_vars.c_str()); 
	if (outfile_vars.is_open()){
    		outfile_vars << "# This is an output file containing what is required to restore a run to its last saved position \n";
    		outfile_vars << "# File number: 1 \n";
    		outfile_vars << "# Contains the last values for the variables vars[0:Nchain-1]. vars_mean denotes averaged values of Nbuffer \n";
		outfile_vars << "# Use this if you wish to: \n";
		outfile_vars << "#       (1) complete a finished job that requires more samples ==> set erase_old_file=0 and do_restore_[X]=1 \n" ;
		outfile_vars << "#       (2) restart a finished job by ignoring old samples (e.g. ignoring a Burn-in) ==> set erase_old_file=1 and do_restore_proposal=1 \n" ;
		outfile_vars << "#       (3) terminate an unfinished job which failed to finished (e.g. due to computer unexpected shutdown) ==> set erase_old_file=0 and do_restore_[X]=1 \n" ;
		outfile_vars << "! Nchains= " << Nchains << "\n";
		outfile_vars << "! Nvars= " << Nvars << "\n";
		outfile_vars << "! iteration=" << Nbuffer * (buf_restore.Ncopy) + buf_restore.counts + buf_restore.Nsamples_sofar << "\n";
		// ---- The variable names
		strg.str(std::string());
		strg << "! variable_names=";
		for (int i=0; i<Nvars; i++){
			strg << buf_params.vars_names[i] << "   ";
		}
		strg << "\n";
		outfile_vars << strg.str().c_str();
		outfile_vars << "! vars= " << "\n";
		for (int i=0; i<Nchains; i++){ // we write the parameters... one line is one chain
			strg.str(std::string());
			strg << buf_restore.vars.row(i)  << "\n";
    			outfile_vars << strg.str().c_str();
		}
		// --- Put the averaged values ---
		strg.str(std::string());
		outfile_vars << "! vars_mean= " << "\n";
		for (int i=0; i<Nchains; i++){ // we write the parameters... one line is one chain
			strg.str(std::string());
			strg << buf_restore.vars_mean.row(i)  << "\n";
    			outfile_vars << strg.str().c_str();
		}
		
		outfile_vars.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_vars.close();
  	}
  	else {
		std::cout << " Unable to open file " << filename_vars  << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	/////// Write the sigmas ////////
	outfile_musigmas.open(filename_musigmas.c_str()); 
	if (outfile_musigmas.is_open()){
    		outfile_musigmas << "# This is an output file containing what is required to restore a run to its last saved position \n";
    		outfile_musigmas << "# File number: 2 \n";
    		outfile_musigmas << "# Contains the last values of (a) sigmas[0:Nchains-1] and (b) mus[0:Nchains-1, 0:Nvars-1].  sigmas_mean and mus_mean denotes averaged values of Nbuffer\n";
		outfile_musigmas << "# Use this if you wish to: \n";
		outfile_musigmas << "#       (1) complete a finished job that requires more samples ==> set erase_old_file=0 and do_restore=1 \n" ;
		outfile_musigmas << "#       (2) restart a finished job by ignoring old samples (e.g. ignoring a Burn-in) ==> set erase_old_file=1 and do_restore=1 \n" ;
		outfile_musigmas << "#       (3) terminate an unfinished job which failed to finished (e.g. due to computer unexpected shutdown) ==> set erase_old_file=0 and do_restore=1 \n" ;
		outfile_musigmas << "! Nchains= " << Nchains << "\n";
		outfile_musigmas << "! Nvars= " << Nvars << "\n";
		outfile_musigmas << "! iteration=" << Nbuffer * (buf_restore.Ncopy) + buf_restore.counts + buf_restore.Nsamples_sofar << "\n";
		// ---- The variable names
		strg.str(std::string());
		strg << "! variable_names=";
		for (int i=0; i<Nvars; i++){
			strg << buf_params.vars_names[i] << "   ";
		}
		strg << "\n";
		outfile_musigmas << strg.str().c_str();
		
		// --- sigmas and mus ---
		outfile_musigmas << "! sigmas= " << buf_restore.sigmas.transpose() << "\n";
		outfile_musigmas << "! mus= " << "\n";
		for (int chain=0; chain<Nchains; chain++){ // we write the buffer
			strg.str(std::string());
			strg << buf_restore.mus.row(chain)  << "\n";
    			outfile_musigmas << strg.str().c_str();
		}
		// --- Put the averaged values ---
		outfile_musigmas << "! sigmas_mean= " << buf_restore.sigmas_mean.transpose() << "\n";
		outfile_musigmas << "! mus_mean= " << "\n";
		for (int chain=0; chain<Nchains; chain++){ // we write the buffer
			strg.str(std::string());
			strg << buf_restore.mus_mean.row(chain)  << "\n";
    			outfile_musigmas << strg.str().c_str();
		}
		
		outfile_musigmas.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_musigmas.close();
  	}
  	else {
		std::cout << " Unable to open file " << filename_musigmas  << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	/////// Write the covariance matrix ////////
	//std::cout << "Write the covarmat" << std::endl;
	outfile_covarmats.open(filename_covarmats.c_str());
	if (outfile_covarmats.is_open()){
		outfile_covarmats << "# This is an output file containing what is required to restore a run to its last saved position \n";
    		outfile_covarmats << "# File number: 3 \n";
    		outfile_covarmats << "# Contains the last value of the covariance matrix covarmats[0:Nchains-1, 0:Nvars-1, 0:Nvars-1]. covarmats_mean denotes the averaged values over Nbuffer\n";
		outfile_covarmats << "# Use this if you wish to: \n";
		outfile_covarmats << "#       (1) complete a finished job that requires more samples ==> set erase_old_file=0 and do_restore=1 \n" ;
		outfile_covarmats << "#       (2) restart a finished job by ignoring old samples (e.g. ignoring a Burn-in) ==> set erase_old_file=1 and do_restore=1 \n" ;
		outfile_covarmats << "#       (3) terminate an unfinished job which failed to finished (e.g. due to computer unexpected shutdown) ==> set erase_old_file=0 and do_restore=1 \n" ;
		outfile_covarmats << "! Nchains= " << Nchains << "\n";
		outfile_covarmats << "! Nvars= " << Nvars << "\n";
		outfile_covarmats << "! iteration=" << Nbuffer * (buf_restore.Ncopy) + buf_restore.counts + buf_restore.Nsamples_sofar << "\n";
		// ---- The variable names ----
		strg.str(std::string());
		strg << "! variable_names=";
		for (int i=0; i<Nvars; i++){
			strg << buf_params.vars_names[i] << "   ";
		}
		strg << "\n";
		outfile_covarmats << strg.str().c_str();
		outfile_covarmats << "! covarmats= " << "\n";
		for (int chain=0; chain<Nchains; chain++){
			outfile_covarmats << "*" << chain << "\n"; // put the indicator for the chain
			for (int i=0; i<Nvars; i++){ // we write the covariance matrix	
				strg.str(std::string());	
				strg << buf_restore.covarmats[chain]->row(i)  << "\n";
	    			outfile_covarmats << strg.str().c_str();
			}
		}
		outfile_covarmats.flush(); // Explicitly specify to flush the data into the disk
		// --- Put the averaged values ---
		
		outfile_covarmats << "! covarmats_mean= " << "\n";
		for (int chain=0; chain<Nchains; chain++){
			outfile_covarmats << "*" << chain << "\n"; // put the indicator for the chain
			for (int i=0; i<Nvars; i++){ // we write the covariance matrix	
				strg.str(std::string());	
				strg << buf_restore.covarmats_mean[chain]->row(i)  << "\n";
	    			outfile_covarmats << strg.str().c_str();
			}
		}
		outfile_covarmats.flush(); // Explicitly specify to flush the data into the disk
		
		strg.str(std::string());	
		outfile_covarmats.close();
  	} // End of the if with is_open()
  	else {
		std::cout << " Unable to open file " << filename_covarmats << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
}

                             /////////////// Methods FOR WRITTING IN BINARY ////////////////

void Outputs::write_bin_prop_params(const long Nrest){

	bool boolval=0;
	size_t size_bool=sizeof(boolval);
	double dblval=0.0;
	size_t size_dbl=sizeof(dblval);

	 bool need_header=1; // Variable for handling the header (which is in ASCII)... this is the metadata.
	 int chain, ind_row, ind_col, index, Ntot, it;
	 std::ostringstream strg, ind_str;
	 std::string filename_sigmas, filename_moves, filename_sigmas_hder, filename_moves_hder, filename_covarmats_hder, filename_mus_hder;
	 std::vector< std::string > filename_covarmats, filename_mus;
	 std::ofstream outfile_sigmas, outfile_moves, outfile_sigmas_hder, outfile_moves_hder, outfile_covarmats_hder, outfile_mus_hder;
	 std::vector<std::ofstream*> outfile_covarmats(Nchains), outfile_mus(Nchains);

	 filename_sigmas=proposal_txtbin_fileout + "_sigmas." + file_ext; //".txt" ;
	 filename_moves=proposal_txtbin_fileout + "_moves." + file_ext; //"txt" ;

	 filename_sigmas_hder=proposal_txtbin_fileout + "_sigmas.hdr"; // Header in ASCII
	 filename_moves_hder=proposal_txtbin_fileout + "_moves.hdr";  // Header in ASCII
 	 filename_covarmats_hder=proposal_txtbin_fileout + "_covarmats" + ".hdr"; // Common Header to all chains in ASCII
 	 filename_mus_hder=proposal_txtbin_fileout + "_mus" ".hdr"; // Common Header to all chains in ASCII

	 for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_covarmats.push_back(proposal_txtbin_fileout + "_covarmats_chain-" + ind_str.str() + "." + file_ext) ; // Binary data
		filename_mus.push_back(proposal_txtbin_fileout + "_mus_chain-" + ind_str.str() + "." + file_ext) ; // Binary data
		ind_str.str(std::string());
	 }

	if(file_exists(filename_moves.c_str()) && (erase_old_files == 0)){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the sigmas ////////
	if (need_header == 1){ // Write the header if requested (should always be on actually)
		outfile_sigmas_hder.open(filename_sigmas_hder.c_str());  // Write the header ASCII file
    		outfile_sigmas_hder << "# This is the header of the BINARY output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
		outfile_sigmas_hder << "# This file contains only values for sigma[0:Nchains-1]\n" ;
		outfile_sigmas_hder << "! Nchains= " << Nchains << "\n";
		outfile_sigmas_hder << "! Nsamples_done=" << Nbuffer * (buf_proposal.Ncopy) + buf_proposal.counts + buf_restore.Nsamples_sofar << "\n";
		outfile_sigmas_hder.close();
	}
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
			outfile_sigmas.open(filename_sigmas.c_str(), std::ofstream::binary); // Overwrite any existing file in BINARY
	 } else {
			//outfile_sigmas.open(filename_sigmas.c_str(), std::ios::app | std::ios::binary);
			outfile_sigmas.open(filename_sigmas.c_str(), std::ofstream::app | std::ofstream::binary);
	 }
	 if (outfile_sigmas.is_open()){
		for (int i=0; i<Ntot; i++){ // we write the buffer in the binary file
			for(int j=0; j<Nchains; j++){
				outfile_sigmas.write(reinterpret_cast<char*>(&buf_proposal.sigmas(i,j)), size_dbl);
			}
		}
		outfile_sigmas.flush(); // Explicitly specify to flush the data into the disk
		outfile_sigmas.close();
  	}
  	else {
		std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << file_ext << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	/////// Write the mus ////////
	if (need_header == 1){ 
		outfile_mus_hder.open((filename_mus_hder).c_str(), std::ofstream::app);    		
		outfile_mus_hder << "# This is the header of the BINARY output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
		outfile_mus_hder << "# This file contains only values for mu[0:Nchains-1][ 0:Nvars-1]. Each matrix is in a different file, indexed by the chain number\n" ;
		outfile_mus_hder << "! Nchains= " << Nchains << "\n";
		outfile_mus_hder << "! Nvars= " << Nvars << "\n";
		outfile_mus_hder << "! Nsamples_done=" << Nbuffer * (buf_proposal.Ncopy) + buf_proposal.counts + buf_restore.Nsamples_sofar << "\n";
		//outfile_mus_hder << "! chain= " << chain << "\n";
		// ---- The variable names
		strg.str(std::string());
		strg << "! variable_names=";
		for (int i=0; i<Nvars; i++){
			strg << buf_params.vars_names[i] << "   ";
		}
		strg << "\n";
		outfile_mus_hder << strg.str().c_str();
		outfile_mus_hder.close();
	}
	for (chain=0; chain<Nchains; chain++){
		 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
			outfile_mus[chain]=new std::ofstream((filename_mus[chain]).c_str(), std::ofstream::binary); // Write in Binary
		} else {
			outfile_mus[chain]=new std::ofstream((filename_mus[chain]).c_str(), std::ofstream::app | std::ofstream::binary);
		 }
		 if (outfile_mus[chain]->is_open()){
		   for (int i=0; i<Ntot; i++){ // we write the buffer
			for(int j=0; j<(*buf_proposal.mus[i]).cols(); j++){
				outfile_mus[chain]->write( reinterpret_cast<char*>(& (*buf_proposal.mus[i])(chain,j) ), size_dbl );
			}
		   }
		   outfile_mus[chain]->flush(); // Explicitly specify to flush the data into the disk
		   outfile_mus[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << file_ext << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

	/////// Write the Pmoves and moveds ////////
	if (need_header == 1){
		outfile_moves_hder.open((filename_moves_hder).c_str(), std::ofstream::app);    	
    		outfile_moves_hder << "# This is the header of the BINARY output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
		outfile_moves_hder << "# This file contains only values for Pmove[0:Nchains-1] (first) and for moved[0:Nchains-1] (second group of Nchain values)\n" ;
		outfile_moves_hder << "! Nchains= " << Nchains << "\n";
		outfile_moves_hder << "! Nsamples_done=" << Nbuffer * (buf_proposal.Ncopy) + buf_proposal.counts + buf_restore.Nsamples_sofar << "\n";
		outfile_moves_hder.close();
	}
	 strg.str(std::string());
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_moves.open(filename_moves.c_str(), std::ofstream::binary); // Overwrite any existing file in BINARY
	 } else {
		outfile_moves.open(filename_moves.c_str(), std::ofstream::app | std::ofstream::binary); // std::app is for append
	 }
	 if (outfile_moves.is_open()){
		for (int i=0; i<Ntot; i++){ // we write the buffer
			for(int j=0; j<Nchains; j++){ // The probabilities for each moves
				outfile_moves.write(reinterpret_cast<char*>(&buf_proposal.Pmoves(i,j)), size_dbl);
			}
			for(int j=0; j<Nchains; j++){ // Booleans telling whether we moved or not
				boolval=buf_proposal.moveds[i][j];
				outfile_moves.write(reinterpret_cast<char*>(&boolval), size_bool);
			}
		}
		outfile_moves.flush(); // Explicitly specify to flush the data into the disk
		outfile_moves.close();
  	}
  	else {
		std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << file_ext << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

	/////// Write the covarmats ////////
	if (need_header == 1){
		outfile_covarmats_hder.open(filename_covarmats_hder.c_str());  // Write the header ASCII file
		outfile_covarmats_hder << "# This is the header of the BINARY output file for the parameters of the proposal law. These may vary if the MALA algorithm is learning.\n";
		outfile_covarmats_hder << "# This file contains only values for covarmat[0:Nchains-1][ 0:Nvars-1][ 0:Nvars-1]. Each matrix is in a different file, indexed by the chain number\n" ;
		outfile_covarmats_hder << "! Nchains= " << Nchains << "\n";
		outfile_covarmats_hder << "! Nvars= " << Nvars << "\n";
		outfile_covarmats_hder << "! Nsamples_done=" << Nbuffer * (buf_proposal.Ncopy) + buf_proposal.counts + buf_restore.Nsamples_sofar << "\n";
		// ---- The variable names
		strg.str(std::string());
		strg << "! variable_names=";
		for (int i=0; i<Nvars; i++){
			strg << buf_params.vars_names[i] << "   ";
		}
		strg << "\n";
		outfile_covarmats_hder << strg.str().c_str();
		outfile_covarmats_hder.close();
	}
	for (chain=0; chain<Nchains; chain++){

		 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
			outfile_covarmats[chain]= new std::ofstream(filename_covarmats[chain].c_str(), std::ofstream::binary);  
		 } else {
			outfile_covarmats[chain]= new std::ofstream(filename_covarmats[chain].c_str(), std::ofstream::app | std::ofstream::binary); // std::app is for append
		 }
		 if (outfile_covarmats[chain]->is_open()){
			for (int i=0; i<Ntot; i++){ // we write the buffer
				for(ind_row=0; ind_row<Nvars; ind_row++){
					for(ind_col=0; ind_col<Nvars; ind_col++){
						outfile_covarmats[chain]->write( reinterpret_cast<char*>(& (*buf_proposal.covarmats[i][chain])(ind_row,ind_col) ), size_dbl );
					}
				}
			}
			outfile_covarmats[chain]->flush(); // Explicitly specify to flush the data into the disk
			outfile_covarmats[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << proposal_txtbin_fileout << "." << file_ext << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
	} // End of the loop on chain

	for( std::vector<std::ofstream*>::iterator it=outfile_mus.begin() ; it != outfile_mus.end(); ++it){
		delete (*it);
	}
	for( std::vector<std::ofstream*>::iterator it=outfile_covarmats.begin() ; it != outfile_covarmats.end(); ++it){
		delete (*it);
	}
	
}


//------------------
void Outputs::write_bin_params(const long Nrest){
	 
         double dblval=0.0;
	 size_t size_dbl=sizeof(dblval);

	 bool need_header=1; // Variable for handling the header (which is in ASCII)... this is the metadata.
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_vars_hder;
	 std::vector< std::string > filename_vars;
	 std::ofstream outfile_vars_hder;
	 std::vector<std::ofstream*> outfile_vars(Nchains);

	filename_vars_hder=params_txtbin_fileout + ".hdr";

	for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_vars.push_back(params_txtbin_fileout + "_chain-" + ind_str.str() + "." + file_ext) ;
		ind_str.str(std::string());
	 }

	if(erase_old_files == 0){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists
	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the header ////////
	if (need_header == 1){
		outfile_vars_hder.open(filename_vars_hder.c_str());  // Write the header ASCII file
	 	outfile_vars_hder << "# This is the header file of the BINARY output file for the model parameters \n";
		outfile_vars_hder << "# This file contains values for vars[0:Nchains-1][ 0:Nvars-1]. Each matrix is in a different file, indexed by the chain number\n" ;
		// ----			
		outfile_vars_hder << "! Nsamples= " << Nsamples << "\n";	
		outfile_vars_hder << "! Nchains= " << Nchains << "\n";
		outfile_vars_hder << "! Nsamples_done=" << Nbuffer * (buf_params.Ncopy) + buf_params.counts + buf_restore.Nsamples_sofar << "\n";
		outfile_vars_hder << "! Nvars= " << Nvars << "\n";
		outfile_vars_hder << "! Ncons= " << Ncons << "\n";
		outfile_vars_hder << "! relax= " << buf_params.relax.transpose() << "\n";
		// ---- The constant names
		strg.str(std::string());
		strg << "! plength= ";
		strg << buf_params.plength.transpose() << "\n";
		outfile_vars_hder << strg.str().c_str();
		// ---- The constant names
		strg.str(std::string());
		strg << "! constant_names= ";
		for (int i=0; i<buf_params.cons.size(); i++){
			strg << buf_params.cons_names[i] << "   ";
		}
		strg << "\n";
		outfile_vars_hder << strg.str().c_str();
		// ---- The constant values
		strg.str(std::string());
		strg << "! constant_values= ";
		if ( buf_params.cons_names[0] == "None"){
			strg << "-1\n";
		} else {
			strg << buf_params.cons.transpose() << "\n";
		}
		outfile_vars_hder << strg.str().c_str();
		// ---- The variable names
		strg.str(std::string());
		strg << "! variable_names=";
		for (int i=0; i<Nvars; i++){
			strg << buf_params.vars_names[i] << "   ";
		}
		strg << "\n";
		outfile_vars_hder << strg.str().c_str();
		outfile_vars_hder.close();
		}
	/////// Write the parameters ////////
	for (chain=0; chain<Nchains; chain++){
		 if (erase_old_files == 1 && buf_params.Ncopy == 0) {
			outfile_vars[chain]=new std::ofstream((filename_vars[chain]).c_str(), std::ofstream::binary);
		 } else {
			outfile_vars[chain]=new std::ofstream((filename_vars[chain]).c_str(), std::ofstream::app | std::ofstream::binary);
		 }
		 if (outfile_vars[chain]->is_open()){
			for (int i=0; i<Ntot; i++){ // we write the buffer	
				for(int j=0; j< Nvars; j++){
					outfile_vars[chain]->write( reinterpret_cast<char*>(& (*buf_params.vars[i])(chain,j) ), size_dbl );
				}
			}
			outfile_vars[chain]->flush(); // Explicitly specify to flush the data into the disk
			strg.str(std::string());
			
			outfile_vars[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << params_txtbin_fileout << "_[chain]." + file_ext << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}

	} // End of the loop on chain

	for( std::vector<std::ofstream*>::iterator it=outfile_vars.begin() ; it != outfile_vars.end(); ++it){
		delete (*it);
	}
}

//------------------
void Outputs::write_bin_parallel_temp_params(const long Nrest){

	bool boolval=0;
	int intval=0;
	double dblval=0.0;
	size_t size_bool=sizeof(boolval);
	size_t size_int=sizeof(intval);
	size_t size_dbl=sizeof(dblval);

	 bool need_header=1; // Variable for handling the header (which is in ASCII)... this is the metadata.
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_parallel_temp, filename_parallel_temp_hder;
	 std::ofstream outfile_parallel_temp, outfile_parallel_temp_hder;

 	filename_parallel_temp=parallel_tempering_txtbin_fileout + "." + file_ext;
	filename_parallel_temp_hder=parallel_tempering_txtbin_fileout + ".hdr";

	if(erase_old_files == 0){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists
	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}

	/////// Write the Header in a separate file /////////
	if (need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
		outfile_parallel_temp_hder.open(filename_parallel_temp_hder.c_str());  // Write the header ASCII file
   		outfile_parallel_temp_hder << "# This is the header of the BINARY output file for the parameters of the parallel tempering.\n";
		outfile_parallel_temp_hder << "# This file contains values for \n" ;
		outfile_parallel_temp_hder << "# Correspondance between chain0=[0:Nchains-1] and temperature Tcoefs[chain] \n" ;
		outfile_parallel_temp_hder << "! Nsamples_done=" << Nbuffer * (buf_parallel_temp.Ncopy) + buf_parallel_temp.counts + buf_restore.Nsamples_sofar << "\n";
		strg.str(std::string());
		strg << "! Tcoefs = ";
		strg << buf_parallel_temp.Tcoefs.transpose() << "\n";
		outfile_parallel_temp_hder << strg.str().c_str();
		outfile_parallel_temp_hder << "! labels= attempt_mixing    chain0    Pswitch    switched \n" ;
		outfile_parallel_temp_hder.close();
	}

	/////// Write the parameters of the parallel tempering ////////
	 if (erase_old_files == 1 && buf_proposal.Ncopy == 0) {
		outfile_parallel_temp.open(filename_parallel_temp.c_str(), std::ofstream::binary); // Overwrite any existing file in BINARY
	 } else {
		outfile_parallel_temp.open(filename_parallel_temp.c_str(), std::ofstream::app | std::ofstream::binary); // std::app is for append
	 }
	 if (outfile_parallel_temp.is_open()){
		for (int i=0; i<Ntot; i++){ // we write the buffer
			boolval=buf_parallel_temp.attempt_mixing[i];
			outfile_parallel_temp.write(reinterpret_cast<char*>(&boolval), size_bool); // Write first the attempt_mixing variable (booleans)
			outfile_parallel_temp.write(reinterpret_cast<char*>(&buf_parallel_temp.chain0s[i]), size_int); // then the chain0s variable (integers)
			dblval=buf_parallel_temp.Pswitchs[i];
			outfile_parallel_temp.write(reinterpret_cast<char*>(&dblval), size_dbl); // them the Pswitchs (doubles)
			boolval=buf_parallel_temp.switcheds[i];
			outfile_parallel_temp.write(reinterpret_cast<char*>(&boolval), size_bool); // finaly write the switcheds (booleans)
		}
		outfile_parallel_temp.flush(); // Explicitly specify to flush the data into the disk
		//strg.str(std::string()); // clear the strg buffer
		outfile_parallel_temp.close();
  	}
  	else {
		std::cout << " Unable to open file " << parallel_tempering_txtbin_fileout << "." + file_ext << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}
}


//------------------
void Outputs::write_bin_models(const long Nrest){

	 double dblval=0.0;
	 size_t size_dbl=sizeof(dblval);

	 bool need_header=1; // Variable for handling the header (which is in ASCII)... this is the metadata.
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_models_hder;
	 std::vector< std::string > filename_models;
	 std::ofstream outfile_models_hder;
	 std::vector<std::ofstream*> outfile_models(Nchains);

	 for(chain=0; chain<Nchains; chain++){
		ind_str << chain;
		filename_models.push_back(model_txtbin_fileout + "_models_chain-" + ind_str.str() + "." + file_ext) ;
		ind_str.str(std::string());
	 }

	if(file_exists(filename_models[0].c_str()) && (erase_old_files == 1)){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}
	/////// Write the Header in a separate file ////////
	if (buf_models.Ncopy == 0 && need_header == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
		outfile_models_hder.open(filename_models_hder.c_str());  // Write the header ASCII file
		outfile_models_hder << "# This is a header for the BINARY output file for the models. \n";
		outfile_models_hder << "! Nchains= " << Nchains << "\n";
		outfile_models_hder << "! Nvars= " << Nvars << "\n";
		outfile_models_hder << "! chain= " << chain << "\n";
		outfile_models_hder.close();
	}

	/////// Write the models ////////
	for (chain=0; chain<Nchains; chain++){
		 	if (erase_old_files == 1 && buf_models.Ncopy == 0) {
				outfile_models[chain]= new std::ofstream(filename_models[chain].c_str(), std::ofstream::binary); // Overwrite any existing file
			 } else {
				outfile_models[chain]= new std::ofstream(filename_models[chain].c_str(), std::ofstream::app | std::ofstream::binary); // std::app is for append
		 	}
		 //} // End of the if with Ncopy
		 if (outfile_models[chain]->is_open()){
			for (int i=0; i<Ntot; i++){ // we write the buffer	
				for(int j=0; j<Ndata_obs; j++){
					outfile_models[chain]->write( reinterpret_cast<char*>(& (*buf_models.models[i])(chain,j) ), size_dbl );
				}
			}
			outfile_models[chain]->flush(); // Explicitly specify to flush the data into the disk
			outfile_models[chain]->close();
  		} // End of the if with is_open()
  		else {
			std::cout << " Unable to open file " << model_txtbin_fileout << "." + file_ext << std::endl;	
			std::cout << " Check that the full path exists" << std::endl;
			std::cout << " The program will exit now" << std::endl;
			exit(EXIT_FAILURE);
		}
	} // End of the loop on chain

	for( std::vector<std::ofstream*>::iterator it=outfile_models.begin() ; it != outfile_models.end(); ++it){
		delete (*it);
	}
}

void Outputs::write_bin_stat_criteria(const long Nrest){

	 double dblval=0.0;
	 size_t size_dbl=sizeof(dblval);

	 bool need_header=1; // Variable for handling the header (which is in ASCII)... this is the metadata.
	 int chain, ind_row, index, Ntot;
	 std::ostringstream strg, ind_str;
	 std::string filename_statcrit, filename_statcrit_hder;
	 std::ofstream outfile_statcrit, outfile_statcrit_hder;
	 std::vector<std::string> labels(3);

	labels[0]= "logLikelihood";
	labels[1]= "logPrior";
	labels[2]= "logPosteriors";
	filename_statcrit=stat_criteria_txtbin_fileout + "." + file_ext ;
	filename_statcrit_hder=stat_criteria_txtbin_fileout + ".hdr"  ;

	if(erase_old_files == 0){need_header=0;} // In case of an append, we do not need of a header only if the file actually exists

	if (Nrest < Nbuffer){ // If the number of remaining samples is smaller than the buffer, write only what remains
		Ntot=Nrest;
	} else{
		Ntot=Nbuffer;
	}
	/////// Write the Header in a separate file ////////
	if (need_header == 1){
		outfile_statcrit_hder.open(filename_statcrit_hder.c_str());  // Write the header ASCII file
    		outfile_statcrit_hder << "# This is the header of the BINARY output file for the statistical information.\n";
		outfile_statcrit_hder << "# This file contains values for the logLikelihood (columns 0:Nchains-1), logPrior (columns Nchains:2*Nchains-1) and logPosterior (columns 2*Nchains:3*Nchains-1),  \n" ;
		outfile_statcrit_hder << "! Nsamples_done=" << Nbuffer * (buf_stat_crit.Ncopy) + buf_stat_crit.counts + buf_restore.Nsamples_sofar << "\n";
		outfile_statcrit_hder << "! Nchains= " << Nchains << "\n";
		outfile_statcrit_hder << "! labels= ";
		for (int k=0; k<labels.size();k++){
			for (int i=0; i<Nchains; i++){			
				outfile_statcrit_hder << labels[k] << "[" << i << "]   ";
			}
		}
		outfile_statcrit_hder << "\n";
		outfile_statcrit_hder.close();
	}

	/////// Write the statistical criteria (logLikelihood, logPrior, logPosterior) ////////
	 if (erase_old_files == 1 && buf_stat_crit.Ncopy == 0) {
		outfile_statcrit.open(filename_statcrit.c_str(), std::ofstream::binary); // Overwrite any existing file
	 } else {
		outfile_statcrit.open(filename_statcrit.c_str(), std::ofstream::app | std::ofstream::binary); // std::app is for append
	 }
	 if (outfile_statcrit.is_open()){
		for (int i=0; i<Ntot; i++){ // we write the buffer
			for(chain=0; chain<Nchains; chain++){
				dblval=buf_stat_crit.logLikelihoods(i,chain);
				outfile_statcrit.write(reinterpret_cast<char*>(&dblval), size_dbl);
			}
			for(chain=0; chain<Nchains; chain++){
				dblval=buf_stat_crit.logPriors(i,chain);
				outfile_statcrit.write(reinterpret_cast<char*>(&dblval), size_dbl);
			}			
			for(chain=0; chain<Nchains; chain++){
				dblval=buf_stat_crit.logPosteriors(i,chain);
				outfile_statcrit.write(reinterpret_cast<char*>(&dblval), size_dbl);
			}
			
		}
		outfile_statcrit.flush(); // Explicitly specify to flush the data into the disk
		
		outfile_statcrit.close();
  	}
  	else {
		std::cout << " Unable to open file " << stat_criteria_txtbin_fileout << "." + file_ext << std::endl;	
		std::cout << " Check that the full path exists" << std::endl;
		std::cout << " The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	}

}
//------------------
//------------------
//------------------

void Outputs::update_buffer_stat_criteria(const VectorXd& logLikelihood, const VectorXd& logPrior, const VectorXd& logPosterior){

 if(get_statcriteria == 1) {
	long rest=Nsamples - buf_restore.Nsamples_sofar - Nbuffer*buf_stat_crit.Ncopy;
	if ((buf_stat_crit.counts < Nbuffer) && (buf_stat_crit.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		buf_stat_crit.logLikelihoods.row(buf_stat_crit.counts)=logLikelihood;
		buf_stat_crit.logPriors.row(buf_stat_crit.counts)=logPrior;
		buf_stat_crit.logPosteriors.row(buf_stat_crit.counts)=logPosterior;
		buf_stat_crit.counts=buf_stat_crit.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text"){
			write_txt_stat_criteria(rest, 0);
			std::cout << "     The output file (ASCII) with Statistical information has been updated" << std::endl;
		}		
		if (file_format == "binary"){
			write_bin_stat_criteria(rest);
			std::cout << "     The output file (BINARY) with Statistical information has been updated" << std::endl;
		}
		if (file_format == "debug"){
			write_bin_stat_criteria(rest);
			write_txt_stat_criteria(rest, 1);
			std::cout << "    Debug Mode: The output file (ASCII + BINARY) with Statistical information has been updated" << std::endl;
		}
		buf_stat_crit.counts=0;
		buf_stat_crit.Ncopy=buf_stat_crit.Ncopy +1;

		// Then we write the new entry on the empty buffer
		buf_stat_crit.logLikelihoods.row(buf_stat_crit.counts)=logLikelihood;
		buf_stat_crit.logPriors.row(buf_stat_crit.counts)=logPrior;
		buf_stat_crit.logPosteriors.row(buf_stat_crit.counts)=logPosterior;

		buf_stat_crit.counts=buf_stat_crit.counts +1;

	}
  }
}

void Outputs::update_buffer_models(const MatrixXd& models_in){

 if(get_models == 1) {
	long rest=Nsamples - buf_restore.Nsamples_sofar - Nbuffer*buf_models.Ncopy;
	if ((buf_models.counts < Nbuffer) && (buf_models.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		set_3dMatrix(buf_models.models, models_in, buf_models.counts); // (Nbuffer, Nchains, Ndata_obs);
		buf_models.counts=buf_models.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text") {
			write_txt_models(rest, 0);
			std::cout << "     The output file (ASCII) with model information has been updated" << std::endl;
		}
		if (file_format == "binary") {
			write_bin_models(rest);
			std::cout << "     The output file (BINARY) with model information has been updated" << std::endl;
		}
		if (file_format == "debug") {
			write_bin_models(rest);
			write_txt_models(rest, 1);
			std::cout << "     DEBUG MODE: The output file (ASCII + BINARY) with model information has been updated" << std::endl;
		}

		buf_models.counts=0;
		buf_models.Ncopy=buf_models.Ncopy +1;

		set_3dMatrix(buf_models.models, models_in, buf_models.counts); // (Nbuffer, Nchains, Ndata_obs);
		buf_models.counts=buf_models.counts +1; 
	}
  }
}

void Outputs::update_buffer_params(const MatrixXd& params_in){

 if(get_params == 1) {
	long rest=Nsamples - buf_restore.Nsamples_sofar - Nbuffer*buf_params.Ncopy;
	if ((buf_params.counts < Nbuffer) && (buf_params.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		set_3dMatrix(buf_params.vars, params_in, buf_params.counts); // (Nbuffer, Nchains, Nvars);
		buf_params.counts=buf_params.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text") {
			write_txt_params(rest, 0);
			std::cout << "     The output file (ASCII) with parameters information has been updated" << std::endl;
		}
		if (file_format == "binary") {
			write_bin_params(rest);
			std::cout << "     The output file (BINARY) with parameters information has been updated" << std::endl;
		}
		if (file_format == "debug") {
			write_bin_params(rest);
			write_txt_params(rest, 1);
			std::cout << "     Debug Mode: The output file (ASCII + BINARY) with parameters information has been updated" << std::endl;
		}
		buf_params.counts=0;
		buf_params.Ncopy=buf_params.Ncopy +1;
		

		set_3dMatrix(buf_params.vars, params_in, buf_params.counts); // (Nbuffer, Nchains, Nvars);
		buf_params.counts=buf_params.counts +1; 

	}
  }
}

void Outputs::update_buffer_ptempering(const bool tempted_mixing, const int chain_A, const long double Probaswitch, const bool bool_switched){

 if(get_parallel_tempering_params == 1) {
	long rest=Nsamples - buf_restore.Nsamples_sofar - Nbuffer*buf_parallel_temp.Ncopy;
	if ((buf_parallel_temp.counts < Nbuffer) && (buf_parallel_temp.counts != (rest - 1))){ // We do not write on file only if the counts < Nbuffer and if the remaining number of samples is not reached 
		buf_parallel_temp.attempt_mixing[buf_parallel_temp.counts]=tempted_mixing; // Tells us if we actually tried to mix the chains (should be the term i%dN_mixing == 0 in MALA::execute())
		buf_parallel_temp.chain0s[buf_parallel_temp.counts]=chain_A;
		buf_parallel_temp.Pswitchs[buf_parallel_temp.counts]=Probaswitch; // Probability of switching chain_A
		buf_parallel_temp.switcheds[buf_parallel_temp.counts]=bool_switched; // Did we switched the chain_A with the chain_A+1 ?
	
		buf_parallel_temp.counts=buf_parallel_temp.counts +1; 
	} else { 
		// We write the data and reinitialize the counts to 0
		if (file_format == "text") {
			write_txt_parallel_temp_params(rest, 0);
			std::cout << "     The output file (ASCII) with information on parallel tempering parameters has been updated" << std::endl;
		}
		if (file_format == "binary") {
			write_bin_parallel_temp_params(rest);
			std::cout << "     The output file (BINARY) with information on parallel tempering parameters has been updated" << std::endl;
		}
		if (file_format == "debug") {
			write_bin_parallel_temp_params(rest);
			std::cout << "   Debug Mode:  The output file (ASCII + BINARY) with information on parallel tempering parameters has been updated" << std::endl;
			write_txt_parallel_temp_params(rest, 1);
		}
		buf_parallel_temp.counts=0;
		buf_parallel_temp.Ncopy=buf_parallel_temp.Ncopy +1;
		

		buf_parallel_temp.attempt_mixing[buf_parallel_temp.counts]=tempted_mixing; // Tells us if we actually tried to mix the chains (should be the term i%dN_mixing == 0 in MALA::execute())
		buf_parallel_temp.chain0s[buf_parallel_temp.counts]=chain_A;
		buf_parallel_temp.Pswitchs[buf_parallel_temp.counts]=Probaswitch; // Probability of switching chain_A
		buf_parallel_temp.switcheds[buf_parallel_temp.counts]=bool_switched; // Did we switched the chain_A with the chain_A+1 ?
	
		buf_parallel_temp.counts=buf_parallel_temp.counts +1;
	}
	
  }
}

void Outputs::update_buffer_proposals(const VectorXd& sigma_chains, const MatrixXd& mu_chains, const VectorXd& Pmove_chains, 
			const std::vector<bool> moved_chains, MatrixXd **covarmat_chains){

  MatrixXd acceptance;
  long rest=Nsamples - buf_restore.Nsamples_sofar - Nbuffer*buf_proposal.Ncopy;

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
		if(get_proposal_params == 1) {
			if(file_format == "text"){
				write_txt_prop_params(rest, 0);
				std::cout << "     The output file (ASCII) with information on proposal parameters has been updated" << std::endl;
			}
			if(file_format == "binary"){
				write_bin_prop_params(rest);
				std::cout << "     The output file (BINARY) with information on proposal parameters has been updated" << std::endl;
			}
			if(file_format == "debug"){
				write_bin_prop_params(rest);
				write_txt_prop_params(rest, 1);
				std::cout << "    Debug Mode: The output file (ASCII + BINARY) with information on proposal parameters has been updated" << std::endl;
			}
		}
		buf_proposal.counts=0;
		buf_proposal.Ncopy=buf_proposal.Ncopy +1;
		
		buf_proposal.sigmas.row(buf_proposal.counts)=sigma_chains;
		set_3dMatrix(buf_proposal.mus, mu_chains, buf_proposal.counts);
		buf_proposal.Pmoves.row(buf_proposal.counts)=Pmove_chains;		

		for(int m=0; m<Nchains; m++){
			buf_proposal.moveds.at(buf_proposal.counts)[m]=moved_chains[m]; // 2D vector of booleans. can be called by doing something like this: resize[counts][m]
			set_4dMatrix(buf_proposal.covarmats, *covarmat_chains[m], buf_proposal.counts, m); // can be called by doing something like this: covarmats[counts][chain]
		}

		// We calculate the acceptance rate
		acceptance=reject_rate(buf_proposal.moveds, Nbuffer, buf_acceptance.Ncopy, Nchains, buf_restore.Nsamples_sofar);
		buf_acceptance.xaxis=acceptance(0,0); // other values are irrelevant because these are all the same
		buf_acceptance.acceptance_rate=acceptance.row(1);
		// We write the acceptance rate on a file
		write_txt_acceptance();
		std::cout << "     The output file (ASCII) with information on the acceptance rate has been updated" << std::endl;
		buf_acceptance.Ncopy=buf_acceptance.Ncopy +1;
		
		buf_proposal.counts=buf_proposal.counts +1; 
	}
  //}

}


void Outputs::update_buffer_restore(const VectorXd& sigma_chains, const MatrixXd& mu_chains, MatrixXd **covarmat_chains, const MatrixXd& params_in){
/* 
 * Handle the output used to restore a job which has finished (e.g. end of a Burn-in phase) or that has failed for some reason (computer failure)
*/
  long rest=Nsamples - buf_restore.Nsamples_sofar - Nbuffer*buf_restore.Ncopy;

  if ((buf_restore.counts == Nbuffer) || (buf_restore.counts == (rest - 1))){ // We write on file only if the i is modulo Nbuffer OR if the remaining number of samples is reached 


		// We write the data and reinitialize the counts to 0
		buf_restore.sigmas=sigma_chains;
		buf_restore.mus=mu_chains;
		buf_restore.covarmats=covarmat_chains;
		buf_restore.vars=params_in;
		
		// Compute the mean by (1) adding the last samples and (2) normalizing by the number of computed samples since the last record (Nbuffer or rest)
		buf_restore.sigmas_mean=(buf_restore.sigmas_mean + sigma_chains)/buf_restore.counts;
		buf_restore.mus_mean=(buf_restore.mus_mean + mu_chains)/buf_restore.counts;
		for (int chain=0; chain<Nchains; chain++){
			//buf_restore.covarmats_mean[chain]=(buf_restore.covarmats_mean[chain] + covarmat_chains[chain])/Nbuffer;
			for (int i=0; i<Nvars; i++){ // we write the covariance matrix	
				buf_restore.covarmats_mean[chain]->row(i) = (buf_restore.covarmats_mean[chain]->row(i) + covarmat_chains[chain]->row(i))/buf_restore.counts;;
			}
		}
		buf_restore.vars_mean=(buf_restore.vars_mean + params_in)/buf_restore.counts;
		

		// Write the buffer using the dedicated function
		write_buffer_restore();
	
		// Initialize the variables with the mean values
		buf_restore.sigmas_mean.setZero();
		buf_restore.mus_mean.setZero();
		for (int chain=0; chain<Nchains; chain++){
			buf_restore.covarmats_mean[chain]->setZero();
		}
		buf_restore.vars_mean.setZero();
		
		// Initialize counters
		buf_restore.counts=0;
		buf_restore.Ncopy=buf_restore.Ncopy +1;
		std::cout << "     The restoration file (ASCII) has been updated" << std::endl;
	
		buf_restore.counts=buf_restore.counts +1; 
   } else {
   		// Do the Sum that enables the computation of the averaged values
   		buf_restore.sigmas_mean=buf_restore.sigmas_mean + sigma_chains;
		buf_restore.mus_mean=buf_restore.mus_mean + mu_chains;

		for (int chain=0; chain<Nchains; chain++){
			for (int i=0; i<Nvars; i++){ // we write the covariance matrix	
				buf_restore.covarmats_mean[chain]->row(i) = (buf_restore.covarmats_mean[chain]->row(i) + covarmat_chains[chain]->row(i));
			}
		}
		buf_restore.vars_mean=buf_restore.vars_mean + params_in;

		buf_restore.counts=buf_restore.counts +1; 
	}
 
}

///////////////// Extra functions /////////////////
bool Outputs::file_exists(const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

MatrixXd Outputs::reject_rate(const std::vector< std::vector<bool> > moveds, const long Nbuffer, const long Ncopy, const int Nchains, const long Nsamples_init){
/*
 * Takes the moves (saved in the buffer for proposals) and calculates 
 * the acceptance rates for each chain for a given slice of size Nbuffer. 
 * Nbuffer and Ncopy are used to create a relevant x-axis.
*/

	long Ntrue;
	MatrixXd acceptance(2, Nchains);

	for (int m=0; m<Nchains; m++){
		//std::cout << "chain = " << m << std::endl;
		Ntrue=count_accepted_vals(moveds, m, Nbuffer);
		acceptance(1, m)=double(Ntrue)/double(Nbuffer);
		acceptance(0, m)=(double(Ncopy)+0.5)*double(Nbuffer) + Nsamples_init; // Average index of the averaged zone for the acceptance
	}

return acceptance;
}
 
long Outputs::count_accepted_vals(const std::vector< std::vector<bool> > moves_chain, const long m, const long Nbuffer){
/*
 * One 1D bool array, counts the number of states which are at True
*/
	long count=0;

	for(long i=0; i<Nbuffer; i++){
		if (moves_chain.at(i)[m] == 1){
			count=count+1;
		}
	}
	
	//std::cout << "count = " << count << std::endl;
	return count;
}

int Outputs::get_Nbuffer(){
	return Nbuffer;
}


