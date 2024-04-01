/*
 * outputs_targz.cpp
 *
 * Contains the class and all kind of functions
 * used to write the data on files that are compressed in tar.gz
 * 
 *  Created on: 02 Nov 2023
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <zlib.h>
#include <tar.h>

//void Outputs::write_bintar_prop_params(const long Nrest);

bool compressFile(const std::string& inputFile, const std::string& outputFile) {
    gzFile outFile = gzopen(outputFile.c_str(), "wb");
    if (!outFile) {
        std::cerr << "Failed to open output file: " << outputFile << std::endl;
        return false;
    }

    std::ifstream inFile(inputFile, std::ios::binary);
    if (!inFile) {
        std::cerr << "Failed to open input file: " << inputFile << std::endl;
        gzclose(outFile);
        return false;
    }

    char buffer[1024];
    while (inFile.good()) {
        inFile.read(buffer, sizeof(buffer));
        gzwrite(outFile, buffer, inFile.gcount());
    }

    inFile.close();
    gzclose(outFile);

    return true;
}


void Outputs::write_bintar_prop_params(const long Nrest){

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

