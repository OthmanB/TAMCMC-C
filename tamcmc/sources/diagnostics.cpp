/*  diagnostics.cpp
 *
 *  Contains procedure to handle graphics (using gnuplot) and other diagnostics outputs
 * 
 *  Created on: 20 Apr 2016
 *      Author: obenomar
 */


#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <stdio.h>
#include "diagnostics.h"

Diagnostics::Diagnostics(Config *cfg){
/*
 *
 * Initialise the configuration for the diagnostics by importing relevant variable into the diagnostic class
 *
*/

    std::stringstream ss;

	ps=1;
	Nchains=cfg->MALA.Nchains;
	Nbuffer=cfg->outputs.Nbuffer;
	Nsamples=cfg->outputs.Nsamples;
	file_out_format=cfg->outputs.file_out_format;

	if (file_out_format== "text"){
		filename_likelihood=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.stat_txt_fileout + ".txt";
		filename_acceptance=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.acceptance_txt_fileout + ".txt";
		filename_params=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.params_txt_fileout + "_chain-0.txt"; // We only care about the coldest chain

		filename_likelihood_hdr="";
		filename_params_hdr="";
	} 
	if (file_out_format== "binary"){
		filename_likelihood=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.stat_txt_fileout + ".bin";
		filename_acceptance=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.acceptance_txt_fileout + ".txt"; // The acceptance is still written in ASCII
		filename_params=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.params_txt_fileout + "_chain-0.bin"; // We only care about the coldest chain

		filename_likelihood_hdr=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.stat_txt_fileout + ".hdr";
		filename_params_hdr=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.params_txt_fileout + ".hdr"; // The header of the parameters
	}
	if (file_out_format== "debug"){
		filename_likelihood=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.stat_txt_fileout + ".dbg";
		filename_acceptance=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.acceptance_txt_fileout + ".txt"; // The acceptance is still written in ASCII
		filename_params=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.params_txt_fileout + "_chain-0.dbg"; // We only care about the coldest chain

		filename_likelihood_hdr=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.stat_txt_fileout + ".hdr";
		filename_params_hdr=cfg->outputs.dir_out + cfg->outputs.output_root_name + cfg->outputs.params_txt_fileout + ".hdr"; // The header of the parameters
	}
	if (file_out_format != "text" && file_out_format != "binary" && file_out_format != "debug" ) {
                std::cout << " Other format than text/binary for the likelihood / acceptance is not supported " << std::endl;
                std::cout << " Please set file_format either to 'text' 'binary' or 'debug' into the config.cfg file" << std::endl;
		std::cout << " The program will exit now " << std::endl;
		exit(EXIT_FAILURE);
	}
	chains_diags=cfg->diags.chains_diags;;
	evidence_diags=cfg->diags.evidence_diags;
	pdfs_diags=cfg->diags.pdfs_diags;
	
	output_dir=cfg->diags.output_dir;
	output_root_name=cfg->diags.output_root_name;

	file_chains_diags=output_dir + output_root_name + cfg->diags.file_chains_diags;
	file_evidence_diags=output_dir + output_root_name +  cfg->diags.file_evidence_diags;
	file_pdfs_diags=output_dir +  output_root_name + cfg->diags.file_pdfs_diags; 

	model_initial_diags=cfg->diags.model_initial_diags;
	model_buffer_diags=cfg->diags.model_buffer_diags;
	model_final_diags=cfg->diags.model_final_diags;

	file_model_init_diags=output_dir +  output_root_name + cfg->diags.file_model_init_diags;
	file_model_buffer_diags=output_dir +  output_root_name + cfg->diags.file_model_buffer_diags;
	file_model_final_diags=output_dir +  output_root_name + cfg->diags.file_model_final_diags;

	data_scoef1=cfg->diags.data_scoef1; 
	data_scoef2=cfg->diags.data_scoef2;
	show_original_data=cfg->diags.show_original_data;
	Nclasses=cfg->diags.Nclasses;

//BENOIT: make temp file unique per running program (ie enable many copies of program to run concurrently)
    	ss.str(std::string());
    	ss << getpid();
        file_data_tmp="tmp/graph_tmp_"+ss.str()+".txt";
        file_data_likelihood_tmp="tmp/likelihood_diag_" + ss.str() +".txt";
        file_hist_tmp="tmp/hist_tmp_" + ss.str() + ".txt";

        evidence_interpolation_factor=cfg->diags.evidence_interpolation_factor;

        firstpass=1;
}

Diagnostics::Diagnostics(){


}

void Diagnostics::gnuplt_chains_diags(int i, const VectorXd& Tcoefs){
/* 
 * Simple function that handle the plots of the likelihood and of the acceptance rate.
 * Here the data are read in two files. These are assumed to exist a priori.
 * If ps is set to 1, then the plot is written in file_out. Otherwise, it is shown
 * on the screen.
 * I use disk-file temporary files to use gnuplot because memory-based calls of gnuplot
 * happenned to be extremely slow/inneficient (surprisingly)
*/

    int m;
    MatrixXd Mat, Mattmp;
    std::string m_str, m2_str, tmp_str, filetmp;
    std::stringstream ss;
    Evidence_out evidence;

      if(chains_diags == 1 && ((i%Nbuffer == 0 && i !=0) || i == Nsamples-1)){
    	m=0; // initialize the index of the chain that will be shown for the likelihood

    	m_str=int_to_str(m+1);
    	m2_str=int_to_str(m);

    	Gnuplot gn;
	
    	if (ps == 0) { 
		gn << "set term X11 \n";
    	} else {
		gn << "set term post eps enhanced color font 'Times-Bold, 15' size 29.7cm,21cm \n";
		gn << "set out '" + file_chains_diags + ".eps'\n";
   	 }
   	 // Setup the common plot configuration
 	   gn << "set datafile commentschars '#!'\n";
 	   gn << "set autoscale \n"; // scale axis automatically
 	   gn << "unset label \n";
 	   gn << "set key default \n"; // restore the default position for the legends
 	   gn << "set key bottom right \n"; // legend will be on bottom right

  	  // The plot for the likelihood of all parallel chains
 	   gn << "set xlabel 'iteration'" << std::endl;
 	   gn << "set ylabel 'Likelihood[" + m_str + "]'" << std::endl;
	   gn << "set xrange [0:" + int_to_str(Nsamples) + "] \n";
	   if (file_out_format == "binary" || file_out_format == "debug"){ // BINARY FILE CASE
                filetmp=file_hist_tmp;
                Mat=read_bin_matrix_dbl(filename_likelihood, 3*Nchains, Nsamples, "dbl"); // Read the Binary file
                write_txt_matrix_dbl(Mat.topLeftCorner(Mat.rows(), Nchains), filetmp); // Write on a temporary ASCII file with only the Likelihoods (ignore the Priors and Posteriors)
                gn << "plot '" + filetmp +  "' using 0:" + m_str + " with lines title 'm=" + m2_str +  "'" ; // Plot
	   } else{ // PLAIN ASCII CASE
		filetmp=filename_likelihood; // Directly read the ASCII outputs
		gn << "plot '" + filetmp +  "' using 0:" + m_str + " with lines title 'm=" + m2_str +  "'" ;
	   }

 	   for (m=1; m< Nchains; m++){
 	  	m_str=int_to_str(m+1);
 	   	m2_str=int_to_str(m);
		tmp_str= tmp_str + ", "; 
		tmp_str= tmp_str + "'" + filetmp +  "' using 0:" + m_str + " with lines title 'm=" + m2_str + "'";
  	  }
  	  gn << tmp_str << std::endl;

 	   // On a split screen, we show the coolest chain and the acceptance rate for all chains
	   m=0; // in the split screen, only the coolest chain is shown...
	   m_str=int_to_str(m+1);
	   m2_str=int_to_str(m);

	   gn << "set multiplot\n"; // activate multiplt with two zones
	   gn << "set lmargin at screen 0.06\n";
	   gn << "set rmargin at screen 0.90\n";
	   
         // Top Area
	 gn << "set xrange [0:" + int_to_str(Nsamples) + "] \n";
	 if (ps == 0){
		gn << "set origin 0, 0.5" << std::endl;
		gn << "set size 1, 0.5" << std::endl;
                gn << "set key bottom right \n"; // legend will be on bottom right
  	  } else{
                gn << "set tmargin at screen 0.90\n";
                gn << "set bmargin at screen 0.52\n";
  	  }

         gn << "set xlabel 'iteration'" << std::endl;
         gn << "set ylabel 'Likelihood[" + m_str + "]'" << std::endl;
         gn << "plot '" + filetmp +  "' using 0:" + m_str + " with lines title 'm=" + m2_str +  "' \n" ;

  	 // Bottom Area
 	 gn << "set xrange [0:" + int_to_str(Nsamples) + "] \n";
  	 if (ps == 0){
		gn << "set origin 0, 0" << std::endl;
		gn << "set size 1, 0.5" << std::endl;
		gn << "set key top right \n"; // legend will be on bottom right
  	 } else{
                gn << "set tmargin at screen 0.48\n";
                gn << "set bmargin at screen 0.06\n";
	 }

 	 gn << "set xlabel 'iteration'" << std::endl;
 	 gn << "set ylabel 'Acceptance[" + m_str + "]'" << std::endl;
	
	 m=0;
	 m_str=int_to_str(m + 2);
	 m2_str=int_to_str(m);
  	 tmp_str="plot '" + filename_acceptance +  "' using 1:" + m_str + " with linespoints pointtype 5 title 'm=" + m2_str + "' " ;
 	 for (int m=1; m< Nchains; m++){
		m_str=int_to_str(m+2);
		m2_str=int_to_str(m);
		tmp_str= tmp_str + ", ";
		tmp_str= tmp_str + "'" + filename_acceptance +  "' using 1:" + m_str + " with linespoints pointtype 5 title 'm=" + m2_str + "'";
  	 } 
 	 gn << tmp_str << std::endl;

  	 gn << "unset multiplot" << std::endl; // come back to single zone

	std::cout << "     The output file with information on the acceptance rate has been updated" << std::endl;
	std::cout << "         - Likelihoods are read from the file: " << filename_likelihood << std::endl;
	std::cout << "         - Acceptance rates are read from the file: " << filename_acceptance << std::endl;

        // Compute the evidence if requested
	if(evidence_diags == 1){
		evidence=evidence_calc(Tcoefs, Mat.topLeftCorner(Mat.rows(), Nchains));
                write_evidence(evidence, i); // Write on a file the evidence as it is calculated:
    	}
	
        firstpass=0;
    /*if(file_out_format == "binary" || file_out_format == "debug") {  // This part of the code does not work... temporary files are not erased
		//std::cout << "In file: " << filetmp << std::endl;
		shell_exec("rm " + filetmp); // erase the temporary file
	}
	*/

	} // endif
	
}

Diagnostics::~Diagnostics(){ // The Destructor

	// Nothing so far because no pointer declared
}

//------------------


void Diagnostics::gnuplt_model_diags(Data *data_struc, const VectorXd& model, std::string phase){
/*
 * Used to plot the model on the top of the data using directly the
 * C++ variables that contain the above-mentionned quantities.
 * scoef1 and scoef2 corresponds to the smooth coefficients in natural units (e.g. microHz)
*/
	if(model_initial_diags == 1 && phase == "init"){
		std::cout << "    - Generating an initial plot of the data and of the model...";
		gnuplt_model_core(data_struc, model, phase);
		std::cout << " Done" << std::endl;
	}
	if(model_buffer_diags == 1 && phase == "buffer"){
		std::cout << "    - Generating a plot using the samples from the last Buffer for the data and for the model...";
		gnuplt_model_core(data_struc, model, phase);
		std::cout << " Done" << std::endl;
	}
	if(model_final_diags == 1 && phase == "final"){
		std::cout << "    - Generating the final plot of the data and of the model using all samples...";
		gnuplt_model_core(data_struc, model, phase);
		std::cout << " Done" << std::endl;
	}
}

void Diagnostics::gnuplt_model_core(Data *data_struc, const VectorXd& model, std::string phase){
/* 
 * Simple function that handle the plots of the model.
 * Here the data are read in a temporary file that is created here. It will contain:
 *    - Header with the smoothing coefficient scoef1 and scoef2
 *    - col1 the x-axis
 *    - col2 the y-axis for the data
 *    - col3 a smooth of the y-axis at a level scoef1
 *    - col4 a smooth of the y-axis at a level scoef2
 *    - col5 the y-axis for the model
 * If ps is set to 1, then the plot is written in file_out. Otherwise, it is shown
 * on the screen.
*/

    bool ps;
    std::stringstream ss;

 	ps=1;

    Gnuplot gn;
	
	// Create a temporary file with the data that we wish to plot
	write_data_tmp(data_struc->x, data_struc->y, model); // This function create a temporary file that is used by gnuplot

   	if (ps == 0) { 
		gn << "set term X11 \n";
    	} else {
		gn << "set term post eps enhanced color font 'Times-Bold, 15'\n";
		if (phase == "init") {gn << "set out '" + file_model_init_diags + ".eps'\n";}
		if (phase == "buffer") {gn << "set out '" + file_model_buffer_diags + ".eps'\n";}
		if (phase == "final") {gn << "set out '" + file_model_final_diags + ".eps'\n";}
   	 }
   	 // Setup the common plot configuration
         gn << "set datafile commentschars '#!'\n";
         gn << "set autoscale \n"; // scale axis automatically
         gn << "unset label \n";
         gn << "set key default \n"; // restore the default position for the legends
         gn << "set key top right \n"; // legend will be on bottom right

  	  // The plot for the likelihood of all parallel chains
        gn << "set xlabel '" + data_struc->xlabel + " " + data_struc->xunit + "'" << std::endl;
       	gn << "set ylabel '" + data_struc->ylabel + " " + data_struc->yunit + "'" << std::endl;
	if(show_original_data == 1){
		gn << "plot '" + file_data_tmp +  "' using 1:2 with lines title 'No smoothing'" + 
		    ",'" + file_data_tmp +  "' using 1:3 with lines title 'scoef="+dbl_to_str(data_scoef1) +"'" +
		    ",'" + file_data_tmp +  "' using 1:4 with lines title 'scoef="+dbl_to_str(data_scoef2) +"'" +
		    ",'" + file_data_tmp +  "' using 1:5 with lines title 'Model'";
	} else{
		gn << "plot '" + file_data_tmp +  "' using 1:3 with lines title 'scoef="+dbl_to_str(data_scoef1) +"'" +
		    ",'" + file_data_tmp +  "' using 1:4 with lines title 'scoef="+dbl_to_str(data_scoef2) +"'" +
		    ",'" + file_data_tmp +  "' using 1:5 with lines title 'Model'";

	}
  	gn << std::endl;
		
    //shell_exec("rm " + file_data_tmp); // erase the temporary file // Not work I dont know why

}

void Diagnostics::write_data_tmp(const VectorXd& x, const VectorXd& y, const VectorXd& z){
/* 
 * Write on a file the data, smooth(data, scoef1), smooth(data,scoef2) and the model.
 *
*/
	VectorXi Nchars(3), precision(3);
	VectorXd ys1, ys2;

	std::ofstream outfile;

	Nchars << 20, 20, 20;
	precision << 10, 10, 7;

	ys1=smooth(x, y, data_scoef1);
	ys2=smooth(x, y, data_scoef2);

	outfile.open(file_data_tmp.c_str());
	if(outfile.is_open()){
	
		for(int i=0; i<x.size(); i++){
			outfile << std::setw(Nchars[0]) << std::setprecision(precision[0]) << x[i];
			outfile << std::setw(Nchars[1]) << std::setprecision(precision[1]) << y[i];
			outfile << std::setw(Nchars[1]) << std::setprecision(precision[1]) << ys1[i];
			outfile << std::setw(Nchars[1]) << std::setprecision(precision[1]) << ys2[i];
			outfile << std::setw(Nchars[2]) << std::setprecision(precision[2]) << z[i];
			outfile << std::endl;
		}
		outfile.close();
	}  
	else {
        	file_error(file_data_tmp, "openfile", "Diagnostics::write_data_tmp");
    }
}


void Diagnostics::write_txt_matrix_dbl(const MatrixXd& Mat, const std::string file_out){
/* 
 * Write on a file a MatrixXd
 *
*/
	const int Nchars=20;
	const int precision=10;
	std::ofstream outfile;
	
	outfile.open(file_out.c_str());
	if(outfile.is_open()){
		outfile << std::setw(Nchars) << std::setprecision(precision) << Mat << std::endl;
		outfile.close();
	}  
	else {
        file_error(file_out, "openfile", "Diagnostics::write_txt_matrix_dbl");
   }

}

std::vector<double> Diagnostics::vectXd_to_vec(const VectorXd& vecXd_in){

   std::vector<double> vec;
   vec.resize(vecXd_in.size());
   VectorXd::Map(&vec[0], vecXd_in.size()) = vecXd_in;
   return vec;
}



VectorXd Diagnostics::smooth(const VectorXd& xin, const VectorXd& in, const double scoef){
/* 
 * Return a boxcar smooth of a vector. The smooth coeficient
 * is in natural unit of the vector (ie, microHz not in bins)
*/

  VectorXd out(in.size());
  double dx;
  double Nbins_d;
  int Nbins;
  dx=std::abs(xin[1]-xin[0]);
  Nbins_d=scoef/dx;
  Nbins=int(Nbins_d);

  if(Nbins > 1){
  	if(Nbins%2 != 0){Nbins=Nbins +1;} // ensure that the Nbins/2 is an integer

         for (int i=Nbins; i<in.size()-Nbins; i++){
		out[i]=in.segment(i- Nbins/2, Nbins).sum()/Nbins; // Centered average
	  }
	  for (int i=0; i<Nbins; i++){
		out[i]=in.segment(0, i+Nbins/2).sum()/(i+Nbins/2); // Right side average
	  }
	  for (int i=in.size()-1-Nbins; i<in.size(); i++){
		out[i]=in.segment(i-Nbins/2, in.size()-1-i).sum()/(in.size()-1-i); // Left side average
	  }
  } else {
	out=in;
  }

return out;
}


void Diagnostics::fwrite_histogram(const Data_Nd samples, const std::string file_out, const size_t Nclasses, const int ind_param){
/*
 * Compute and write the histogram of a given parameter of index ind_param
*/

    std::string vars;
    Data_Nd tmppdf;

    // configure the histogram 
#ifdef TAMCMC_WITH_GSL
    gsl_histogram* hist = gsl_histogram_alloc(Nclasses);

    assert( hist != NULL );
    gsl_histogram_set_ranges_uniform(hist, samples.data.col(ind_param).minCoeff(), samples.data.col(ind_param).maxCoeff());

    // add the samples to the histogram
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);

    for(int j = 0; j < samples.data.rows(); j++ ) {
      gsl_histogram_increment(hist, samples.data(j,ind_param));            
    }

    // Create the pdf for the structure 'histogram'
    gsl_histogram_pdf *MyHistPdf = gsl_histogram_pdf_alloc(Nclasses);
    assert( MyHistPdf != NULL );
    int status = gsl_histogram_pdf_init(MyHistPdf, hist);
    assert( status != GSL_EDOM );

    FILE* tmpout;
    tmpout=fopen(file_hist_tmp.c_str(), "w");
    int stat = gsl_histogram_fprintf(tmpout, hist, "%f", "%f");
    fclose(tmpout);
    tmppdf=read_txt_output_params(file_hist_tmp, Nclasses, 0); // read the temporary file

    std::ofstream fileout_stream;
    fileout_stream.open(file_out.c_str());
    if(fileout_stream.is_open()){
	fileout_stream << "! variable_name= " << samples.labels[ind_param] << std::endl;
	for(int i=0; i<tmppdf.data.rows();i++){
		fileout_stream << (tmppdf.data(i,0) + tmppdf.data(i,1))/2 << "  " << tmppdf.data(i, 2) << std::endl; // simplifies the writting by using a two column format
	}
    } else{
        file_error(file_out, "openfile", "Diagnostics::fwrite_histogram");
    }
    fileout_stream.close();
#endif
    //shell_exec("rm " + tmpfile); // erase the temporary file

    // These Deletes make crash the code! Investigate how to avoid this (small) memory leaks here...
    //delete hist;
    //delete MyHistPdf;
    //delete tmpout;
    //delete T;
    //----------------------------------------------------------------------------------

}

Data_Nd Diagnostics::read_txt_output_params(const std::string file, const int Nlines, const bool verbose){
/*
	This function is intended to read the pdfs output text files, created by fwrite_histogram
	This is necessary because it is more efficient to plot from a file with gnuplot than
	plotting variables in memory (and it is more portable).
*/
  int cpt=0;
  size_t pos=0;

  int Ncols;
  Eigen::MatrixXd params_out;
  std::string line;
  std::vector<std::string> keyword, varnames;
  std::ifstream myfile;
  Eigen::VectorXd words;

  Data_Nd data;

  myfile.open(file.c_str());
  if (myfile.is_open()){
    if(verbose == 1) {std::cout << "File: " << file << "..." << std::endl;}
    while(getline(myfile,line)){
       if((pos = line.find_first_of("!#")) == std::string::npos){
		words=str_to_Xdarr(line, " \t"); // split the line when we encounter a tab or a white space
		if(cpt ==0){
			Ncols=words.size();
			data.data.resize(Nlines, Ncols);
			data.data.row(cpt)=words;
		} else {
			data.data.row(cpt)=words;
		}	
	cpt=cpt+1;
	} else {
	  keyword=strsplit(line, "=");
	  if(keyword[0]== "! variable_names"){
		varnames=strsplit(keyword[1], " ");
		for(int j=0; j<varnames.size(); j++) {data.labels.push_back(strtrim(varnames[j]));}
	  } else {
		data.header.push_back(strtrim(line)); // any comment or extra parameters are put in the header
	  }
	}
    }
    myfile.close();
  } else {
    file_error(file, "openfile", "Diagnostics::read_txt_output_params");
  }

return data;
}


bool Diagnostics::get_model_buffer_diags(){
	return model_buffer_diags;
}

bool Diagnostics::get_model_final_diags(){
	return model_final_diags;
}


std::string Diagnostics::formated_int_to_str(const int ind_param){

	std::string out;
	if (ind_param <10) out="00" + int_to_str(ind_param);
	if ((ind_param <100) && (ind_param >= 10)) out="0" + int_to_str(ind_param);
	if ((ind_param <1000) && (ind_param >= 100)) out= int_to_str(ind_param);
	if (ind_param >= 1000){
		std::cout << "Error when formating the pdf output files" << std::endl;
		std::cout << " The maximum number of parameters that can be handled is 1000" << std::endl;
		std::cout << " If the problem is bigger than that, you might consider editing the diagnostics.cpp file" << std::endl;
		std::cout << " The program will stop now" << std::endl;
		exit(EXIT_FAILURE);
	}
	return out;
}

void Diagnostics::gnuplt_pdfs_diags(const std::string dir_pdfs_files, const std::vector<std::string> varnames){
// 
// * Simple function that handle the plots of the probability density functions.
// * Here the data are read in a file per parameter that must be created before
// * execution. 
// * Here ind_var is the index of the variable that will be plotted. 
// *
  
    int m, Ncount;
	std::string postop, posbottom, posleft, posright;
	
    Gnuplot gn;
	
	// --- Initialize the plot positionning variables ------
	Ncount=0; // Initialise the counter of plots per window
	postop="0.90";
	posbottom="0.52";
	posleft="0.06";
	posright="0.41";
	// -----------------------------------------------------
	
	gn << "set term post eps enhanced color font 'Times-Bold, 16' size 29.7cm,21cm \n";
	gn << "set out '" + file_pdfs_diags + ".eps'\n";

   	 // Setup the common plot configuration
 	 gn << "set datafile commentschars '#!'\n";
 	 gn << "set autoscale \n"; // scale axis automatically
 	 gn << "unset label \n";
 	 gn << "set key default \n"; // restore the default position for the legends
 	 gn << "set key top right \n"; // legend will be on bottom right
	 gn << "set multiplot\n";
         // On a split screen, we show the pdfs for each of the variables
         for(int ind_var=0; ind_var<varnames.size(); ind_var++){
	   	// Move to the next corner in a cycle, from top to bottom and from left to right
	   	if(Ncount == 1){ // Set the position to the second plot position
			posleft="0.47";
			posright="0.90";
		}
                if(Ncount == 2){ // Set position to the third plot position
	   		postop="0.46";
			posbottom="0.06";
			posleft="0.06";
			posright="0.41";
	   	}
                if(Ncount == 3){ // Set position to the third plot position
			posleft="0.47";
			posright="0.90";
	   	}
                if(Ncount == 4){ // Set the position to the first plot
	   		gn << "set multiplot\n";
	   		postop="0.90";
			posbottom="0.52";
			posleft="0.06";
			posright="0.41";
   			Ncount=0;
	   	}
	 	gn << "set tmargin at screen " + postop + "\n";
                gn << "set bmargin at screen " + posbottom + "\n";
	 	gn << "set lmargin at screen " + posleft + "\n";
                gn << "set rmargin at screen " + posright + "\n";
 	  	gn << "set xlabel '" + varnames[ind_var] + "'" << std::endl;
 	  	gn << "set ylabel 'PDF (Counts)'" << std::endl;
		gn << "plot '" + dir_pdfs_files + output_root_name + "pdf_" + formated_int_to_str(ind_var) + ".txt' using 1:2" + " with boxes lc rgb'green' notitle\n" ;
	   	   
	   	 Ncount=Ncount+1;
          }
	
        gn << "unset multiplot" << std::endl; // come back to single zone
        std::cout << "     The output files with information on the PDF has been updated" << std::endl;
        std::cout << "         -  Histograms are stored in: " << dir_pdfs_files << std::endl;
        std::cout << "         -  Plots of the histogram are in : " << file_pdfs_diags << ".eps" << std::endl;
        
}

void Diagnostics::gnuplt_pdfs_diags_main(const int i){
// 
// * The main function that handle the plots of the probability density functions.
// * Here the data are read in a file per parameter that must be created before
// * execution. 
// * Here ind_var is the index of the variable that will be plotted. 
// *

   bool verbose=0;  
   std::string file_out;
  
   std::string dir_pdfs_files=output_dir + "pdfs/";
   Data_Nd data_samples;
   Eigen::MatrixXd data_array;
   Params_hdr hdr;

   if(pdfs_diags == 1 && ((i%Nbuffer == 0 && i !=0) || i == Nsamples-1)){
   	// Load the samples
   	if(file_out_format == "text"){
		data_samples=read_txt_output_params(filename_params, Nsamples, verbose);
	} 
	else{
		hdr=read_params_header(filename_params_hdr); // Get the metadata from the header file 
		data_array=read_bin_matrix_dbl(filename_params, hdr.Nvars, Nsamples, "dbl");
		data_samples.data=data_array;
		data_samples.labels=hdr.variable_names;
	}
    //Compute and write the histogram on a output file... Do it only if GSL is available
    #ifdef TAMCMC_WITH_GSL
        for(int ind_param=0; ind_param<data_samples.data.cols(); ind_param++){
            file_out=dir_pdfs_files + output_root_name + "pdf_" + formated_int_to_str(ind_param) + ".txt";
            fwrite_histogram(data_samples, file_out, Nclasses, ind_param);
        }
        gnuplt_pdfs_diags(dir_pdfs_files, data_samples.labels);
    #else
        std::cout << "         - Program compiled without GSL support ==> Skip histogram rendering" << std::endl;
    #endif
    }
}

// -----------------------------------------------------------------------------------
// ------------------------ READING FUNCTIONS FOR THE OUTPUTS ------------------------
// -----------------------------------------------------------------------------------
Eigen::MatrixXd Diagnostics::read_bin_matrix_dbl(const std::string binfile, const long Ncols, const long Nrows, const std::string type){
/*
 * Function that read the outputs file that contains inputs in double and in Matricial format (Columns - Rows).
 * It requires as input:
 * 	- The name of the file with the data
 * 	- The number of columns in the matrix 
 *	- The number of rows in the matrix
 *      - The output format (string): type="dbl" (double), type="ldbl" (long double)
*/
	long Nread_rows;
	double val_dbl=0;
	long double val_ldbl=0;
	size_t size_ldbl=sizeof(val_ldbl);
	size_t size_dbl=sizeof(val_dbl);

	Eigen::MatrixXd vals(Nrows, Ncols);
	std::ifstream file;

	file.open(binfile.c_str(), std::ios::binary); // Open the binary file in read only
	if(file.is_open()){
	   Nread_rows=0;
	   if(type == "ldbl"){
		  while(!file.eof() && (Nread_rows<Nrows)){
			for(int i=0; i<Ncols; i++){
				file.read(reinterpret_cast<char*>(&vals(Nread_rows,i)), size_ldbl);
			}
			Nread_rows=Nread_rows+1;
		}
	   }
	   if(type == "dbl"){
		  while(!file.eof() && (Nread_rows<Nrows)){
			for(int i=0; i<Ncols; i++){
				file.read(reinterpret_cast<char*>(&vals(Nread_rows,i)), size_dbl);
			}
			Nread_rows=Nread_rows+1;
		}
	   }
		file.close();
	} else {
        file_error(binfile, "openfile", "Diagnostics::read_bin_matrix_dbl");
	}

	vals.conservativeResize(Nread_rows-1, Ncols);
	return vals;
}


Ptempering_out Diagnostics::read_bin_parallel_temp_params(const std::string binfile, const long Nrows){
/*
 * Function that read the outputs file that contains inputs as writen in the file that contains
 * the information about the parallel tempering. Return a structure of type Buffer_parallel_tempering
 * It requires as input:
 * 	- The name of the file with the data
 *	- The number of rows in the matrix
*/

	char data[1]; // used to read the booleans
	bool val_bool=0;
	int val_int=0;
	double val_dbl=0;
	size_t size_int=sizeof(val_int);
	size_t size_dbl=sizeof(val_dbl);
	std::ifstream file;
	Ptempering_out ptemp_read;

	ptemp_read.attempt_mixing.resize(Nrows); 
	ptemp_read.chain0s.resize(Nrows); 
	ptemp_read.Pswitchs.resize(Nrows); 
	ptemp_read.switcheds.resize(Nrows); 

	file.open(binfile.c_str(), std::ios::binary); // Open the binary file in read only
	 if (file.is_open()){
		for (int i=0; i<Nrows; i++){ // we read the buffer
			file.read(&data[0], 1); // Read first the attempt_mixing variable (booleans)
			ptemp_read.attempt_mixing[i]=(data[0] >> val_bool);
			file.read(reinterpret_cast<char*>(&val_int), size_int); // then the chain0s variable (integers)
			ptemp_read.chain0s[i]= val_int;
			file.read(reinterpret_cast<char*>(&val_dbl), size_dbl); // them the Pswitchs (doubles)
			ptemp_read.Pswitchs[i]=val_dbl;
			file.read(&data[0], 1); // finaly read the switcheds (booleans)
			ptemp_read.switcheds[i]=(data[0] >> val_bool);
			std::cout << ptemp_read.attempt_mixing[i] << "  " << ptemp_read.chain0s[i] << "  " << ptemp_read.Pswitchs[i] << "  " << ptemp_read.switcheds[i] << std::endl;			
			
		}
		file.close();
  	}
  	else {
        file_error(binfile, "openfile", "Diagnostics::read_bin_parallel_temp_params");
    }
return ptemp_read;
}


Moves_out Diagnostics::read_bin_moves(const std::string binfile, const long Nchains, const long Nrows){
/*
 * Function that read the outputs file that contains inputs as writen in the file that contains
 * the information about the moves. Return a structure of type Buffer_parallel_tempering
 * It requires as input:
 * 	- The name of the file with the data
 *	- The number of rows in the matrix
*/
	char data[1];
	bool val_bool=0;
	double val_dbl=0;
	size_t size_bool=sizeof(val_bool);
	size_t size_dbl=sizeof(val_dbl);
	std::ifstream file;
	Moves_out outmoves;

	outmoves.Pmoves.resize(Nrows, Nchains); // A matrix of double
	outmoves.moveds.resize( Nrows , std::vector<bool>( Nchains , val_bool ) ); // A matrix of booleans

	file.open(binfile.c_str(), std::ios::binary); // Open the binary file in read only
	if (file.is_open()){
		for (int i=0; i<Nrows; i++){ // we read the buffer
			for(int j=0; j<Nchains; j++){ // The probabilities for each moves
				file.read(reinterpret_cast<char*>(&outmoves.Pmoves(i,j)), size_dbl);
			}
			for(int j=0; j<Nchains; j++){ // Booleans telling whether we moved or not
				file.read(&data[0], 1); // Read first the attempt_mixing variable (booleans)
				outmoves.moveds[i][j]=(data[0] >> val_bool);
			}
		}
		file.close();
  	}
  	else {
        file_error(binfile, "openfile", "Diagnostics::read_bin_moves");
    }
return outmoves;
}


Params_hdr Diagnostics::read_params_header(const std::string file){

  int cpt=0;
  size_t pos=0;

  long varcase;
  std::string line;
  std::vector<std::string> keyword, varnames;
  std::ifstream myfile;

  Params_hdr hdr;

  myfile.open(file.c_str());
  if (myfile.is_open()){
    while(getline(myfile,line)){
	keyword=strsplit(line, "=");
	varcase=9; // Defaut case is header case
	if(keyword[0]== "! Nchains"){
		varcase=0;
	}
	if(keyword[0]== "! Nvars"){
		varcase=1;
	}
	if(keyword[0]== "! Ncons"){
		varcase=2;
	}
	if(keyword[0]== "! relax"){
		varcase=3;
	}
	if(keyword[0]== "! plength"){
		varcase=4;
	}
	if(keyword[0]== "! constant_names"){
		varcase=5;
	}
	if(keyword[0]== "! constant_values"){
		varcase=6;
	}
	if(keyword[0]== "! variable_names"){
		varcase=7;
	}
	if(keyword[0]== "! Nsamples"){
		varcase=8;
	}
	//for(int el=0; el<keyword.size();el++){std::cout << keyword[el] << std::endl;}
	//std::cout << " -----" << std::endl;
	switch(varcase){
		case 0: 
		  hdr.Nchains=0;
		  varnames=strsplit(keyword[1], " ");
		  std::stringstream(strtrim(varnames[0])) >> hdr.Nchains;
		  break;
		case 1: 
		  hdr.Nvars=0;
		  varnames=strsplit(keyword[1], " ");
		  std::stringstream(strtrim(varnames[0])) >> hdr.Nvars;
		  break;
		case 2: 
		  hdr.Ncons=0;
		  varnames=strsplit(keyword[1], " ");
		  std::stringstream(strtrim(varnames[0])) >> hdr.Ncons;
		  break;
		case 3: 
		  varnames=strsplit(keyword[1], " ");
		  hdr.relax.resize(hdr.Ncons + hdr.Nvars);
		  hdr.relax.setConstant(-1);
		  for(int j=0; j<varnames.size(); j++) {
			std::stringstream(strtrim(varnames[j])) >> hdr.relax[j];
		  }
		  break;
		case 4: 
		  varnames=strsplit(keyword[1], " ");
		  hdr.plength.resize(varnames.size());	
		  hdr.plength.setConstant(-1);
		  for(int j=0; j<varnames.size(); j++) {
			std::stringstream(strtrim(varnames[j])) >> hdr.plength[j];
		  }
		  break;
		case 5: 
		  varnames=strsplit(keyword[1], " ");
		  for(int j=0; j<varnames.size(); j++) {
			hdr.constant_names.push_back(strtrim(varnames[j]));
		  }
		  break;
		case 6: 
		  varnames=strsplit(keyword[1], " ");
		  hdr.constant_values.resize(varnames.size());	
		  for(int j=0; j<varnames.size(); j++) {
			std::stringstream(strtrim(varnames[j])) >> hdr.constant_values[j];
		  }
		  break;
		case 7: 
		  varnames=strsplit(keyword[1], " ");
		  for(int j=0; j<varnames.size(); j++) {
			hdr.variable_names.push_back(strtrim(varnames[j]));
		  }
		  break;
		case 8: 
		  hdr.Nsamples=0;
		  varnames=strsplit(keyword[1], " ");
		  std::stringstream(strtrim(varnames[0])) >> hdr.Nsamples;
		  break;
		default:
		  hdr.header.push_back(strtrim(line)); // any comment or extra parameters are put in the header
	}
    }
    myfile.close();
  } else {
    file_error(file, "openfile", "Diagnostics::read_params_header");
  }

return hdr;
}


Evidence_out Diagnostics::evidence_calc(const VectorXd& Tcoefs, const MatrixXd& Likelihoods){
/*
 * Function that use the temperatures and the likelihoods in order to compute de evidence of the model.
 * Returns all values in a structure of type Evidence_out
*/
	int Npts;
	Evidence_out results;
	
	Npts=evidence_interpolation_factor*Tcoefs.size();

	results.interpolation_factor=evidence_interpolation_factor;

	results.beta.resize(Tcoefs.size());
	results.L_beta.resize(Tcoefs.size());
	
	results.beta_interp.resize(Npts);
	results.L_beta_interp.resize(Npts);

	// Compute beta and L_beta == P<P(D|M,I)>_beta
	for(int i=0; i<Tcoefs.size(); i++){
		results.beta[i]=1./Tcoefs[i];
		results.L_beta[i]=Likelihoods.col(i).sum()/Likelihoods.rows();
	}

	// Interpolate beta and P<P(D|M,I)>_beta
	results.beta_interp=quad_interpol( results.beta, Npts );
	results.L_beta_interp=quad_interpol( results.L_beta, Npts );

	// Compute the Integral, which is the average of all the interpolated values
    results.evidence=results.L_beta_interp.sum()/Npts;

   return results;
}

void Diagnostics::write_evidence(const Evidence_out& evidence, const int i){

     std::ostringstream strg;
     std::ofstream outfile_evidence;
     std::string file;

     const int Nchars=20;
     const int precision=10;

    file= file_evidence_diags + ".txt";
	/////// Write the evidence ////////

	 if (firstpass == 1) {
		outfile_evidence.open(file.c_str()); // Overwrite any existing file in PLAIN ASCII
	 } else {
		outfile_evidence.open(file.c_str(), std::ofstream::app); // std::app is for append
	 }
	 if (outfile_evidence.is_open()){
		if (firstpass == 1){ // Write the header only if it is the first time that we write OR if an append of an existing file was requested and that the file exists
    			outfile_evidence << "# This is an output file for the evidence. Evidence is calculated after a quadratic interpolation of L_beta. \n";
			outfile_evidence << "# This file contains values for the L_beta[0:Nchains-1], the evidence calculated at each time the buffer was written \n" ;
			outfile_evidence << "# col(1): Number of samples used to compute the evidence \n";
			outfile_evidence << "# col(2:2+Nchains): averaged probability <P(D|M,I)> over the samples of each chain \n";
			outfile_evidence << "# col(2+Nchains+1): Evidence P(M|D, I) computed by (1) interpolation and (2) averaging \n";
			outfile_evidence << "! beta=" << evidence.beta.transpose() << "\n";
			outfile_evidence << "! interpolation_factor=" << evidence.interpolation_factor << "\n";
		}
		strg.str(std::string());
		strg << i; // The current index of the samples (which is then also the total number of samples so far)
		strg << " "; 
		strg << std::setw(Nchars) << std::setprecision(precision) << evidence.L_beta.transpose(); // The L_beta table (not interpolated)
		strg << " "; 
		strg << std::setw(Nchars) << std::setprecision(precision) << evidence.evidence << "\n"; // The evidence value

		outfile_evidence  << strg.str().c_str();
		outfile_evidence.flush(); // Explicitly specify to flush the data into the disk
		strg.str(std::string()); // clear the strg buffer
		
		outfile_evidence.close();
  	}
  	else {
        file_error(file, "openfile", "Diagnostics::write_evidence");
	}

}

std::string Diagnostics::shell_exec(const std::string cmd){ // For some reason, this seems to not work...
/*
 * Small program that execute a given shell command and return the result of that command
 * Taken from stackoverflow.com
*/

	char buffer[128];
	std::string result="";

	FILE* pipe = popen(cmd.c_str(), "r");
	if(!pipe){ throw std::runtime_error("popen() failed!");}
	try{
		while(!feof(pipe)){
			if(fgets(buffer, 128, pipe) != NULL){
				result+= buffer;
			}
		}
	} catch (...){
		pclose(pipe);
		throw;
	}
	pclose(pipe);
	return result;
}

int Diagnostics::file_error(const std::string file, const std::string error_type, const std::string fct_name){
/*
 * Function that show an error message depending on error_type
 * It receives the filename, fonction name and returns an error code
 * A negative error code correspond to a fatal error. A positive
 * code correspond to a warning.
 * At the moment, the return value is not used: Automatic exit after the error code
*/

    if(error_type == "openfile"){
        std::cout << "Fatal Error in " << fct_name << std::endl;
        std::cout << "Could not open the file: " << file << std::endl;
        std::cout << "Check that the destination path exists" << std::endl;
        std::cout << "The progran will exit now" << std::endl;
    }
    if(error_type == ""){
        std::cout << "Unknown Fatal Error in " << fct_name << std::endl;
        std::cout << "Could not open the file: " << file << std::endl;
        std::cout << "Neeed debuging..." << std::endl;
        std::cout << "The progran will exit now" << std::endl;
    }

    exit(EXIT_FAILURE);
    return -1;
}

