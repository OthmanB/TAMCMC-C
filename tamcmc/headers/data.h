/*
 * data.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "gnuplot-iostream.h"
#include "../../external/Alm/Alm_cpp/data.h" // To load interpolation grids

//using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::VectorXd;

// The general structure that contain the input data (spectrum, lightcurve, ...) that has to be analysed
struct Data{
		VectorXd x; // In the case of a 1D fit (e.g. fit of Mass, [Fe/H], ...), this variable will be ignored by the model function.
		VectorXd xrange; // Contains the min and max of x... in practice, it is used to limit the data range if requested by the cfg file (e.g. freq_range option in the .MCMC)
		VectorXd y;
		VectorXd sigma_y;
		long Nx; // Ny is not checked but should be as long as x
		std::string xlabel; // label for x-axis. In ASCII file, the axis labels are identified by the symbol ! at the begining of the line
		std::string ylabel;
		std::string xunit;
		std::string yunit;
		std::vector<std::string> header; // Any kind of information that is useful. e.g. if it is a test, model info + model inputs. Any header (In ASCII, marked by #), will be put there.
};


// Used to keep in memory all of the information from that tabulated priors provided by the user (if any)
struct tabpriors{
	bool empty = true; // defined if it was used or not
    int depth; // number prior tables
    VectorXd Nrows; 
    VectorXd Ncols;
    MatrixXd** data_3d; // Contain submatrices. Each of them was read from a file and give the prior table (either 1D or Nd)
	std::vector<std::vector<std::string>> headers;
	std::vector<std::vector<std::string>> labels;
};


struct Input_Data{
	std::string model_fullname; // The fullname of the model that is going to be processed
	std::vector<std::string> inputs_names;
	std::vector<std::string> priors_names;
	VectorXi priors_names_switch; // Used instead of priors_names in loops. This allows the use of the switch(x) - case statements instead of if/else.
	VectorXd inputs;
	VectorXi relax;
	MatrixXd priors;
	VectorXi plength;
	VectorXd extra_priors; // Contains extra parameters that could be used for priors
	tabpriors tabulated_priors; // Contains all table of priors given by the users in the *.priors files
};

// A Generic structure that helps to encapsulate a Matrix of information along with some metadata
struct Data_Nd{
	MatrixXd data; // 2D array in which each columns are expected to contain the values of a given parameter
	std::vector<std::string> header;
	std::vector<std::string> labels;
	std::vector<std::string> units;
};

// Header of the parameters
// Useful to read the stand alone ASCII header written when dealing with BINARY OUTPUTS.
struct Params_hdr{
	std::vector<std::string> header; // Any comment
	int Nsamples; // Total number of samples (no necessarily the actual number of written samples)
	int Nchains;
 	int Nvars;
	int Ncons;
	//std::vector<bool> relax;
	VectorXi relax;
	VectorXi plength;
	std::vector<std::string> constant_names; // The name of the constants
	VectorXd constant_values;
	std::vector<std::string> variable_names; // The name of the variables
};

struct MCMC_files{
	std::string ID;
	double Dnu;
	double numax;
	double err_numax;
	double C_l;
	VectorXi els;
	VectorXd freq_range;
	std::vector<std::string> param_type;
	//std::vector<double> freqs_ref;conservativeResize
	VectorXd freqs_ref;
	std::vector<bool> relax_freq, relax_gamma, relax_H;

	//std::vector<double> hyper_priors;
	MatrixXd hyper_priors; // Col 0 is the fref_bias, Col 1-5 Are for the prior parameters
	std::vector<std::string> hyper_priors_names; // This gets the prior name (U, GU, G, etc..). Check io_ms_global to see how this is read
	MatrixXd eigen_params;
	VectorXd noise_params;
	MatrixXd noise_s2;

	std::vector<std::string> common_names;
	std::vector<std::string> common_names_priors;
	MatrixXd modes_common;
	//std::string filter_type; // For Alm models only
};

struct aj_files{
	std::string filter_type;
	double Dnu;
	VectorXd els;
	VectorXd nu_nl;
	bool do_a2;
	bool do_a4;
	bool do_a6;
	bool do_CFonly;
};


// Structure that keep information of the derivatives
struct Deriv_out{
	VectorXd xderiv;
	VectorXd deriv;
	VectorXd error;
};

struct gnuplt_Data {
/*
 * This is an encapsulator for data when ploting with gnuplot-iostream.h
*/
    double x;  // x axis value
    double y1;             // y axis series 1
    double y2;             // y axis series 2
    double y3;             // y axis series 3
};

typedef std::vector<gnuplt_Data> gnuplt_Dataset;

namespace gnuplotio {
    template<>
    struct TextSender<gnuplt_Data> {
        static void send(std::ostream &stream, const gnuplt_Data &v) {
           // TextSender<std::string>::send(stream, v.x);
            stream << " ";
            TextSender<double>::send(stream, v.x);
            stream << " ";
            TextSender<double>::send(stream, v.y1);
            stream << " ";
            TextSender<double>::send(stream, v.y2);
            stream << " ";
            TextSender<float>::send(stream, v.y3);
        }
    };
}

struct Data_Basic{
	std::vector<std::string> strarr; // Any comment
	VectorXi vecXi; // Case number
};

struct external_data{ // A structure designed as a container for any kind of additional static data that can be called by models
	gsl_funcs Alm_interp_gate; // Table for Alm interpolations for band of activity
	gsl_funcs Alm_interp_gauss; // Table for Alm interpolations for gauss activity zone
	gsl_funcs Alm_interp_triangle; // Table for Alm interpolations for triangular activity zone
};
