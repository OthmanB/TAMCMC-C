/*
 * pmode.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;


struct Pmode{
		int el; // degree of the mode
		int en; // index or radial order of the mode
		double nu; // central frequency
		double width; // FWHM
		double height; // Height
		double a1; // average splitting
		double a2; // effect on splitting of the centrifugal distorsion
		double a3; // effect on splitting of latitudinal rotation
		double inc; // stellar inclination
		int relax_nu;
		int relax_width;
		int relax_height;
		int relax_a1;
		int relax_a2;
		int relax_a3;
		int relax_inc;
		char common_params[]; // lists name of parameters in common between modes (e.g. {'inc', 'a1', 'a2', 'a3'})

};





