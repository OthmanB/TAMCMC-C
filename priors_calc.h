/*
 * priors_apply.h
 *
 *  Contains the headers of the function that calculate the penalization
 *  according to presets. So far, the only priors preset is for
 *  a main sequence star: priors_MS_Global()
 * 
 *  Created on: 02 Mar 2016
 *      Author: obenomar
 */
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::VectorXi;

double priors_MS_Global(VectorXd params, VectorXi param_length, VectorXd priors_params);
double priors_Test_Gaussian(VectorXd params, VectorXi param_length, VectorXd priors_params);
