/*
 * priors_calc.h
 *
 *  Contains the headers of the function that calculate the penalization
 *  according to presets. So far, the only priors preset is for
 *  a main sequence star: priors_MS_Global()
 * 
 *  Created on: 02 Mar 2016
 *      Author: obenomar
 */
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "stats_dictionary.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

long double priors_local(const VectorXd& params, const VectorXi& param_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch, const VectorXd& extra_priors);
long double priors_MS_Global(const VectorXd& params, const VectorXi& param_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch, const VectorXd& extra_priors);
long double priors_asymptotic(const VectorXd& params, const VectorXi& params_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch, const VectorXd& extra_priors);
long double priors_Test_Gaussian(const VectorXd& params, const VectorXi& param_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch);
long double apply_generic_priors(const VectorXd& params, const MatrixXd& priors_params, const VectorXi& priors_names_switch);
long double priors_Harvey_Gaussian(const VectorXd& params, const VectorXi& param_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch);
long double priors_ajfit(const VectorXd& params, const VectorXi& param_length, const MatrixXd& priors_params, const VectorXi& priors_names_switch);