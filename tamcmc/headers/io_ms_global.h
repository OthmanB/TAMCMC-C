#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

#include "data.h" // contains the structure Data
#include "string_handler.h"
#include "io_models.h" // Contains the function that create the final vector

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

MCMC_files read_MCMC_file_MS_Global(const std::string cfg_model_file, const bool verbose);
Input_Data build_init_MS_Global(const MCMC_files inputs_MS_global, const bool verbose, const double resol);
short int set_noise_params(Input_Data *Noise_in, const MatrixXd& noise_s2, const VectorXd& noise_params);
short int fatalerror_msg_io_MS_Global(const std::string varname, const std::string param_type, const std::string syntax_vals, const std::string example_vals);
double getnumax(const VectorXd& fl, const VectorXd& Hl);
Input_Data set_width_App2016_params_v1(const double numax, Input_Data width_in);
Input_Data set_width_App2016_params_v2(const double numax, Input_Data width_in);
VectorXd settings_aj_splittings(const int i, const MCMC_files inputs_MS_global, Input_Data* Snlm_in, const int aj_switch);