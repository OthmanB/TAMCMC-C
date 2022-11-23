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

MCMC_files read_MCMC_file_local(const std::string cfg_model_file, const int slice_ind,  const bool verbose);
Input_Data build_init_local(const MCMC_files inputs_MS_global, const bool verbose, const double resol);
short int set_noise_params_local(Input_Data *Noise_in, const MatrixXd& noise_s2, const VectorXd& noise_params, const VectorXd& freq_range);
short int fatalerror_msg_io_local(const std::string varname, const std::string param_type, const std::string syntax_vals, const std::string example_vals);
