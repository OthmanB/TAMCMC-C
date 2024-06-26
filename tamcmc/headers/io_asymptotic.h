#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "data.h" // contains the structure Data
#include "io_models.h"
#include "function_rot.h"
#include "io_ms_global.h"
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;


VectorXd settings_aj_splittings_RGB(const int i, const MCMC_files inputs_MS_global, Input_Data* Snlm_in); 
MCMC_files read_MCMC_file_asymptotic(const std::string cfg_model_file, const bool verbose);
Input_Data build_init_asymptotic(const MCMC_files inputs_MS_global, const bool verbose, const double resol);
