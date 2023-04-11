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

Input_Data build_init_ajfit(const aj_files i_ajfit, const double a1_obs);
//Data set_observables_ajfit(aj_files i_ajfit, Data data);
Data_Nd set_observables_ajfit(const aj_files i_ajfit, const Data_Nd data, const int x_col, const int y_col, const int ysig_col);
aj_files read_ajfile(const std::string cfg_model_file);