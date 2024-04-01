/*
 *  Created on: 02 Mar 2016
 */
#pragma once
#include "build_lorentzian.h"
#include <vector>
#include "../../external/Alm/Alm_cpp/data.h"
#include "data.h"

using Eigen::VectorXi;
using Eigen::VectorXd;

VectorXd model_RGB_asympt_aj_AppWidth_HarveyLike_v4_parallel(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_RGB_asympt_aj_CteWidth_HarveyLike_v4_parallel(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
