/*
 *  Created on: 02 Mar 2016
 */
#pragma once
//#include "noise_models.h"
#include "build_lorentzian.h"
//#include "interpol.h"
//#include "function_rot.h"

using Eigen::VectorXi;
using Eigen::VectorXd;

double eta0_fct(const VectorXd& fl0_all);

VectorXd model_MS_Global_a1nl_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1n_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1l_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); // Added on 10 Feb 2020
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); // Added on 10 Feb 2020
VectorXd model_MS_Global_a1acta3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1etaa3_Harvey1985(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1acta3_Harvey1985(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);

VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);

VectorXd model_MS_local_basic(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_MS_local_Hnlm(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); // Added on 14 Feb 2020 : Handles local fit without fitting inclination... instead used Hnlm

VectorXd model_Test_Gaussian(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);
VectorXd model_Harvey_Gaussian(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);

VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); 
VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); 
VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); 
VectorXd model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false);

VectorXd model_MS_Global_a1n_a2a3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); // Added on 18 Jan 2021: Handles the a2 coefficient with n free
VectorXd model_MS_Global_a1nl_a2a3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); // Added on 18 Jan 2021: Handles the a2 coefficient with n,l free (but only n dependence accounted for)
VectorXd model_MS_Global_a1a2a3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); // Added on 18 Jan 2021: Handles the a2 coefficient with 2 terms : a constant term + 1 slope term in nu

VectorXd model_MS_Global_aj_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams=false); // Added on 18 Jan 2021: Handles the a2 coefficient with 2 terms : a constant term + 1 slope term in nu

VectorXd model_MS_Global_a1etaAlma3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams); // Added on 31 Mar 2021

void write_star_params(const VectorXd& spec_params, const VectorXd& raw_params, const VectorXi& plength, const MatrixXd& mode_params, const MatrixXd& noise_params, const std::string file_out, const std::string modelname, const std::string name_params);

bool debug(const VectorXd& model, const long double Hl, const long double fl, const long double a1, const long double eta, const long double a3,
               const long double asym, const long double Wl, const long double el, const long double step, const double inclination, const VectorXd& ratios,
               const long double trunc_c, bool exit_c=true);
bool debug_solver(const VectorXd& x, const VectorXd& fl1_all, const VectorXd& fl0_all, const int el, const long double delta0l, 
                  const long double DPl, const long double alpha_g, const long double q_star, const long double sigma_p_l1);
