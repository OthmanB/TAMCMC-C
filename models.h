/*
 *  Created on: 02 Mar 2016
 */

using Eigen::VectorXi;
using Eigen::VectorXd;

//#include "noise_models.h"
//#include "build_lorentizan.h"
//#include "interpol.h"
//#include "function_rot.h"

//double lin_interpol(VectorXd x, VectorXd y, double x_int);
//VectorXd amplitude_ratio(int l, double beta);
//VectorXd optimum_lorentzian_calc(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V, double step);
//VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);
VectorXd model_MS_Global(VectorXd params, VectorXi params_length, VectorXd x);
VectorXd model_Test_Gaussian(VectorXd params, VectorXi params_length, VectorXd x);

