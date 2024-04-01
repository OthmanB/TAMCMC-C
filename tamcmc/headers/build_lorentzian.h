#pragma once
#include <math.h>
#include <Eigen/Dense>
#include "function_rot.h"
#include "../../external/Alm/Alm_cpp/data.h"
# include <iostream>
# include <iomanip>
#include <string>

using Eigen::VectorXd;
using Eigen::VectorXi;

struct Optim_L{
    // To store the local zone where the computation was made using optimum_lorentzian_calc_* functions
    VectorXd y;
    int i0; // first element where the calculated block should go
    int N; // Number of data points
};

double Qlm(const int l, const int m);
VectorXi set_imin_imax(const VectorXd& x, const int l, const double fc_l, const double gamma_l, const double f_s, const double c, const double step);

VectorXd build_l_mode_a1l_etaa3(const VectorXd& x_l, const double H_l,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V);
VectorXd build_l_mode_a1etaa3(const VectorXd& x_l,  const double H_l,  const double fc_l,  const double f_s,  const double eta, const double a3,  const double asym,  const double gamma_l, const int l, const VectorXd& V);
VectorXd build_l_mode_a1etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l);
VectorXd build_l_mode_a1l_etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym, double gamma_l, const int l);
VectorXd build_l_mode_a1l_a2a3(const VectorXd& x_l, const double H_l,  const double fc_l,  const double f_s1,  const double f_s2,  const double a2,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V);
VectorXd build_l_mode_a1a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V);

VectorXd optimum_lorentzian_calc_a1l_etaa3(const VectorXd& x, const VectorXd& y, const  double H_l, const  double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);
VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd& x, const VectorXd& y,  const double H_l,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const VectorXd& V,  const double step, const double c);

VectorXd optimum_lorentzian_calc_a1l_etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s1,  const double f_s2,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l, double step, const double c);
VectorXd optimum_lorentzian_calc_a1etaa3_v2(const VectorXd& x,  const VectorXd& y,  const VectorXd& H_lm,  const double fc_l,  const double f_s,  const double eta,  const double a3,  const double asym,  const double gamma_l, const int l,  const double step, const double c);

VectorXd optimum_lorentzian_calc_a1l_a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);

// Generalized form of a-coefficients
VectorXd build_l_mode_aj(const VectorXd& x_l, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta, const double asym, const double gamma_l, const int l, const VectorXd& V);
//VectorXd optimum_lorentzian_calc_aj(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, 
//        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
//        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);
Optim_L optimum_lorentzian_calc_aj(const VectorXd& x, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c);


VectorXd build_l_mode_ajAlm(const VectorXd& x_l, const double H_l, const double fc_l, const double a1, const double a3, const double a5, 
    const double eta0, const double epsilon_nl, const VectorXd& thetas, 
    const double asym, const double gamma_l, const int l, const VectorXd& V, 
    const std::string filter_type, const gsl_funcs);
Optim_L optimum_lorentzian_calc_ajAlm(const VectorXd& x,  const double H_l,  const double fc_l,  const double a1, const double a3, const double a5, 
		const double eta0, const double epsilon_nl, const VectorXd& thetas,  
                const double asym,  const double gamma_l, const int l,  const VectorXd& V,  
                const double step, const double c, const std::string filter_type, const gsl_funcs);

