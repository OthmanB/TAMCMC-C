/*
 * build_lorentzian.cpp
 *
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
//# include <iostream>
//# include <iomanip>

using Eigen::VectorXd;

VectorXd build_l_mode(VectorXd x_l, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V);
VectorXd optimum_lorentzian_calc(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V, double step);

VectorXd build_l_mode(const VectorXd x_l, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V){

	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), result(Nxl);
	double Qlm, clm;

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for a2
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
			profile=(x_l - tmp.setConstant(fc_l*(1. + a2*Qlm) + m*f_s + clm*a3)).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
	}

return result;
}

VectorXd optimum_lorentzian_calc(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V, double step){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
*/
	const double c=16.;
	double pmin, pmax;
	long imin, imax;
	VectorXd m0, x_l;
	//VectorXd mall(x.size());

	if(gamma_l >= 1 && f_s >= 1){
		if(l != 0){
			pmin=fc_l - c*(l*f_s + gamma_l);
			pmax=fc_l + c*(l*f_s + gamma_l);
		} else{
			pmin=fc_l - c*gamma_l*2.2;
			pmax=fc_l + c*gamma_l*2.2;
		}
	}
	if(gamma_l <= 1 && f_s >= 1){
		if(l !=0){
			pmin=fc_l -c*(l*f_s + 1);
			pmax=fc_l + c*(l*f_s + 1);
		} else{
			pmin=fc_l - c*2.2;
			pmax=fc_l + c*2.2;
		}
	}
	if(gamma_l >= 1 && f_s <= 1){
		if(l != 0){
			pmin=fc_l - c*(l + gamma_l);
			pmax=fc_l + c*(l + gamma_l);
		} else{
			pmin=fc_l - c*2.2*gamma_l;
			pmax=fc_l + c*2.2*gamma_l;
		}
	}
	if(gamma_l <= 1 && f_s <= 1){
		if(l !=0){
			pmin=fc_l - c*(l+1);
			pmax=fc_l + c*(l+1);
		} else{
			pmin=fc_l -c*2.2;
			pmax=fc_l +c*2.2;
		}
	}
	// ---------- Handling boundaries ----------
	//     case when the proposed values lead to (pmax - step) < min(x)....
	if( (pmax - step) < x.head(1)(0)){ 
		pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
	}
	//     case when the proposed values lead to (pmin + step) >= max(x)...
	if( (pmin + step) >= x.tail(1)(0)){
		pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
	}

	imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
	imax=ceil((pmax-x.head(1)(0))/step);

	if(imin < 0){imin=0;} // ensure that we always stay within the vector x
	if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x

	
	//std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
	//std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
	
	//std::cout << "step=" << step << std::endl;
	//std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
	//std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

	x_l=x.segment(imin, imax-imin);
 
	m0=build_l_mode(x_l, H_l, fc_l, f_s, a2, a3, gamma_l, l, V);
	//mall.setZero();
	//mall.segment(imin, imax-imin)=m0;
    y.segment(imin, imax-imin)= y.segment(imin, imax-imin) + m0;
return y;
}



