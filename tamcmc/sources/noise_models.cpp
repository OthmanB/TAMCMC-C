/*
 * noise_models.cpp
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>

using Eigen::VectorXd;

VectorXd harvey_like(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey){
	/* This function calculate a sum of harvey like profile + a white noise and adds 
	   these profiles to a initial input vector y. The function assumes that the 
	   inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx), y_out(Nx); //, y_out(Nx);
	
	white_noise.setConstant(noise_params.tail(1)(0));
	y_out=y;
	int i;
	for(i=0; i<Nharvey;i++){
		if(noise_params(cpt+1) != 0){
		  		tmp=((1e-3)*noise_params(cpt+1)*x).array().pow(noise_params(cpt+2)); 
				tmp=noise_params(cpt)*(tmp + ones.setConstant(1)).cwiseInverse();
				y_out= y_out + tmp;
		} 	
		cpt=cpt+3;
	}
	y_out=y_out + white_noise;

return y_out;
}

VectorXd harvey1985(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey){
	/* This function calculate a sum of harvey profile + a white noise and adds 
	 * these profiles to a initial input vector y. The Harvey profile differ from the
	 * Harvey-like by the fact that H0_1985= H0_like * tc is correlated to the timescale tc. 
	 * There is also a 2pi factor in the denominator. The function assumes that the 
	 * inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

	const long double pi = 3.141592653589793238L;
	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx), y_out(Nx);
	
	white_noise.setConstant(noise_params.tail(1)(0));
 	y_out=y;

	for(long i=0; i<Nharvey;i++){
		if(noise_params(cpt+1) != 0){
			tmp=((1e-3)*2*pi*noise_params(cpt+1)*x).array().pow(noise_params(cpt+2)); // Denominator
			tmp=noise_params(cpt)*noise_params(cpt+1)*(tmp + ones.setConstant(1)).cwiseInverse(); // Numerator/Denominator
			y_out= y_out + tmp; // Numerator/Denominator + white noise
		}
		cpt=cpt+3;
	}
	y_out=y_out + white_noise;

return y_out;
}

double get_ksinorm(const double b, const double c, const Eigen::VectorXd& x) {
    double integral = 0.0;
    double h = x(1) - x(0); // assuming x is equally spaced
    
    for (int i = 0; i < x.size(); i++) {
        double term = 1.0 / (1.0 + std::pow(x(i) / b, c));
        
        if (i == 0 || i == x.size() - 1) {
            integral += 0.5 * term;
        } else {
            integral += term;
        }
    }
    
    integral *= h;
    double ksi = b/integral;
    
    return ksi;
}


VectorXd eta_squared_Kallinger2014(const VectorXd& x){
	const double x_nyquist=x.maxCoeff();
	VectorXd eta(x.size());
	eta=sin(0.5*M_PI*x.array()/x_nyquist)/(0.5*M_PI*x.array()/x_nyquist);
	if (x[0] == 0){ // This to avoid the Division by 0
		eta[0]=1;
	}
	return eta.array().square();
}


VectorXd Kallinger2014(const double numax, const double mu_numax, const VectorXd& noise_params,const VectorXd& x, const VectorXd& y){
/*
	Using notations from Table 2 of Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
	Not that here we assume the instrumental noise to be Pinstrument(nu) = 0
	The noise_params must have parameters in this order:
	    - Granulation noise: 4 parameters to create a0 and b0
		- Noise a1, a2 : ka, sa, t
		- Noise b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
		- Noise b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
	Such that at the end we have: [ka,sa,t,k1,s1,c1, k2,s2,c2, N0]
*/
	const long Nx=x.size();
	VectorXd ones(Nx), white_noise(Nx), tmp0(Nx), tmp1(Nx), tmp2(Nx), y0(Nx), Power(Nx);
	ones.setOnes();
	// Compute the Leakage effect as a sinc function (Eq 1 of Kallinger+2014)
	const VectorXd eta_squared=eta_squared_Kallinger2014(x);
	// Compute b1, b2 and a
	const double a0=std::abs(noise_params[0]*std::pow(std::abs(numax),noise_params[1])); // Very Low frequencies
	const double b0=std::abs(noise_params[2]*std::pow(std::abs(numax + mu_numax),noise_params[3]));
	const double c0=std::abs(noise_params[4]);
	const double a1=noise_params[5];
	const double a2=noise_params[6];
	const double b1=std::abs(noise_params[7]*std::pow(std::abs(numax + mu_numax),noise_params[8])); // Intermediate
	const double b2=std::abs(noise_params[10]*std::pow(std::abs(numax + mu_numax),noise_params[11])); // Below numax
	const double c1=std::abs(noise_params[9]);
	const double c2=std::abs(noise_params[12]);
	const double N0=std::abs(noise_params[13]);
	// Compute the normalisation constants ksi1 and ksi2
	const double ksi0=get_ksinorm(b0, c0, x);
	const double ksi1=get_ksinorm(b1, c1, x);
	const double ksi2=get_ksinorm(b2, c2, x);
	y0=y;
	// White noise first, added to the input y-vector
	white_noise.setConstant(N0);
 	Power=y + white_noise;
	// Granulation SuperLorentzian
	tmp0=(x/b0).array().pow(c0); // Denominator
	tmp1=(eta_squared * ksi0 * std::pow(a0,2)/b0).cwiseProduct((tmp0 + ones).cwiseInverse()); // Numerator/Denominator
	Power=Power + tmp1;
	// First SuperLorentzian
	tmp0=(x/b1).array().pow(c1); // Denominator
	tmp1=(eta_squared * ksi1 * std::pow(a1,2)/b1).cwiseProduct((tmp0 + ones).cwiseInverse()); // Numerator/Denominator
	Power=Power + tmp1;
	// Second SuperLorentzian
	tmp0=(x/b2).array().pow(c2); // Denominator
	tmp2= (eta_squared *ksi2*std::pow(a2,2)/b2).cwiseProduct((tmp0 + ones).cwiseInverse()); // Numerator/Denominator
	Power=Power + tmp2;
	return Power;

}