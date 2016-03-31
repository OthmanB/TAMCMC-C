/*
 * noise_models.cpp
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;

VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);

VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey){
	/* This function calculate a sum of harvey like profile + a white noise and adds 
	   these profiles to a initial input vector y. The function assumes that the 
	   inputs are in the following order: [H0, tc0, p0, ..., Hn, tcn, pn, N0]
	*/

	const long Nx=x.size();
	int cpt=0;
	VectorXd ones(Nx), white_noise(Nx), tmp(Nx); //, y_out(Nx);
	
	white_noise.setConstant(noise_params.tail(1)(0));
	//y_out.setZero();
 
	for(long i=0; i<Nharvey;i++){
		tmp=((1e-3)*noise_params(cpt+1)*x).array().pow(noise_params(cpt+2)); 
		tmp=noise_params(cpt)*(tmp + ones.setConstant(1)).cwiseInverse();
		y= y + tmp;
		cpt=cpt+3;
	}
	y=y + white_noise;

return y;
}

