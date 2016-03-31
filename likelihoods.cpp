/*
 * likelihoods.cpp
 *
 * This file contains the functions that generate the likelihood
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;

long double likelihood_chi22p(VectorXd y, VectorXd model, long p){

	long double f;
//  exact formula :
//	b=(p-1)*alog(p)-alog(gamma(p)+(p-1)*total(alog(y))-p*total(alog(model_final))-p*total(y/model_final)
//  constant terms are dropped as they cancel out in the logarithmic space
	f=(y.transpose()*model.cwiseInverse()).sum() + (model.array().log()).sum();
	f=-p*f;
return f;
}


long double likelihood_chi_square(VectorXd y, VectorXd model, VectorXd sigma){

	long double f;
// Usual chi square definition
	f=-((y - model).array().square() * sigma.array().square().cwiseInverse()).sum();

return f;
}

