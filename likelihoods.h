/*
 * likelihoods.h
 *
 * This file contains the header for functions that generate the likelihood
 *
 *  Created on: 28 Feb 2016
 *      Author: obenomar
 */
#include <Eigen/Dense>

using Eigen::VectorXd;

long double likelihood_chi22p(VectorXd y, VectorXd model, long p);
long double likelihood_chi_square(VectorXd y, VectorXd model, VectorXd sigma);
