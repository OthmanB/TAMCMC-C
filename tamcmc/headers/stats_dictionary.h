/*
 * priors_dictionary.h
 *
 *  Contains the headers for the function that define the priors.
 *  Currently set priors: 
 *       - Uniform
 *       - Gaussian
 *       - Jeffrey
 *       - Uniform-Gaussian
 *       - Gaussian-Uniform
 *       - Gaussian-Uniform-Gaussian
 *
 *  Created on: 09 Apr 2016
 *      Author: obenomar
 */
#pragma once
#include <iostream>
#include <iomanip>
#include <math.h>
#include <Eigen/Dense>
 
 long double logP_gaussian(const long double mean, const long double sigma, const long double x);
 long double logP_gaussian_truncated(const long double mean, const long double sigma, const long double x, const long double bmin, const long double bmax);
 long double logP_uniform(const long double b_min, const long double b_max, const long double x);
 long double logP_uniform_abs(const long double b_min, const long double b_max, const long double x);
 long double logP_uniform_cos(const long double b_min, const long double b_max, const long double x);
 long double logP_multivariate_gaussian( const Eigen::VectorXd& mean, const Eigen::MatrixXd& Matrix, const Eigen::VectorXd& x);
 long double logP_jeffrey(const long double hmin, const long double hmax, const long double h);
long double logP_jeffrey_abs(const long double hmin, const long double hmax, const long double h);
 long double logP_uniform_gaussian(const  long double b_min, const long double b_max, 
			const long double sigma, const long double x);
 long double logP_gaussian_uniform( const long double b_min, const  long double b_max, 
			const long double sigma, const long double x);
 long double logP_gaussian_uniform_gaussian(const  long double b_min, const long double b_max, 
			const long double sigma1, const long double sigma2,const  long double x);
