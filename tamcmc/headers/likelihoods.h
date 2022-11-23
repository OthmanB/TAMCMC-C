/*
 * likelihoods.h
 *
 * This file contains the header for functions that generate the likelihood
 *
 *  Created on: 28 Feb 2016
 *      Author: obenomar
 */
#pragma once
#include <Eigen/Dense>

using Eigen::VectorXd;

long double likelihood_chi22p(const VectorXd& y, const VectorXd& model, const long p);
long double likelihood_chi_square(const VectorXd& y, const VectorXd& model, const VectorXd& sigma);
