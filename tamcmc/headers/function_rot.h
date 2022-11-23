/*
 * function_rot.h
 *
 *  Created on: 11 Feb 2016
 *      Author: obenomar
 */
#pragma once
#include <math.h>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

double combi(const int n, const int r);
double dmm(const int l, const int m1, const int m2, const double beta);
int factorial(const int n);
MatrixXd function_rot( const int l, const double beta);
VectorXd amplitude_ratio(const int l, const double beta);

