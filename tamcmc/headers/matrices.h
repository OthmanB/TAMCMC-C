/*
 * matrix3d.h
 *
 *  Created on: 16 Mar 2016
 *      Author: obenomar
 */
#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd** initialize_3dMatrix(const int depth, const int Nrows, const int Ncols);
MatrixXd** initset_3dMatrix(const MatrixXd& m_init, const int depth);
void set_3dMatrix(MatrixXd** matrix3d, const MatrixXd& m_in, const int position);

MatrixXd*** initialize_4dMatrix(const int depth, const int Nchains, const int Nrows, const int Ncols, const long double value=0.0);
MatrixXd*** initset_4dMatrix(const MatrixXd& m_init, const int depth1, const int depth2);
void set_4dMatrix(MatrixXd*** matrix4d, const MatrixXd& m_in, const int position1, const int position2);
