/*
 * matrix3d.h
 *
 *  Created on: 16 Mar 2016
 *      Author: obenomar
 */
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd** initialize_3dMatrix(const int depth, const int Nrows, const int Ncols);
MatrixXd** initset_3dMatrix(MatrixXd m_init, int depth);
void set_3dMatrix(MatrixXd** matrix3d, MatrixXd m_in, int position);

MatrixXd*** initialize_4dMatrix(int depth, int Nchains, int Nrows, int Ncols, long double value=0.0);
MatrixXd*** initset_4dMatrix(MatrixXd m_init, int depth1, int depth2);
void set_4dMatrix(MatrixXd*** matrix4d, MatrixXd m_in, int position1, int position2);
