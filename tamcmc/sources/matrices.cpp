/*
 * matrix3d.cpp
 *
 *  Created on: 16 Mar 2016
 *      Author: obenomar
 */
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using Eigen::MatrixXd;
using Eigen::VectorXd;


////////////// 3D Matrix ////////////////


MatrixXd** initialize_3dMatrix(const int depth, const int Nrows, const int Ncols){
/* 
 * Initialize to 0 an array of matrix 
*/
    //std::cout << "depth = " << depth << "   Nrows = " << Nrows << "   Ncols = " << Ncols << std::endl;
    MatrixXd** matrix = new MatrixXd* [depth];
    for (int d1 = 0; d1 < depth; d1++) {
            matrix[d1] = new (MatrixXd)(MatrixXd::Zero(Nrows,Ncols)); //(MatrixXd::Zero(Nrows,Ncols).array());
    }
    return matrix;
} 

MatrixXd** initset_3dMatrix(const MatrixXd& m_init, const int depth){
/* 
 * Initialize an array of matrix using an initial matrix m_init
*/
    MatrixXd** matrix = new MatrixXd* [depth];
    for (int d1 = 0; d1 < depth; d1++) {
            matrix[d1] = new (MatrixXd)(m_init); //(m_init.array());
    } 
    return matrix;
} 

MatrixXd** copy_3dMatrix(MatrixXd** m3d_in, const int depth){
/* 
 * Hard copy an array of matrix from m3d_in into m3d_out 
*/
    MatrixXd** m3d_out = new MatrixXd* [depth];
    for (int d1 = 0; d1 < depth; d1++) {
	    //std::cout << "d1= " << d1 << std::endl;
	    //std::cout << "*m3d_in[d1]= " << *m3d_in[d1] << std::endl;
            m3d_out[d1] = new (MatrixXd)(*m3d_in[d1]);
    } 
    return m3d_out;
} 

void set_3dMatrix(MatrixXd** matrix3d, const MatrixXd& m_in, const int position){
/* 
 * Replace the matrix which is within an array 
 * of matrix at the position p. 
*/
    //std::cout << "position = " << position << std::endl;
    delete matrix3d[position];
    matrix3d[position] =new (MatrixXd)(m_in); //(m_in.array()); 
} 

////////////// 4D Matrix ////////////////


MatrixXd*** initialize_4dMatrix(const int depth, const int Nchains, const int Nrows, const int Ncols, const long double value=0.0) {
    MatrixXd*** matrix = new MatrixXd** [depth];
    for (int d1 = 0; d1 < depth; d1++) {
        matrix[d1] = new MatrixXd* [Nchains];
        for (int d2 = 0; d2 < Nchains; d2++) 
            matrix[d1][d2] = new (MatrixXd)(MatrixXd::Ones(Nrows,Ncols).array() * value);
    }
    return matrix;
} 


MatrixXd*** initset_4dMatrix(const MatrixXd& m_init, const int depth1, const int depth2){
/* 
 * Initialize an array of matrix using an initial matrix m_init
*/
    MatrixXd*** matrix = new MatrixXd** [depth1];
    for (int d1 = 0; d1 < depth1; d1++) {
        matrix[d1] = new MatrixXd* [depth2];
        for (int d2 = 0; d2 < depth2; d2++) 
            matrix[d1][d2] = new (MatrixXd)(m_init);
    }
    return matrix;
} 

void set_4dMatrix(MatrixXd*** matrix4d, const MatrixXd& m_in, const int position1, const int position2){
/* 
 * Replace the matrix which is within an array 
 * of matrix at the position p1 and p2. 
*/
    //std::cout << "position = " << position << std::endl;
    delete matrix4d[position1][position2];
    matrix4d[position1][position2] = new (MatrixXd)(m_in); //(m_in.array()); 
} 

// DOES NOT WORK YET
//MatrixXd*** copy_4dMatrix(MatrixXd** m4d_in, int depth1, int depth2){
///* 
// * Hard copy an array of matrix from m4d_in into m4d_out.
// * Here, we work by 2d plans.
//*/
//    MatrixXd*** m4d_out = new MatrixXd** [depth1];
//   for (int d1 = 0; d1 < depth1; d1++) {
//	m4d_out[d1] = new MatrixXd* [depth2];
//	for (int d2 = 0; d2 < depth2; d2++) {
//            m4d_out[d1][d2] = new (MatrixXd)(*m4d_in[d1][d2]);
//	}
//    } 
//    return m4d_out;
//} 

