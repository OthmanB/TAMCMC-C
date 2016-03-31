/*
 * function_rot.cpp
 *
 *  Created on: 11 Feb 2016
 *      Author: obenomar
 */
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main ( int argc, char *argv[] );
double combi(int n, int r);
double dmm(int l, int m1, int m2, double beta);
int factorial(int n);
MatrixXd function_rot( int l, double beta);
VectorXd amplitude_ratio(int l, double beta);

int main ( int argc, char *argv[] ){
    int l=2;
    double beta=30, sum;
    MatrixXd mat;
    VectorXd heights;

    cout << "\n";
    cout << "Test factorial...\n";
    for(int i=0; i<5; i++){
	cout << "factorial(" << i << ")=" << factorial(i) << "\n";
    }
    cout << "\n";
    cout << "Testing dmm calling combi and factorial... \n";
    
    for(int m1=-l; m1<l; m1++){
	for(int m2=-l; m2<l; m2++){
		cout << "l=" << l << " m1=" << m1 << " m2=" << m2 << "\n";
       		sum=dmm(l, m1, m2, beta);
		cout << "sum=" << sum << "\n";       		
	}
    }
    cout << "\n";
    cout << "Testing function_rot calling dmm... \n";
    mat=function_rot(l, beta);
    std::cout << mat << std::endl;       		

    cout << "\n";
    cout << "Testing amplitude_ratio which run everything... \n";
    
    for(l=0; l<=2; l++){ 
	for(beta=0; beta<=90; beta=beta+20){
    	heights=amplitude_ratio(l, beta);
    	cout << "l=" << l << ", i=" << beta << ":" <<endl;
        cout << heights << endl;
	}
    }
}

