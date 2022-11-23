/*
Date: 16 Nov 2021
Testing of the acoefs.cpp function present in [root]/tamcmc/sources and its header 
Adapted from acoefs.py from the acoefs_check project (see github)
This code verifies creates an executable that allow to verify that when
you provide a central frequency nu_nl, along with aj coefficients j={1,...,6}
(1) the decomposition in acoeficient is properly done and along the way
(2) Tnlm is properly calculated
(3) Snlm is properly calculated
*/
#include <math.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "../../tamcmc/headers/acoefs.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;

void usage(int argc, char* argv[]);
int  options(int argc, char* argv[]);
void test_acoefs(const int l, const long double nu_nl, const long double a1, const long double a2, const long double a3, const long double a4, const long double a5, const long double a6);

void test_acoefs(const int l, const long double nu_nl, const long double a1, const long double a2, const long double a3, const long double a4, const long double a5, const long double a6){
	/*
		With this function you can see whether by putting some input
		a-coefficients, you do retrieve them after you measure directly
		the symetrical splittings (Tnlm) and asymetrical splitting (Snlm)
	*/
	VectorXd nu_nlm, Tn, Sn, anl;

	nu_nlm=nunlm_from_acoefs(nu_nl, l, a1, a2, a3,a4,a5,a6);
	std::cout << " nu_nlm:" << std::endl; 
	std::cout << "    " << nu_nlm.transpose() << std::endl;

	std::cout << " Tnlm : " << std::endl;
	Tn=Tnlm(nu_nlm, l);
	for (int m=0; m<Tn.size();m++){
		std::cout << "    Tn" << l << m+1 << " = " << Tn[m] << std::endl;	
	}
	std::cout << " Snlm : " << std::endl;
	Sn=Snlm(nu_nlm, l);
	for (int m=0; m<Sn.size();m++){
		std::cout << "    Sn"<< l << m+1 << " = " << Sn[m] << std::endl;
	}
	std::cout << " Input a-coefficients:" << std::endl;
	std::cout << "    a1 = " << a1 << std::endl;
	std::cout << "    a2 = " << a2 << std::endl;
	std::cout << "    a3 = " << a3 << std::endl;
	std::cout << "    a4 = " << a4 << std::endl;
	std::cout << "    a5 = " << a5 << std::endl;
	std::cout << "    a6 = " << a6 << std::endl;
	
	std::cout << " Retrieved a-coefficients from the function eval_acoefs():" << std::endl;
	anl=eval_acoefs(l, nu_nlm);
	std::cout << "    a1 = " << anl[0] << std::endl;
	std::cout << "    a2 = " << anl[1] << std::endl;
	std::cout << "    a3 = " << anl[2] << std::endl;
	std::cout << "    a4 = " << anl[3] << std::endl;
	std::cout << "    a5 = " << anl[4] << std::endl;
	std::cout << "    a6 = " << anl[5] << std::endl;
}

void usage(int argc, char* argv[]){

			std::cout << "Unrecognized arguments" << std::endl;
			std::cout << "     - To execute: " << argv[0] << " <mode_degree> <nu_nl> <a1> <a2> <a3> <a4> <a5> <a6>  "<< std::endl;
			exit(EXIT_FAILURE);
}

int options(int argc, char* argv[]){
	// Check that we have the correct number of arguments
	int val=-1;	
	if(argc == 9){
		val=1; 
	} 
	if (val == -1){ // Error code
		usage(argc, argv);
	} 
	if (val > 0 ){
		return val; // Execution code val
	} else{
		return -1; // Default value is to return an error code
	}
}

int main(int argc, char* argv[]){
	int msg_code;
	int l;
	long double nu_nl, a1,a2,a3,a4,a5,a6;

	msg_code=options(argc, argv);
	if(msg_code == -1){
		std::cout << "Error detected in options. Cannot proceed. Debug required." << std::endl;
		std::cout << "The program will exit now" << std::endl;
		exit(EXIT_FAILURE);
	} else{
		std::istringstream(argv[1]) >> l;
		std::istringstream(argv[2]) >> nu_nl;	
		std::istringstream(argv[3]) >> a1;
		std::istringstream(argv[4]) >> a2;
		std::istringstream(argv[5]) >> a3;
		std::istringstream(argv[6]) >> a4;
		std::istringstream(argv[7]) >> a5;
		std::istringstream(argv[8]) >> a6;
		//test_acoefs(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8])
		test_acoefs(l, nu_nl, a1, a2, a3, a4, a5, a6);
	}
}

