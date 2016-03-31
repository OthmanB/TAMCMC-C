/*
 * likelihoods_test.cpp
 *
 * This file contains the functions that generate the likelihood
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#include <Eigen/Dense>
# include <iostream>
# include <iomanip>
#include "likelihoods.h"

using Eigen::VectorXd;

int main(){

	long p, Nx=10;
	long double f;
	VectorXd y(Nx), model(Nx), sigma(Nx);
	
	y.setConstant(1);
	sigma.setConstant(0.1);
	model.setConstant(1.1);
	p=1;
	
	std::cout << "Testing the likelihood_chi22p function..." << std::endl;
	std::cout << "  degree of freedoms p = " << 2*p << std::endl;
	std::cout << "  y  / model(y) " << std::endl;	
	for(int i=0; i< model.size(); i++){std::cout << "  " << y(i) << "    " << model(i) << std::endl;}
	f=likelihood_chi22p(y, model, p);
	
	std::cout << "likelihood_chi22p(y, model, p) = " << f << std::endl;
	std::cout << "----------------------------------------" << std::endl;

	std::cout << "Testing the likelihood_chi_square function..." << std::endl;
	std::cout << "  y  / sigma / model(y) " << std::endl;	
	for(int i=0; i< model.size(); i++){std::cout << "  " << y(i) << "    " << sigma(i) << "   " << model(i) << std::endl;}
	f=likelihood_chi_square(y, model, sigma);
	std::cout << "likelihood_chi_square(y, model, sigma) = " << f << std::endl;
	std::cout << "Expected result: -10" << std::endl;
	std::cout << "----------------------------------------" << std::endl;


}
