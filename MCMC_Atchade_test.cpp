/*
 *  MCMC_Atchade.cpp
 *
 *  This file contains the core program for the MCMC
 *  as it is defined in Atchade Y.F., 2006, Meth. Comp. In Applied Probab. 8, 235
 *  This MCMC is an adaptive truncated-Langevin Metropolis Hasting algorithm.
 *  The strategy involves a Robins-Monroe stochastic algorithm during the learning
 *  process and optionaly can use the local hyperspace gradient to enhance the sampling
 *  (MALA option). 
 *  The adpative part of the algorithm tries to minimize differences between successive 
 *  samples. The minimized quantities are:
 *      - The first moment (mean). It should be stable and converge towards mu when
 *        Nsamples ---> infinity
 *      - The second moment (covarmat) should be stable and converge towards 
 *        The covariance matrix near the global maxima
 *      - A scaled variance (mean variance). It is here to ease the search
 *        for the best covariance matrix. It should be stable when N ---> infinity
 *  WARNING: The current implementation does not use the Gradient.
 *  Remark: The algorithm proposed by Atchade does not implement parallel tempering
 *          Here, we added the parallel tempering which greatly improve the performance.
 *
 *  Created on: 02 Mar 2016
 *      Author: obenomar
 */

#include <math.h>
//#include <time.h>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using Eigen::VectorXd;
using Eigen::MatrixXd;





/*void update_proposal(var_cur, acceptance, m_ind){
	

}
*/



/*// --------------------------------------------------
int main(){
	long double likelihood=-10.;
	long double lambda=1.7;
	long j=8, m=12;
	long Nx=4;
	VectorXd x1(Nx), x2(Nx), tmp(Nx);
	MatrixXd M(Nx,Nx);

	srand(time(NULL));

	std::cout << "Testing apply_temperature()..." << std::endl;
	std::cout << "  likelihood=" << likelihood << std::endl;
	std::cout << "  lambda=" << lambda << std::endl;
	std::cout << "  j=" << j << std::endl;
	std::cout << "  ===> T=lambda^(j-1)" <<  std::endl;
	std::cout << "  ===> likelihood/T=" << apply_temperature(likelihood, lambda, j) << std::endl;
	
	
	std::cout << "Testing the random discrete number generator random_int_vals..." << std::endl;
	std::cout << "   number of parallel chain m = " << m << std::endl;
  	std::cout << "   Generating 20 random discrete number in [1,m]... " << std::endl;
	for(int i=0; i<20; i++){
		std::cout <<  "      " << random_int_vals(m) << std::endl;;
	}

	std::cout << "Testing the scalar product <x1|M|x2> using eigen..." << std::endl;
	x1.setLinSpaced(Nx, 0, 3);
	x2.setLinSpaced(Nx, 0, 3);
	x2=x2/4 + tmp.setConstant(1);
	M.setIdentity();
	M(0,1)=0.5;
	M(0,2)=1.5;
	M(1,0)=0.15;
	M(1,3)=0.25;
	M(2,0)=0.5;
	M(2,1)=0.1;
	M(3,2)=0.5;
	std::cout << " Vector x1 =" << x1.transpose() << std::endl;
	std::cout << " Vector x2 =" << x2.transpose() << std::endl;
	std::cout << "    Matrix M " << std::endl;
	std::cout << M << std::endl;
	std::cout << "Calculated <x1|M|x2> = " << x1.transpose()*M*x2 << std::endl;
	std::cout << "Expected   <x1|M|x2> = 13.587500" << std::endl;

	std::cout << "Testing the product |x1><x2| using eigen..." << std::endl;
	std::cout << x1*x2.transpose() << std::endl;
return 0;
}
*/
