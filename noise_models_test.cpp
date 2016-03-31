/*
 * noise_models.cpp
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#include <iostream>
#include <iomanip>

#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;

VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);
int main();

int main(){

	const long Nx=21, Nharvey=2;
	VectorXd x(Nx), y(Nx), noise_params(Nharvey*3+1);
	
	noise_params(0)=10.;
	noise_params(1)=2.;
	noise_params(2)=0.7;
	noise_params(3)=0.5;
	noise_params(4)=0.7;
	noise_params(5)=1.4;
	noise_params(6)=0.05;
	
	x.setLinSpaced(Nx, 0, 10);
	y.setZero();

	y=harvey_like(noise_params, x, y, Nharvey);

	std::cout << "x     y" << std::endl;
	for(int i=0; i<Nx; i++){
		std::cout << x(i) << "  " << y(i) << std::endl;
	}

return 0;
}
