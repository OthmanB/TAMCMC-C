/*
 * build_lorentzian_test.cpp
 *
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

//using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd build_l_mode(VectorXd x_l, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V);
VectorXd optimum_lorentzian_calc(const VectorXd x, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V, double step);
int main();

int main(){

	long Nxl;
	int l=1;
	double fc_l, gamma_l, f_s, H_l, a2,a3, step;
	VectorXd x_l(Nxl), V(2*l+1), x_l2(Nxl);
	VectorXd r, r2;

	std::ofstream file_out;

	clock_t begin_time = clock();


	//std::cout << "Values for x_l:" << std::endl;
	//std::cout << x_l << std::endl;
	//std::cout << "---------------" << std::endl;

	Nxl=(4000. - 500.)*(1460e0*86400e-6); //4 years observation
	x_l.setLinSpaced(Nxl, 500, 4000);
	std::cout << "  - Spectra generated using 4 years observations ==> Number of data points = " << Nxl << std::endl;
	std::cout << "  - Spectra generated in the range 500 - 4000 microHz = " << std::endl;
	std::cout << "  - step = " << step << std::endl;
	step=x_l(1)-x_l(0);

	l=1;
	H_l=1;
	fc_l=1500.;
	f_s=0.5;
	a2=1e-5;
	a3=0.01;
	gamma_l=0.75;
	V(0)=0.25;
	V(1)=0.5;
	V(2)=0.25;
	
	std::cout << "Calculating model without approximation.." << std::endl;
	r=build_l_mode(x_l, H_l, fc_l, f_s, a2, a3, gamma_l, l, V);
	
	//std::cout << "Values for Lorentzian:" << std::endl;
	//std::cout << "x_l   r"  << std::endl;
	//for(int i=0; i<Nxl; i++){
	//	std::cout << x_l(i) << "     " << r(i) << std::endl;
	//}
	std::cout << "    Calculation finished in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

	std::cout << "Calculating model with approximation.." << std::endl;
	begin_time = clock();	
	r2=optimum_lorentzian_calc(x_l, H_l, fc_l, f_s, a2, a3, gamma_l, l, V, step);
	std::cout << "    Calculation finished in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

	std::cout << "Writting on a file results so that they can be directly compared using plots..." << std::endl;

	file_out.open("build_lorentzian_outputs.txt");
	file_out << "x    full model    trunc model " << std::endl;
	for(long i=0; i<Nxl; i++){
		file_out << x_l(i) << "    " << r(i) << "    " << r2(i) << std::endl;
	}
	file_out.close();
	std::cout << " Success. File closed." << std::endl;
return 0;
}
