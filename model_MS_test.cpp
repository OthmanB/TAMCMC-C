# include <iostream>
# include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>

using Eigen::VectorXd;

VectorXd model_MS_Global(const VectorXd params, const VectorXd params_length, const VectorXd x);
double lin_interpol(VectorXd x, VectorXd y, double x_int);

int main(){

	const long Nx=10;
	VectorXd x(Nx), y(Nx);
	double x_int, y_int;
	
	long Nmax=10, Nx2;
	double Dnu, a1, a2, a3, inc;
	int l, lmax=2;
	VectorXd fl0(Nmax), fl1(Nmax), fl2(Nmax), Hl0(Nmax), Wl0(Nmax), Vl(lmax), tmp(Nmax);
	VectorXd noise_params(7), params(Nmax + Nmax*(lmax+1) + lmax + 3 + Nmax + 7 + 1), params_length(10);
	VectorXd x_m, m;
	std::ofstream file_out;

	std::cout << "-------------" << std::endl;
	std::cout << "    Testing the linear interpolation..." << std::endl;
	std::cout << "-------------" << std::endl;


	//x.setLinSpaced(Nx, 0, 9);
	for(int i=0; i<Nx; i++){
		x(i)=i*1.6577;
	}
	y = 2*x.array().square();

	std::cout << "-------------" << std::endl;
	std::cout << "x, y" << std::endl;
	for(int i=0; i<Nx; i++){
		std::cout << x[i]  << " , " << y[i] << std::endl;
	}
	std::cout << "-------------" << std::endl;

	x_int=10;
	y_int=lin_interpol(x, y, x_int);

	std::cout << "x_int=" << x_int << "   " << "y_int=" << y_int << std::endl;

	std::cout << "-------------" << std::endl;
	std::cout << "   Testing the model_MS_Global function..." << std::endl;
	std::cout << "-------------" << std::endl;
	
	Nx2=(4000. - 500.)*(1460e0*86400e-6); //4 years observation
	x_m.setLinSpaced(Nx2, 500, 4000);
	std::cout << "  - Spectra generated using 4 years observations ==> Number of data points = " << Nx2 << std::endl;
	std::cout << "  - Spectra generated in the range 500 - 4000 microHz = " << Nx2 << std::endl;
	Dnu=85.;
	l=0;
	fl0.setLinSpaced(Nmax, 10, 19);
	fl0=fl0*Dnu - tmp.setConstant(-l*Dnu/2 + l*(l+1)*2.*Dnu/100);
	l=1;
	fl1.setLinSpaced(Nmax, 10, 19);
	fl1=fl1*Dnu - tmp.setConstant(-l*Dnu/2 + l*(l+1)*2.*Dnu/100);
	l=2;
	fl2.setLinSpaced(Nmax, 10, 19);
	fl2=fl2*Dnu - tmp.setConstant(-l*Dnu/2 + l*(l+1)*2.*Dnu/100);
	//l=3;
	//fl3.setLinSpaced(Nmax, 10, 19)*Dnu - l*(l+1)*1.5*Dnu/100;
	Hl0.setConstant(1.);
	Vl(0)=1.5;
	Vl(1)=0.5;
	Wl0.setConstant(1.);
	a1=1.5;
    a2=3e-3;
    a3=0.1; //0.05;
	inc=45;
	noise_params.setConstant(1.);

	params.head(Nmax)=Hl0;
	params.segment(Nmax, Vl.size())=Vl;
	params.segment(Nmax + lmax, fl0.size())=fl0;
	params.segment(Nmax + lmax + fl0.size(), fl1.size())=fl1;
	params.segment(Nmax + lmax + fl0.size() + fl1.size(), fl2.size())=fl2;
	//params.segment(Nmax+lmax+fl0.size()+fl1.size()+fl2.size(), fl3.size())=fl3;
	params(Nmax + lmax + fl0.size() + fl1.size() + fl2.size())=a1;
	params(Nmax + lmax + fl0.size() + fl1.size() + fl2.size() + 1)=a2;
	params(Nmax + lmax + fl0.size() + fl1.size() + fl2.size() + 2)=a3;
	params.segment(Nmax + lmax + fl0.size() + fl1.size() + fl2.size() + 3, Wl0.size())=Wl0;
	params.segment(Nmax + lmax + fl0.size() + fl1.size() + fl2.size() + 3 + Wl0.size(), noise_params.size())=noise_params;
	params.tail(1)(0)=inc;

	params_length(0)=Nmax;
	params_length(1)=Vl.size();
	params_length(2)=fl0.size();
	params_length(3)=fl1.size();
	params_length(4)=fl2.size();
	params_length(5)=0; // no l3 in this test
	params_length(6)=3; // splitting
	params_length(7)=Wl0.size();
	params_length(8)=noise_params.size();
	params_length(9)=1;

	for(long i=0; i<params.size();i++){
		std::cout << params(i) << std::endl;
	}
	clock_t begin_time = clock();

	m=model_MS_Global(params, params_length, x_m);
	
	std::cout << "    Calculation finished in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

	std::cout << "Writting on a file results so that they can be directly compared using plots..." << std::endl;

	file_out.open("model_MS_Global_outputs.txt");
	file_out << "x    trunc c++ model " << std::endl;
	for(long i=0; i<Nx2; i++){
		file_out << x_m(i) << "    " << m(i) << std::endl;
	}
	file_out.close();
	std::cout << " Success. File closed." << std::endl;


return 0;
}




