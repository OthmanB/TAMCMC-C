#include <Eigen/Dense>
#include <cmath>
//#include "interpol.h"

using Eigen::VectorXd;

double lin_interpol(VectorXd x, VectorXd y, double x_int){
	/* Very simple function that linearly interpolate at the position x_int
	   a curve (x,y). Here we assume that x.size() = y.size()
	*/

	long i=0, Nx=x.size();
	double a=0,b=0; 
	VectorXd xtmp, ytmp;
	
	// ---- case of an actual interpolation -----
	//std::cout << "x_int=" << x_int << std::endl;

	if(x_int >= x.head(1)(0) && x_int <= x.tail(1)(0)){ 
		while((x_int < x[i] || x_int > x[i+1]) && i<Nx-1){
			//std::cout << "x_int < x[i] || x_int > x[i+1]) && i<Nx-1 =" << ((x_int < x[i] || x_int > x[i+1]) && i<Nx-1) << std::endl;
			i=i+1;
		}
		if(i==0){i=i+1;} // case where we never passed by the loop because x_int < x[0] || x_int > x[1] = True	
		a=(y[i] - y[i-1])/(x[i] - x[i-1]); // slope
		b=y[i-1] - a*x[i-1]; // ordinate at origin
	}
	// ---- case of an extrapolation toward the lower edge ----
	if(x_int  < x.head(1)(0)){
		//std::cout << "Warning: we have to perform an extrapolation toward the lower edge" << std::endl;
		a=(y[1] - y[0])/(x[1] - x[0]); // slope
		b=y[0] - a*x[0]; // ordinate at origin 
	}
	// ---- case of an extrapolation toward the upper edge ----
	if(x_int  > x.tail(1)(0)){
		//std::cout << "Warning: we have to perform an extrapolation toward the upper edge" << std::endl;
		ytmp=y.tail(2);
		xtmp=x.tail(2);
		a=(ytmp[1] - ytmp[0])/(xtmp[1] - xtmp[0]); // slope
		b=ytmp[0] - a*xtmp[0]; // ordinate at origin 
	}
return a*x_int+b;
}
