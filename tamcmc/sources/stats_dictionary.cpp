/*
 * priors_dictionary.cpp
 *
 *  Contains the function that define the priors
 *  Currently set priors: 
 *       - Uniform
 *       - Gaussian
 *       - Multivariate-Gaussian
 *       - Jeffrey
 *       - Uniform-Gaussian
 *       - Gaussian-Uniform
 *       - Gaussian-Uniform-Gaussian
 *
 *  Created on: 09 Apr 2016
 *      Author: obenomar
 */

#include <iostream>
#include <iomanip>
#include <math.h>
#include <Eigen/Dense>

# define PIl          3.141592653589793238462643383279502884L /* pi */

//using Eigen::MatrixXd;
//using Eigen::VectorXd;
//using Eigen::VectorXi;

// Note: These priors are handled using the 'case' statement. CASE 0 IS THE USER-DEFINED CASE (MEANS NO PRIOR USED)

long double logP_uniform(const long double b_min, const long double b_max, const long double x){
/* 
 * Calculates the log probability for a uniform prior, given the value of x
 * CORRESPONDS TO CASE 1 IN THE MAIN PROGRAM
 *
*/

long double logP;

//if ((x <= b_max) && (x >= b_min)){ // Change on 17 March 2020
if ((x <= b_max) && (x >= b_min)){
	logP=-log(std::abs(b_max-b_min));
} else{
	logP=-INFINITY;
}
return logP;
}

long double logP_uniform_abs(const long double b_min, const long double b_max, const long double x){
/* 
 * Calculates the log probability for a uniform prior, given the value of |x|
 * CORRESPONDS TO CASE 8 IN THE MAIN PROGRAM
 *
*/

long double logP;

//if ((std::abs(x) <= b_max) && (std::abs(x) >= b_min)){ // Change on 17 March 2020
if ((std::abs(x) <= b_max) && (std::abs(x) >= b_min)){
	logP=-log(std::abs(b_max-b_min));
} else{
	logP=-INFINITY;
}
return logP;
}

long double logP_uniform_cos(const long double b_min, const long double b_max, const long double x){
/* 
 * Calculates the log probability for a uniform prior, given the value of cos(x)
 * x is assumed to be given in degree ==> conversion made in this prior
 * CAREFUL: THIS IS NOT GIVING A UNIFORM PRIOR IN COS(I) DUE TO THE MISSING JACOBIAN 
 * CORRESPONDS TO CASE 9 IN THE MAIN PROGRAM
 *
*/

long double logP;
long double x_rad;

x_rad=PIl * x / 180.;
//if ((cos(x_rad) <= b_max) && (cos(x_rad) >= b_min)){ // Change on 17 March 2020
if ((cos(x_rad) < b_max) && (cos(x_rad) > b_min)){
	logP=-log(std::abs(b_max-b_min));
} else{
	logP=-INFINITY;
}
return logP;
}

// -------------

long double logP_gaussian(const long double mean, const long double sigma,const  long double x){
/* 
 * Calculates the log probability for a gaussian at a given value x
 * CORRESPONDS TO CASE 2 IN THE MAIN PROGRAM
 *
*/

return -log(sqrt(2*PIl)*sigma)-0.5*pow((x- mean)/sigma, 2.);
}


// -------------

long double logP_multivariate_gaussian(const Eigen::VectorXd& mean,const  Eigen::MatrixXd& Matrix, const Eigen::VectorXd& x){
/* 
 * Calculates the log probability for a multivariate gaussian for a given vector x
 * CORRESPONDS TO CASE 3 IN THE MAIN PROGRAM
 *
*/
long double logP;

long double vnorm=(x - mean).transpose() * Matrix.inverse() * (x - mean);  // CHECK THAT IT GIVES A SCALAR

logP=-0.5*vnorm + log( pow((2*PIl), (x.size()/2)) * sqrt(Matrix.determinant())  );
return logP;
}

// -------------

long double logP_jeffrey(const long double hmin, const long double hmax, const long double h){
/* 
 * Calculates the log probability for a Jeffreys-truncated prior at a given value h
 * CORRESPONDS TO CASE 4 IN THE MAIN PROGRAM
 *
*/

long double logP;

if (h < hmax && h > 0){
	long double prior=1./(h+hmin);
	long double norm=log((hmax+hmin)/hmin);

	logP=log(prior/norm);
} else {
	logP=-INFINITY;
}
return logP;
}

// -------------

long double logP_jeffrey_abs(const long double hmin, const long double hmax, const long double h){
/* 
 * Calculates the log probability for a Jeffreys-truncated prior at a given value of abs(h)
 * CORRESPONDS TO CASE 4 IN THE MAIN PROGRAM
 *
*/

long double logP;

//if (std::abs(h) < hmax && h > 0){
if (std::abs(h) < hmax){
	long double prior=1./(std::abs(h)+hmin);
	long double norm=log((hmax+hmin)/hmin);

	logP=log(prior/norm);
} else {
	logP=-INFINITY;
}
return logP;
}

// -------------


long double logP_uniform_gaussian( const long double b_min, const long double b_max, 
			const long double sigma,const  long double x){
/* 
 * Calculates the log probability for a Uniform-Gaussian prior at a given value x
 * CORRESPONDS TO CASE 5 IN THE MAIN PROGRAM
 *
*/

long double logP, C;

if (x < b_min){
	logP=-INFINITY;
}
if ((x <= b_max) && (x >= b_min)){
	logP=0; // this is log(1)
}
if ( x > b_max) {
	logP=-0.5*pow((x - b_max)/sigma, 2.);
}
C=log(std::abs(b_max-b_min)+0.5*sqrt(2*PIl)*sigma); // The normalisation constant

return logP-C;
}

// -------------


long double logP_gaussian_uniform(const  long double b_min, const long double b_max, 
			const long double sigma, const long double x){
/* 
 * Calculates the log probability for a Gaussian-Uniform prior at a given value x
 * CORRESPONDS TO CASE 6 IN THE MAIN PROGRAM
 *
*/

long double logP, C;

if (x > b_max){
	logP=-INFINITY;
}
if ((x <= b_max) && (x >= b_min)){
	logP=0; // this is log(1)
}
if ( x < b_min) {
	logP=-0.5*pow((x - b_min)/sigma, 2.);
}
C=log(std::abs(b_max-b_min)+0.5*sqrt(2*PIl)*sigma); // The normalisation constant

return logP-C;
}

// -------------

long double logP_gaussian_uniform_gaussian(const  long double b_min, const long double b_max, 
			const long double sigma1, const long double sigma2, const long double x){
/* 
 * Calculates the log probability for a Gaussian-Uniform-Gaussian prior at a given value x
 * sigma1 is the standard deviation for x<b_min. sigma2 is the standard deviation for x>b_max.
 * CORRESPONDS TO CASE 7 IN THE MAIN PROGRAM
 *
*/

long double logP, C;

if (x < b_min){
	logP=-0.5*pow((x - b_min)/sigma1, 2.);
}
if ((x <= b_max) && (x >= b_min)){
	logP=0; // this is log(1)
}
if ( x > b_max) {
	logP=-0.5*pow((x - b_max)/sigma2, 2.);
}
C=log(std::abs(b_max-b_min)+0.5*sqrt(2*PIl)*(sigma1 + sigma2)); // The normalisation constant

return logP-C;
}


/*
//// Test function
int main(){

std::cout << " Test of all the prior functions" << std::endl;
std::cout << " Expected values were calculated using the IDL parent program, when necessary" << std::endl;
std::cout << " ---------------";

std::cout << "Test of logP_gaussian()..." << std::endl;
std::cout << "    - logP_gaussian(1.L, 0.5L, 4.L)= ";
std::cout << logP_gaussian(1L, 0.5L, 4L) << std::endl;
std::cout << "  ... Expected : -18.225791" << std::endl;
std::cout << " ---------------";

std::cout << "Test of logP_uniform()..." << std::endl;
std::cout << "    - logP_uniform(0, 1, -2)= ";
std::cout << logP_uniform(0, 1, -2) << std::endl;
std::cout << "  ... Expected : -INFINITY" << std::endl;
std::cout << "    - logP_uniform(0, 1, 0.3)= ";
std::cout << logP_uniform(0, 1, 0.3) << std::endl;
std::cout << "  ... Expected : 0" << std::endl;
std::cout << " ---------------";

std::cout << "Test of logP_multivariate_gaussian()..." << std::endl;
VectorXd MVect(2), VVect(2);
MatrixXd Mat(2,2);
Mat << 2, 0.1,
	  0.1, 4;
MVect << 1, 6.5;
VVect << 0.1, 8;
std::cout << "    - Input Covariance Matrix = " << std::endl;
std::cout << Mat << std::endl;
std::cout << "    - Input Mean Vector = " << MVect.transpose() << std::endl;
std::cout << "    - Input Value Vector = " << VVect.transpose() << std::endl;

std::cout << "    - logP_multivariate_gaussian(Mean Vector, Matrix, Value Vector)= ";
std::cout << logP_multivariate_gaussian( MVect, Mat, VVect) << std::endl;
std::cout << "    - Expected : 2.3757209" << std::endl;
std::cout << " ---------------";

std::cout << "Test of logP_jeffrey()" << std::endl;
std::cout << "    - logP_jeffrey(0.1, 1, -2)= ";
std::cout << logP_jeffrey(0.1, 1, -2) << std::endl;
std::cout << "  ... Expected : -INFINITY" << std::endl;
std::cout << "    - logP_jeffrey(0.1, 1, 0.5)= ";
std::cout << logP_jeffrey(0.1, 1, 0.5) << std::endl;
std::cout << "  ... Expected : -0.363766" << std::endl;
std::cout << " ---------------";

std::cout << "Test of logP_uniform_gaussian()" << std::endl;
std::cout << "    - logP_uniform_gaussian(bmin=1, bmax=2, sigma=0.1, -1)= ";
std::cout << logP_uniform_gaussian(1, 2, 0.1, -1) << std::endl;
std::cout << "  ... Expected : -INFINITY" << std::endl;
std::cout << "    - logP_uniform_gaussian(bmin=1, bmax=2, sigma=0.1, 1.5)= ";
std::cout << logP_uniform_gaussian(1, 2, 0.1, 1.5) << std::endl;
std::cout << "  ... Expected : -0.11807759" << std::endl;
std::cout << "    - logP_uniform_gaussian(bmin=1, bmax=2, sigma=0.1, 3)= ";
std::cout << logP_uniform_gaussian(1, 2, 0.1, 3) << std::endl;
std::cout << "  ... Expected : -50.118074" << std::endl;
std::cout << " ---------------";

std::cout << "Test of logP_gaussian_uniform()" << std::endl;
std::cout << "    - logP_gaussian_uniform(bmin=1, bmax=2, sigma=0.1, -1)= ";
std::cout << logP_gaussian_uniform(1, 2, 0.1, -1) << std::endl;
std::cout << "  ... Expected : -200.11806" << std::endl;
std::cout << "    - logP_gaussian_uniform(bmin=1, bmax=2, sigma=0.1, 1.5)= ";
std::cout << logP_gaussian_uniform(1, 2, 0.1, 1.5) << std::endl;
std::cout << "  ... Expected : -0.11807759" << std::endl;
std::cout << "    - logP_gaussian_uniform(bmin=1, bmax=2, sigma=0.1, 3)= ";
std::cout << logP_gaussian_uniform(1, 2, 0.1, 3) << std::endl;
std::cout << "  ... Expected : -INFINITY" << std::endl;
std::cout << " ---------------";

std::cout << "Test of logP_gaussian_uniform_gaussian()" << std::endl;
std::cout << "    - logP_gaussian_uniform_gaussian(bmin=1, bmax=2, sigma1=0.1, sigma2=0.4 -1)= ";
std::cout << logP_gaussian_uniform_gaussian(1, 2, 0.1, 0.4, -1) << std::endl;
std::cout << "  ... Expected : -200.48651" << std::endl;
std::cout << "    - logP_gaussian_uniform_gaussian(bmin=1, bmax=2, sigma1=0.1, sigma2=0.4, 1.5)= ";
std::cout << logP_gaussian_uniform_gaussian(1, 2, 0.1, 0.4, 1.5) << std::endl;
std::cout << "  ... Expected : -0.48652704" << std::endl;
std::cout << "    - logP_gaussian_uniform_gaussian(bmin=1, bmax=2, sigma1=0.1, sigma2=0.4, 3)= ";
std::cout << logP_gaussian_uniform_gaussian(1, 2, 0.1, 0.4, 3) << std::endl;
std::cout << "  ... Expected : -3.6115268" << std::endl;
std::cout << " ---------------";

}
*/
