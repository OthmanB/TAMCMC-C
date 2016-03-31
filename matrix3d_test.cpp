#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>

int main(){
	long Nmat=2;
	long Nch=3;
	MatrixXd mat0(Nmat,Nmat), mat1(Nmat,Nmat), mat2(Nmat,Nmat);
//	MatrixXd **mat3d;
//	matrix3d matr3d(Nch, Nmat, Nmat);	
	
	mat0(0,0)=0.1;
	mat0(0,1)=1;
	mat0(1,0)=0.1;
	mat0(1,1)=1;

	mat1(0,0)=0.5;
	mat1(0,1)=2;
	mat1(1,0)=0.5;
	mat1(1,1)=2;

	mat2(0,0)=0.75;
	mat2(0,1)=3;
	mat2(1,0)=0.75;
	mat2(1,1)=3;


//	mat3d=initialize_3dMatrix(Nch, Nmat, Nmat);
	mat3d=matr3d.initialize_3dMatrix(Nch, Nmat, Nmat);
	matr3d.set_3dMatrix(mat3d, mat0, 0);
	matr3d.set_3dMatrix(mat3d, mat1, 1);
	matr3d.set_3dMatrix(mat3d, mat2, 2);

	std::cout << "Testing that the values of the submatrix are not passed by reference..." << std::endl;
	std::cout << "  - The value of m[0](*,*): " << std::endl;
	std::cout << *mat3d[0] << std::endl;
	std::cout << "  - we change mat1(0,0)=0.5 into mat1(0,0)=10 " << std::endl;
	mat1(0,0)=10.;
	std::cout << "  - The value of m[0](*,*) after the change (should be identical to before the change if values are hard copied): " << std::endl;
	std::cout << *mat3d[0] << std::endl;
	
	std::cout << "Testing that refreshing one submatrix does work" << std::endl;
	std::cout << "  - Setting m[0](*,*) to mat1 " << std::endl;
	matr3d.set_3dMatrix(mat3d, mat1, 0);
	std::cout << "  - New value of m[0](*,*): " << std::endl;
	std::cout << *mat3d[0] << std::endl;

	std::cout << " Values of the 3d matrix were reinitialised to: " << std::endl;
	matr3d.set_3dMatrix(mat3d, mat0, 0);
	matr3d.set_3dMatrix(mat3d, mat1, 1);
	matr3d.set_3dMatrix(mat3d, mat2, 2);
	for(int m=0; m<Nch; m++){
		std::cout << "pointer mat3d[" << m << "]=" << mat3d[m] << std::endl;
		std::cout << "Matrix mat3d[" << m << "]=" << std::endl;
		std::cout << *mat3d[m] << std::endl;
		std::cout << "-----" << std::endl;
	}
	std::cout << "-----" << std::endl;

}
