#include <fstream>
#include <iostream>
//#include <string>
//#include <vector>
//#include <Eigen/Dense>

class class_test{
	public:
		int value1;
		double value2;
		class_test(int v0, double v00);
		void fct1(int *v1, double *v2); // pass by reference
		void fct2(int v1, double v2); // pass by value
};

class_test::class_test(int v0, double v00){
	value1=v0;
	value2=v00;
}


void class_test::fct1(int *v11, double *v22){
	value1=*v11;
	value2=*v22;
}

void class_test::fct2(int v1, double v2){
	value1=v1;
	value2=v2;
}

int main(){
	int v01;
	double v02;
	int v011;
	double v022;

	class_test myclass(1, 3.3);

	std::cout << myclass.value1 << "  " << myclass.value2 << std::endl;
	
	v011=3;
	v022=1.55;
	myclass.fct1(&v011, &v022);
	std::cout << myclass.value1 << "  " << myclass.value2 << std::endl;

	v01=30;
	v02=15.5;
	myclass.fct2(v01, v02);
	std::cout << myclass.value1 << "  " << myclass.value2 << std::endl;
	
}
