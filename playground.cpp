#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
//#include "stdlib.h"
#include "config.h"

using namespace std;

vector<string> tiny_test(vector<string> word1, string word2, string word3);

/*
// Read a list of words and store them into a vector of strings
int main() {
    vector<string> list;

    string word;
    ifstream wordfile("dic.dat");
    while (wordfile >> word) {
        list.push_back(word);
    }

    for (unsigned n=0; n<list.size(); ++n) {
        cout << list.at( n ) << endl;
    }


    return 0;
}
*/

int main() {
    vector<string> list, list2;

    string word1, word2, word3;

   word1="Hello";
   word2="How are you today?";
   word3="I_am_Fine_Thank_you";

   list.push_back(word1);
   list.push_back(word2);
   list.push_back(word3);

 //   for (unsigned n=0; n<list.size(); ++n) {
 //       cout << list.at( n ) << endl;
 //   }

   // Second test
   word2="More for this list";
   list2=tiny_test(list, word2, "I_agree");
    for (unsigned n=0; n<list2.size(); ++n) {
        cout << list2.at( n ) << endl;
    }

    // Another test on resizing VectorXd type
    cout << "Test on resizing vectorXd type..." << endl;
    int Nparams=10;
    Eigen::VectorXd vars(Nparams);
    for (int i=0; i<5; i++){
	vars[i]=1.;
    }
    vars.conservativeResize(5);
    for (unsigned n=0; n<5; ++n) {
        cout << vars[n] << endl;
    }

    // Test on resizing vector type... the way to proceed in order to keep the elements after resizing is different than for VectorXd
    cout << "Test on resizing vector type... the way to proceed in order to keep the elements after resizing is different than for VectorXd" << endl;
    vector<double> vars2;
    vars2.resize(10);
    for (int i=0; i<5; i++){
	vars2[i]=1.;
    }
    vars.resize(5);
    for (unsigned n=0; n<5; ++n) {
        cout << vars2[n] << endl;
    }
    cout << " ---- " << endl;

    // Another test
    Eigen::MatrixXd mat(2,2), mat2(1,2);

    mat.row(0) << 2,2;
    mat.row(1) << 0.1, 0.1;
    cout << " ---- " << endl;
    cout << mat << endl;
    cout << " ---- " << endl;
    cout << "mat.coeff(0,0)=" << mat.coeff(0,0) << endl;
    mat(0,0)=3;
    cout << "mat.coeff(0,0)=" << mat.coeff(0,0) << endl;
    cout << " ---- " << endl;
    cout << " minimum of a vector " << endl;
    Eigen::VectorXd v(3), v_from_mat;
    v << 4, 2, 3;
    cout << v.minCoeff() << endl;
    cout << " ---- " << endl;
    v_from_mat=mat.row(0);

    cout << "operator %: What is it doing in C++?" << endl;
    cout << " 10. / 3. = " << 10./3. << endl;
    cout << " 10 / 3 = " << 10/3 << endl;

    cout << " 10 % 1  = " << 10%1 << endl;
    cout << " 10 % 2  = " << 10%1 << endl;
    cout << " 10 % 3  = " << 10%3 << endl;
    cout << " 10 % 4  = " << 10%4 << endl;
    cout << " 10 % 5  = " << 10%5 << endl;
    cout << " 10 % 6  = " << 10%6 << endl;
    cout << " 10 % 7  = " << 10%7 << endl;
    cout << " 10 % 8  = " << 10%8 << endl;
    cout << " 10 % 9  = " << 10%9 << endl;
    cout << " 10 % 10  = " << 10%10 << endl;
    cout << " Note that the % operator can only be used with integers. In addition, the denominator cannot be 0" << endl;
    cout << " Returns >0 if the result has a numerand" << endl;


    cout << " --- Test when reading files with header ---" << endl;
    string file_in_data="/home/obenomar/Dropbox/Temporary/Cpp-Playground/data_tests/test-gauss.txt";
    Config config(file_in_data, "","");
    //char *line0;
    string line0;

    


    return 0;
}

vector<string> tiny_test(vector<string> list1, string word2, string word3){

	//vector<string> list;
	list1.push_back(word2);
        list1.push_back(word3);

return list1;
}

