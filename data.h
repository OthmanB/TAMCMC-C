
/*
 * data.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#include <Eigen/Dense>
#include <string>
#include <vector>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct Data{
		VectorXd x; // In the case of a 1D fit (e.g. fit of Mass, [Fe/H], ...), this variable will be ignored by the model function.
		VectorXd y;
		VectorXd sigma_y;
		long Nx; // Ny is not checked but should be as long as x
		//vector<string> xlabel; // label for x-axis. In ASCII file, the axis labels are identified by the symbol ! at the begining of the line
		//vector<string> ylabel;
		//vector<string> xunit;
		//vector<string> yunit;
		string xlabel; // label for x-axis. In ASCII file, the axis labels are identified by the symbol ! at the begining of the line
		string ylabel;
		string xunit;
		string yunit;
		vector<string> header; // Any kind of information that is useful. e.g. if it is a test, model info + model inputs. Any header (In ASCII, marked by #), will be put there.
};

struct Data_Nd{
	MatrixXd data; // 2D array in which each columns are expected to contain the values of a given parameter
	vector<string> header;
	vector<string> labels;
	vector<string> units;
};
