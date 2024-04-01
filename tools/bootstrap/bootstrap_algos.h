#pragma once
#include <iostream>
#include <vector>
#include <random>
//#include <cstdlib> // For rand and srand functions
//#include <ctime>   // For time function
#include <Eigen/Dense>

using Eigen::VectorXd;
using namespace std;

VectorXd blockBootstrap(const VectorXd& samples, int blockSize);
VectorXd Bootstrap(const VectorXd& samples);