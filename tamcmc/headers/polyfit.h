#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/QR>

Eigen::VectorXd polyfitXd(const Eigen::VectorXd &t, const Eigen::VectorXd &v, const int order);
std::vector<double> polyfit(const std::vector<double> &t, const std::vector<double> &v, const int order);
