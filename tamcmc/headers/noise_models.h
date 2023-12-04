/*
 * noise_models.cpp
 *
 *  Created on: 24 Feb 2016
 *      Author: obenomar
 */
#pragma once
#include <math.h>
#include <Eigen/Dense>

using Eigen::VectorXd;

VectorXd harvey_like(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey);
VectorXd harvey1985(const VectorXd& noise_params, const VectorXd& x, const VectorXd& y, const int Nharvey);
/**
 * @brief Compute eta squared using definition from Kallinger+2014 (Eq. 1)
 * 
 * 
 * @param x 
 * @return VectorXd containing eta^2
 */
VectorXd eta_squared_Kallinger2014(const VectorXd& x);


/**
 * @brief Calculates the norm ksi associated to the noise profile as per defined by Kallinger+2014.
 *
 * @param b The b parameter of the noise model.
 * @param c The c parameter of the noise model.
 * @param x The x-coordinates of the data points.
 * @return double The norm ksi.
 */
double get_ksinorm(const double b, const double c, const Eigen::VectorXd& x);

/**
 * @brief Calculates the noise background as in Kallinger+2014
 * 
 * This function calculates the noise background as in Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
 * It uses notations similar to their Table 2
 * Note that here we assume the instrumental noise to be Pinstrument(nu) = 0
 * Also note that one would need to add eta^2 to the Gaussian envelope if this noise is used for fitting
 * The noise_params must be have the parameters provided in this order:
 * 		- Granulation Amplitude a : ka, sa, t
 *      - Characteristic frequency ~ 1/timescale b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
 *      - Characteristic frequency ~ 1/timescale b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
 * Such that at the end we have: [ka,sa,t,k1,s1,c1, k2,s2,c2, N0]
 * @param numax The frequency at maximum power of the modes 
 * @param Mass The stellar mass of the star. If not known, you must set it to Mass=1 and also set t=1
 * @param noise_params The noise parameters structured as explained in the function help
 * @param x A vector that contains the frequencies of the data points
 * @param y An initial vector of the same size as x and that contains the power. This vector could be initialised to 0 or be initialised using an instrumental noise
 * @return VectorXd of the same size as x and y that contains the power of the noise background added to the initial y vector.
 * 
 * @note The function returns the noise background.
 */
VectorXd Kallinger2014(const double numax, const double mu_numax, const double Mass, const VectorXd& noise_params,const VectorXd& x, const VectorXd& y);
VectorXd Kallinger2014_V2(const double numax, const double mu_numax, const VectorXd& noise_params,const VectorXd& x, const VectorXd& y);