/*
 * MALA.h
 *
 * Header file that contains all kind of class/structures
 * used to process and/or encapsulate data
 * 
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */

#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "model_def.h"
#include "outputs.h"
#include "config.h"
#include "diagnostics.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

//using namespace std;

int main_standalone();

class MALA{
	private:
		int seed; // for generating random numbers
		long initial_i;
		long Nsamples;
		long Nchains;
		long Nvars;
		std::vector<int> Nt_learn;
		std::vector<int> periods_learn;
		int dN_mixing;
		long double epsilon1;
		MatrixXd epsilon2; // small quantity that is added to a matrix
		long double A1;
		long double delta;
		long double delta_x;
		long double c0; // This is what defines the initial gamma
		long double gamma;
		long double lambda_temp; // defines how high in temperature we have to go. Used ony by set_temperature()
		long double target_acceptance;
		bool use_drift; // if set to 0, we do not use the MALA scheme but just the MH scheme
	public:
		VectorXd sigma; //, sigma_prev; //(Nchains);
		MatrixXd mu; //, mu_prev; //(Nchains,Nparams);
		VectorXd Tcoefs; // list of all temperatures
		MatrixXd **covarmat; //, **covarmat_prev; //(Nchains, Nvars, Nvars);

		MALA(Config *cfg); // the constructor
		~MALA();
		void destroy_3dMatrix(MatrixXd** m3d, const int depth);
		long double multinormal_logpdf(const VectorXd& Deltavars, const VectorXd& drift1, const MatrixXd& covmat1);
		long random_int_vals(int Nmin, int Nmax); // Generate a random number between Nmin and Nmax
		void init_proposal(const VectorXd, const std::vector<std::string> s, const std::vector<std::string> s_inerror, const VectorXd fracerr, const VectorXd offseterr);
		void restore_proposal(const VectorXd& vars, Config *cfg_class); // Restore the proposal parameters from a previous file
		void update_proposal(const VectorXd& vars, long double acceptance, int m); // generate a new proposal and update the class variables
		int parallel_tempering(Model_def *model);
		VectorXd new_prop_values(const VectorXd& vars, int m);
		VectorXd D_MALA();
		long double p1_fct(long double x);
		MatrixXd p2_fct(const MatrixXd& x);
		VectorXd p3_fct(const VectorXd& x);
		void read_proposal(std::string readwhat, int m); // Used to read either current values (readwhat="current") of previous values (readwhat="previous") of mu, sigma and covarmat
		void update_position_MH(Model_def *model_class, Model_def *model_propos, Data *data_struc, Config *cfg_class, int m);
		void execute(Model_def *model_current, Model_def *model_propose, Data *data_struc, Config *cfg_class, Outputs *out, Diagnostics *diags); // main function of the class. Execute a MALA MCMC process.
		int msg_handler(const std::string file, const std::string error_type, const std::string fct_name, const std::string arguments, const bool fatal);
		
};

