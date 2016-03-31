/*
 * MALA.h
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
#include "model_def.h"
#include "outputs.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using namespace std;

class MALA{
	private:
		int seed; // for generating random numbers
		long Nsamples;
		long Nchains;
		long Nvars;
		vector<int> Nt_learn;
		vector<int> periods_learn;
		int dN_mixing;
		double epsilon1;
		MatrixXd epsilon2; // small quantity that is added to a matrix
		double A1;
		double delta;
		double delta_x;
		double c0; // This is what defines the initial gamma
		double gamma;
		double lambda_temp; // defines how high in temperature we have to go. Used ony by set_temperature()
		double target_acceptance;
		bool use_drift; // if set to 0, we do not use the MALA scheme but just the MH scheme
	public:
		VectorXd sigma, sigma_prev; //(Nchains);
		MatrixXd mu, mu_prev; //(Nchains,Nparams);
		VectorXd Tcoefs; // list of all temperatures
		MatrixXd **covarmat, **covarmat_prev; //(Nchains, Nvars, Nvars);
		MALA(long Nspl, long Nchs, VectorXd params, vector<string> params_names, VectorXi relax); // the constructor
		long double multinormal_logpdf(VectorXd Deltavars, VectorXd drift1, MatrixXd covmat1);
		long random_int_vals(int Nmin, int Nmax); // Generate a random number between Nmin and Nmax
		void init_proposal(VectorXd, vector<string> s);
		void update_proposal(VectorXd vars, double acceptance, int m); // generate a new proposal and update the class variables
		int parallel_tempering(Model_def *model);
		VectorXd new_prop_values(VectorXd vars, int m);
		VectorXd D_MALA();
		double p1_fct(double x, double epsilon, double A);
		MatrixXd p2_fct(MatrixXd x, double A);
		VectorXd p3_fct(VectorXd x, double A);
		void read_proposal(string readwhat, int m); // Used to read either current values (readwhat="current") of previous values (readwhat="previous") of mu, sigma and covarmat
		void update_position_MH(Model_def *model_class, Model_def *model_propos, Data *data_struc, Config *cfg_class, int m);
		void execute(Model_def *model_current, Model_def *model_propose, Data *data_struc, Config *cfg_class, Outputs *out); // main function of the class. Execute a MALA MCMC process.
};

