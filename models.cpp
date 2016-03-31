# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
//#include <cmath>

using Eigen::VectorXd;
using Eigen::VectorXi;

double lin_interpol(VectorXd x, VectorXd y, double x_int);
VectorXd amplitude_ratio(int l, double beta);
VectorXd optimum_lorentzian_calc(const VectorXd x, VectorXd y, double H_l, double fc_l, double f_s, double a2, double a3, double gamma_l, const int l, VectorXd V, double step);
VectorXd harvey_like(const VectorXd noise_params, VectorXd x, VectorXd y, const int Nharvey);


VectorXd model_MS_Global(VectorXd params, VectorXi params_length, VectorXd x){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation

	int Nmax=params_length[0]; // Number of Heights
	int lmax=params_length[1]; // number of degree - 1
	int Nfl0=params_length[2]; // number of l=0 frequencies
	int Nfl1=params_length[3]; // number of l=1 frequencies
	int Nfl2=params_length[4]; // number of l=2 frequencies
	int Nfl3=params_length[5]; // number of l=3 frequencies
	int Nsplit=params_length[6]; // number of splitting parameters. Should be 3 for a global MS model (a1,a2,a3)
	int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model

	int Nf=Nfl0+Nfl1+Nfl2+Nfl3;

	double inclination;

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,a2,a3;

	int Nharvey;
	long cpt;

	std::ofstream file_out;

	/* -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	inclination=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise];

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(lmax >=2){
		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

	}
	if(lmax >=3){
		Vl3=std::abs(params[Nmax+2]);
		ratios_l3=amplitude_ratio(3, inclination);
	}

    
	fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
	Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

	a1=std::abs(params[Nmax + lmax + Nf]);
	a2=params[Nmax + lmax + Nf + 1];
	a3=params[Nmax + lmax + Nf + 2];

	model_final.setZero();
	
	//std::cout << "Here 2 " << std::endl;
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
				
		fl0=fl0_all[n];
		Hl0=std::abs(params[n]);
		Wl0=std::abs(Wl0_all[n]);
		//std::cout << "Hl0 " << Hl0 << std::endl;
		//std::cout << "fl0 " << fl0 << std::endl;
		//std::cout << "a1 " << a1 << std::endl;
		//std::cout << "a2 " << a2 << std::endl;
		//std::cout << "a3 " << a3 << std::endl;
		//std::cout << "Wl0 " << Wl0 << std::endl;
		//std::cout << "ratios_l0 " << ratios_l0 << std::endl;
		//std::cout << "step " << step << std::endl;
		model_final=optimum_lorentzian_calc(x, model_final, Hl0, fl0, a1, a2, a3, Wl0, 0, ratios_l0, step);
        //model_final=model_final + model_l0;
		
		//file_out.open("model_l0_MS_outputs.txt");
		//file_out << "x    trunc c++ model " << std::endl;
		//for(long i=0; i<x.size(); i++){
		//	file_out << x(i) << "    " << model_l0(i) << std::endl;
		//}
		//file_out.close();
		//std::cout << " Success. File closed. Stoping program here (press Enter to continue)" << std::endl;

		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Hl1=params[n]*Vl1;
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			model_final=optimum_lorentzian_calc(x, model_final, Hl1, fl1, a1, a2, a3, Wl1, 1, ratios_l1, step);
			//model_final=model_final + model_l1;
            //std::cout << "Hl1 " << Hl1 << std::endl;
            //std::cout << "fl1 " << fl1 << std::endl;
            //std::cout << "a1 " << a1 << std::endl;
            //std::cout << "a2 " << a2 << std::endl;
            //std::cout << "a3 " << a3 << std::endl;
            //std::cout << "Wl1 " << Wl1 << std::endl;
            //std::cout << "ratios_l1 " << ratios_l1 << std::endl;
            //std::cout << "step " << step << std::endl;
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Hl2=params[n]*Vl2;
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			model_final=optimum_lorentzian_calc(x, model_final, Hl2, fl2, a1, a2, a3, Wl2, 2, ratios_l2, step);
			//model_final=model_final + model_l2;
            //std::cout << "Hl2 " << Hl2 << std::endl;
            //std::cout << "fl2 " << fl2 << std::endl;
            //std::cout << "a1 " << a1 << std::endl;
            //std::cout << "a2 " << a2 << std::endl;
            //std::cout << "a3 " << a3 << std::endl;
            //std::cout << "Wl2 " << Wl2 << std::endl;
            //std::cout << "ratios_l2 " << ratios_l2 << std::endl;
            //std::cout << "step " << step << std::endl;
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Hl3=params[n]*Vl3;
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			model_final=optimum_lorentzian_calc(x, model_final, Hl3, fl3, a1, a2, a3, Wl3, 3, ratios_l3, step);
			//model_final=model_final + model_l3;
            //std::cout << "Hl3 " << Hl3 << std::endl;
            //std::cout << "fl3 " << fl3 << std::endl;
            //std::cout << "a1 " << a1 << std::endl;
            //std::cout << "a2 " << a2 << std::endl;
            //std::cout << "a3 " << a3 << std::endl;
            //std::cout << "Wl3 " << Wl3 << std::endl;
            //std::cout << "ratios_l3 " << ratios_l3 << std::endl;
            //std::cout << "step " << step << std::endl;
		}		
	}
    //std::cin.ignore();


	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	noise_params.array().abs();
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params, x, model_final, Nharvey); // this function increment the model_final with the noise background

	return model_final;
}

VectorXd model_Test_Gaussian(VectorXd params, VectorXi params_length, VectorXd x){
/*
 * A model in which we fit a Gaussian + White noise.
 * Parameters are assumed to be in that order: Maximum Height, variance/width, central frequency, White noise
*/
	VectorXd nu0(x.size()), model_final(x.size()), tmp(x.size());

	tmp.setConstant(params[3]);

	model_final= -0.5 * (x - nu0.setConstant(params[2])).array().square() /pow(params[1],2);
	model_final= params[0]*model_final.array().exp();
	model_final= model_final + tmp;
	return model_final;
}


