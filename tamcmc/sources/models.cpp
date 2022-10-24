# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
#include "models.h"
#include "noise_models.h"
#include "../../external/ARMM/solver_mm.h"
#include "../../external/ARMM/bump_DP.h"
#include "../../external/integrate/activity.h"
#include "acoefs.h"
#include "interpol.h"
#include "linfit.h"
#include "polyfit.h"
#include "../../external/spline/src/spline.h"
//#include "linspace.h"
//#include "Qlm.h"
//#include <cmath>

using Eigen::VectorXd;
using Eigen::VectorXi;

//double lin_interpol(VectorXd x, VectorXd y, double x_int);
VectorXd amplitude_ratio(int l, double beta);

VectorXd model_MS_Global_a1l_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = M_PI; //3.141592653589793238462643383279502884L;

    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1(l=1),eta,a3, magb, magalfa, asym, a1(l=2))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    double inclination;
    double trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a11, a12,eta0, a3, asym;
    
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=10; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    
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
    
    a11=std::abs(params[Nmax + lmax + Nf]);
    a12=std::abs(params[Nmax + lmax + Nf+6]);
    //eta0=params[Nmax + lmax + Nf + 1]; // 7 Dec 2021: eta0 is not anymore a parameter
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];

    
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1l_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11 << std::endl;
    //std::cout << "a12=" << a12 << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    //outparams=1;
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }	
        model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl0, fl0, a11, a12, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        
        if(lmax >=1){
 			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=lin_interpol(fl0_all, Wl0_all, fl1);
            Wl1=std::abs(Wl1);
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
            if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a11, a12, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
           model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl1, fl1, a11, a12, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        }
        if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a11, a12, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl2, fl2, a11, a12, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        }
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a11, a12, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl3, fl3, a11, a12, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
  if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a11 / a12 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    return model_final;
}

VectorXd model_MS_Global_a1n_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be Nmax+6 for a ALL global MS model (<a1>0=0,eta,a3, magb, magalfa, asym, a1(n))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double inclination, trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a11, a12; //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym, a1_l1, a1_l2;
    
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=10; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
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
    
    a11=params.segment(Nmax + lmax + Nf + 6, Nmax);
    a12=a11;
    //eta0=params[Nmax + lmax + Nf + 1]; // 7 Dec 2021: eta0 is not anymore a parameter
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1n_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11.transpose() << std::endl;
    //std::cout << "a12=" << a12.transpose() << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
		a1_l1=std::abs(a11[n]);
		a1_l2=std::abs(a12[n]);
        model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl0, fl0, a1_l1, a1_l2, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   

        if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl1, fl1, a1_l1, a1_l2, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1, a1_l2, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
        }
        if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl2, fl2, a1_l1, a1_l2, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l1, a1_l2, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        //std::cout << "lmax=" << lmax << std::endl;
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl3, fl3, a1_l1, a1_l2, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l1, a1_l2, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1_l1 / a1_l2 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }    
    return model_final;
}

VectorXd model_MS_Global_a1n_a2a3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
        
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be Nmax+6 for a ALL global MS model (<a1>0=0,eta,a3, magb, magalfa, asym, a1(n))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double inclination, trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a11, a12, a2_terms(3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a2, a3, asym, a1_l1, a1_l2;
    
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=10; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
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
    
    a11=params.segment(Nmax + lmax + Nf + 6, Nmax);
    a12=a11;
    a2_terms=params.segment(Nmax + lmax + Nf + 6 + Nmax, Nmax); 
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1n_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11.transpose() << std::endl;
    //std::cout << "a12=" << a12.transpose() << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
            //std::cout << "[0] do conversion" << std::endl;
        } else{
            Hl0=std::abs(params[n]);
        }       
        a1_l1=std::abs(a11[n]);
        a1_l2=std::abs(a12[n]);
        a2=0;
        model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl0, fl0, a1_l1, a1_l2, 0, 0, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
                //std::cout << "[1] do conversion" << std::endl;
            } else{
                Hl1=std::abs(params[n]*Vl1);
            }
            //a2=a2_terms[0] + a2_terms[1]*fl1;
            model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl1, fl1, a1_l1, a1_l2, a2_terms[n], a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1, a1_l2, a2_terms[n], a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            if(do_amp){
                Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //std::cout << "[2] do conversion" << std::endl;
            } else{
                Hl2=std::abs(params[n]*Vl2);
            }  
            //a2=a2_terms[0] + a2_terms[1]*fl2;
            model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl2, fl2, a1_l1, a1_l2, a2_terms[n], a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l1, a1_l2, a2_terms[n], a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            if(do_amp){
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
            } else{
                Hl3=std::abs(params[n]*Vl3);            
            }      
            //a2=a2_terms[0] + a2_terms[1]*fl3; 
            model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl3, fl3, a1_l1, a1_l2, a2_terms[n], a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l1, a1_l2, a2_terms[n], a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

  if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1_l1 / a1_l2 / a2 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    std::cout << "    model_MS_Global_a1n_a2a3_HarveyLike() not tested yet" << std::endl;
    std::cout << "    Note also that the fact the we use a2_terms[n] for all l could be problematic because we approximate the nu ~ n  while it is nu ~ (n,l)" << std::endl;
    std::cout << "    An alternative is that the parameters are a2 at nu(l=0) (or any other reference frequency set of length Nmax" << std::endl;
    std::cout << "    The program will exit now" << std::endl;
    exit(EXIT_SUCCESS);
     
    return model_final;
}


VectorXd model_MS_Global_a1l_a2a3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 2*Nmax+6 for a ALL global MS model (<a1>0=0,eta,a3, magb, magalfa, asym, a1(n,l))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    double inclination, trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a2_terms(lmax*3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,a2,a3, asym, a11, a12;
    
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=10; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
 
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
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
    
    a11=std::abs(params[Nmax + lmax + Nf]);
    a12=std::abs(params[Nmax + lmax + Nf+6]);
    a2_terms=params.segment(Nmax + lmax + Nf + 7, lmax*3 ); // 3 parameters for each l... one constant term, one linear and one quadratic
    a3=params[Nmax + lmax + Nf + 2];
    
    asym=params[Nmax+lmax + Nf + 5];
 
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1nl_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11.transpose() << std::endl;
    //std::cout << "a12=" << a12.transpose() << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}				
        model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl0, fl0, 0, 0, 0, 0, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }           
        if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}	
            a2=a2_terms[0] + a2_terms[1]*fl1 + a2_terms[2]*fl1*fl1;			
            model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl1, fl1, a11, a12, a2, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a11, a12, a2, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}
            a2=a2_terms[0] + a2_terms[1]*fl2 + a2_terms[2]*fl2*fl2; 	
            model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl2, fl2, a11, a12, a2, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a11, a12, a2, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        if(lmax >=3){
 			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
           a2=a2_terms[0] + a2_terms[1]*fl3 + a2_terms[2]*fl3*fl3;  
           model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl3, fl3, a11, a12, a2, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a11, a12, a2, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    

        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

  if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a11 / a12 / a2 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
 
    std::cout << "    model_MS_Global_a1l_a2a3_HarveyLike() not tested yet" << std::endl;
    std::cout << "    The program will exit now" << std::endl;
    exit(EXIT_SUCCESS);
   
    return model_final;
}

VectorXd model_MS_Global_a1nl_a2a3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
        
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 2*Nmax+6 for a ALL global MS model (<a1>0=0,eta,a3, magb, magalfa, asym, a1(n,l))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    double inclination, trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a11, a12, a2_terms; //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,a2,a3, asym, a1_l1, a1_l2;
    
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=10; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
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
    
    a11=params.segment(Nmax + lmax + Nf + 6, Nmax);
    a12=params.segment(Nmax + lmax + Nf + 6 + Nmax, Nmax);
    a2_terms=params.segment(Nmax + lmax + Nf + 6 + Nmax + Nmax , Nmax); // after a11 and a12
    a3=params[Nmax + lmax + Nf + 2];
    
    asym=params[Nmax+lmax + Nf + 5];
 
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1nl_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11.transpose() << std::endl;
    //std::cout << "a12=" << a12.transpose() << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
            //std::cout << "[0] do conversion" << std::endl;
        } else{
            Hl0=std::abs(params[n]);
        }       
        a1_l1=std::abs(a11[n]);
        a1_l2=std::abs(a12[n]);
        
        model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl0, fl0, 0, 0, 0, 0, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }           
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
                //std::cout << "[1] do conversion" << std::endl;
            } else{
                Hl1=std::abs(params[n]*Vl1);
            }               
            model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl1, fl1, a1_l1, a1_l2, a2_terms[n], a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1, a1_l2, a2_terms[n], a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            if(do_amp){
                Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //std::cout << "[2] do conversion" << std::endl;
            } else{
                Hl2=std::abs(params[n]*Vl2);
            }   
            model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl2, fl2, a1_l1, a1_l2, a2_terms[n], a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l1, a1_l2, a2_terms[n], a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            if(do_amp){
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
            } else{
                Hl3=std::abs(params[n]*Vl3);            
            }       
           model_final=optimum_lorentzian_calc_a1l_a2a3(x, model_final, Hl3, fl3, a1_l1, a1_l2, a2_terms[n], a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l1, a1_l2, a2_terms[n],  a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1_l1 / a1_l2 / a2 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
   
    std::cout << "    model_MS_Global_a1nl_a2a3_HarveyLike() not tested yet" << std::endl;
    std::cout << "    Note also that the fact the we use a2_terms[n] for all l could be problematic because we approximate the nu ~ n  while it is nu ~ (n,l)" << std::endl;
    std::cout << "    An alternative is that the parameters are a2 at nu(l=0) (or any other reference frequency set of length Nmax" << std::endl;
    std::cout << "    The program will exit now" << std::endl;
    exit(EXIT_SUCCESS);

    return model_final;
}


VectorXd model_MS_Global_a1nl_etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax.
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */
    
    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
        
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 2*Nmax+6 for a ALL global MS model (<a1>0=0,eta,a3, magb, magalfa, asym, a1(n,l))
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    double inclination, trunc_c;
    
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a11, a12;; //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym, a1_l1, a1_l2;
    
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=10; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    /*
     -------------------------------------------------------
     ------- Gathering information about the modes ---------
     -------------------------------------------------------
     */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
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
    
    a11=params.segment(Nmax + lmax + Nf + 6, Nmax);
    a12=params.segment(Nmax + lmax + Nf + 6 + Nmax, Nmax);
    //eta0=params[Nmax + lmax + Nf + 1]; // 7 Dec 2021: eta0 is not anymore a parameter
    eta0=eta0_fct(fl0_all);

    a3=params[Nmax + lmax + Nf + 2];
    
    asym=params[Nmax+lmax + Nf + 5];
 
    model_final.setZero();
    
    //std::cout << "model_MS_Global_a1nl_etaa3_HarveyLike" << std::endl;
    //std::cout << "a11=" << a11.transpose() << std::endl;
    //std::cout << "a12=" << a12.transpose() << std::endl;
    /* -------------------------------------------------------
     --------- Computing the models for the modes  ---------
     -------------------------------------------------------
     */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
            //std::cout << "[0] do conversion" << std::endl;
        } else{
            Hl0=std::abs(params[n]);
        }       
        a1_l1=std::abs(a11[n]);
        a1_l2=std::abs(a12[n]);
        
        model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl0, fl0, a1_l1, a1_l2, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }           
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
                //std::cout << "[1] do conversion" << std::endl;
            } else{
                Hl1=std::abs(params[n]*Vl1);
            }               
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl1, fl1, a1_l1, a1_l2, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1, a1_l2, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            if(do_amp){
                Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //std::cout << "[2] do conversion" << std::endl;
            } else{
                Hl2=std::abs(params[n]*Vl2);
            }   
            model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl2, fl2, a1_l1, a1_l2, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l1, a1_l2, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            if(do_amp){
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
            } else{
                Hl3=std::abs(params[n]*Vl3);            
            }       
           model_final=optimum_lorentzian_calc_a1l_etaa3(x, model_final, Hl3, fl3, a1_l1, a1_l2, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l1, a1_l2, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
    }
    //std::cin.ignore();
    
    /* -------------------------------------------------------
     ------- Gathering information about the noise ---------
     -------------------------------------------------------
     */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    
    /* -------------------------------------------------------
     ---------- Computing the mode of the noise ------------
     -------------------------------------------------------
     */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1_l1 / a1_l2 / eta0 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }    
    return model_final;
}


VectorXd model_MS_Global_a1etaa3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

	double inclination;
	

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;

	int Nharvey;
	long cpt;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited

	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	
  	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);
    
    inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
    inclination=inclination*180./pi;

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

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

	//a1=std::abs(params[Nmax + lmax + Nf]);
    //eta0=params[Nmax + lmax + Nf + 1]; // 7 Dec 2021: eta0 is not anymore a parameter
    eta0=eta0_fct(fl0_all);

	a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
		
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
			
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);

		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				//Hl1=std::abs(params[n]/(pi*Wl1));
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
            if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //Hl2=std::abs(params[n]/(pi*Wl2));
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				//Hl3=std::abs(params[n]/(pi*Wl3));
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}	
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

  if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }

	return model_final;
}

VectorXd model_MS_Global_a1a2a3_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = M_PI;    
    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double inclination;
    

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise), a2_terms(3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,a2,a3, asym;

    int Nharvey;
    long cpt;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited

    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    
    a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);
    
    inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
    inclination=inclination*180./pi;

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

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

    //a1=std::abs(params[Nmax + lmax + Nf]);
    a2_terms=params.segment(Nmax + lmax + Nf + 6,3); // a2 terms are after the a3 and asym terms  //two terms: one constant term + one linear in nu + one quadratic in nu
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    //std::cout << "[" << Nmax + lmax + Nf + 6 << "]  a2_terms = " << a2_terms.transpose() << std::endl;
    /*
    std::cout << "inclination = " << inclination << std::endl;
    std::cout << "a2_terms = " << a2_terms << std::endl;
    std::cout << "a3 = " << a3 << std::endl;
    std::cout << "asym = " << asym << std::endl;
    */
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    cpt=0;
    for(long n=0; n<Nmax; n++){
        
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);
            
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
        } else{
            Hl0=std::abs(params[n]);
        }       
        model_final=optimum_lorentzian_calc_a1a2a3(x, model_final, Hl0, fl0, 0, 0, 0, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            if(do_amp){
                //Hl1=std::abs(params[n]/(pi*Wl1));
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
            } else{
                Hl1=std::abs(params[n]*Vl1);
            }   
            a2=a2_terms[0] + a2_terms[1]*(fl1*1e-3) + a2_terms[2]*(fl1*fl1*1e-6); //two terms: one constant term + one linear in nu, after a11, set in mHz            
            model_final=optimum_lorentzian_calc_a1a2a3(x, model_final, Hl1, fl1, a1, a2, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, a2, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            } 
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            if(do_amp){
                Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //Hl2=std::abs(params[n]/(pi*Wl2));
            } else{
                Hl2=std::abs(params[n]*Vl2);
            }   
            a2=a2_terms[0] + a2_terms[1]*(fl2*1e-3) + a2_terms[2]*(fl2*fl2*1e-6); //two terms: one constant term + one linear in nu, after a11 
            model_final=optimum_lorentzian_calc_a1a2a3(x, model_final, Hl2, fl2, a1, a2, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, a2, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }      
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            if(do_amp){
                //Hl3=std::abs(params[n]/(pi*Wl3));
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
            } else{
                Hl3=std::abs(params[n]*Vl3);            
            }       
            a2=a2_terms[0] + a2_terms[1]*(fl3*1e-3) + a2_terms[2]*(fl3*fl3*1e-6); //two terms: one constant term + one linear in nu, after a11  
            model_final=optimum_lorentzian_calc_a1a2a3(x, model_final, Hl3, fl3, a1, a2, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, a2, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            } 
        }       
    }
    //std::cin.ignore();

    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    //std::cout << "noise_params = " << noise_params.transpose() << std::endl;
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
  if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / a2 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }

    return model_final;
}



VectorXd model_MS_Global_aj_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = M_PI;    
    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    //const std::std::vector<int> Nfl{params_length[2], params_length[3], params_length[4], params_length[5]};
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. 2 parameters for each aj. j={1,2,3,4,5,6} >> 12 aj parameters + 1 asym
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    VectorXd a1_terms(2), a2_terms(2), a3_terms(2),a4_terms(2),a5_terms(2),a6_terms(2);

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise);
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,a2,a3, a4, a5,a6, eta0, asym;

    int Nharvey;

    const int Nrows=200, Ncols=13; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=params[Nmax + lmax + Nf+Nsplit+Nwidth+Nnoise]; 

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
    Hl0_all=params.segment(0, Nmax);

    a1_terms=params.segment(Nmax + lmax + Nf, 2);
    a2_terms=params.segment(Nmax + lmax + Nf + 2,2); // aj and asym terms  //two terms: one constant term + one linear in nu
    a3_terms=params.segment(Nmax + lmax + Nf + 4,2);
    a4_terms=params.segment(Nmax + lmax + Nf + 6,2);
    a5_terms=params.segment(Nmax + lmax + Nf + 8,2);
    a6_terms=params.segment(Nmax + lmax + Nf + 10,2);
    asym=params[Nmax+lmax + Nf + 13];

    if (params[Nmax+lmax + Nf + 12] == 1){
        eta0=eta0_fct(fl0_all); 
    } else{
        eta0=0; // We decide here to include eta0 inside a2 term. eta0 effect can always be removed a posteriori when doing the aj decomposition
    }    
    model_final.setZero();
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    //outparams=1;    
    for(long n=0; n<Nfl0; n++){        
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);     
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
        } else{
            Hl0=std::abs(params[n]);
        }       
        model_final=optimum_lorentzian_calc_aj(x, model_final, Hl0, fl0, 0, 0, 0, 0, 0, 0, 0, asym, Wl0, 0, ratios_l0, step, trunc_c);
       if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0,  0, 0, 0, 0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }
    }   
    for(long n=0; n<Nfl1; n++){        
        fl1=params[Nmax+lmax+Nfl0+n];
        Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
        if(do_amp){
            //Hl1=std::abs(params[n]/(pi*Wl1));
            //Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
            Hl1=std::abs(lin_interpol(fl0_all, Hl0_all, fl1)/(pi*Wl1)*Vl1); // Modified on 1 Mar 2022
        } else{
            //Hl1=std::abs(params[n]*Vl1);
            Hl1=std::abs(lin_interpol(fl0_all, Hl0_all, fl1)); // Modified on 1 Mar 2022
        }
        a1=a1_terms[0] + a1_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        a2=a2_terms[0] + a2_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a3=a3_terms[0] + a3_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a4=a4_terms[0] + a4_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a5=a5_terms[0] + a5_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a6=a6_terms[0] + a6_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        model_final=optimum_lorentzian_calc_aj(x, model_final, Hl1, fl1, a1, a2, a3,a4, a5,a6,eta0,asym, Wl1, 1, ratios_l1, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, a2, a3, a4, a5,a6, eta0, asym, inclination;// mode_vec;
            Line=Line+1;
        } 
    }
    for(long n=0; n<Nfl2; n++){        
        fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
        Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
        if(do_amp){
            Hl2=std::abs(lin_interpol(fl0_all, Hl0_all, fl2)/(pi*Wl2)*Vl2);  // Modified on 1 Mar 2022
            //Hl2=std::abs(params[n]/(pi*Wl2));
        } else{
            //Hl2=std::abs(params[n]*Vl2);
            Hl2=std::abs(lin_interpol(fl0_all, Hl0_all, fl2)*Vl2); // Modified on 1 Mar 2022
        }   
        a1=a1_terms[0] + a1_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        a2=a2_terms[0] + a2_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a3=a3_terms[0] + a3_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a4=a4_terms[0] + a4_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a5=a5_terms[0] + a5_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a6=a6_terms[0] + a6_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        model_final=optimum_lorentzian_calc_aj(x, model_final, Hl2, fl2, a1, a2, a3,a4,a5,a6, eta0,asym, Wl2, 2, ratios_l2, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, a2, a3, a4, a5,a6, eta0, asym, inclination;// mode_vec;
            Line=Line+1;
        }      
    }
    for(long n=0; n<Nfl3; n++){        
        fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
        Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
        if(do_amp){
            //Hl3=std::abs(params[n]/(pi*Wl3));
            //Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
            Hl3=std::abs(lin_interpol(fl0_all, Hl0_all, fl3)/(pi*Wl3)*Vl3);  // Modified on 1 Mar 2022
        } else{
            //Hl3=std::abs(params[n]*Vl3);            
            Hl3=std::abs(lin_interpol(fl0_all, Hl0_all, fl3)*Vl3); // Modified on 1 Mar 2022
        }       
        a1=a1_terms[0] + a1_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        a2=a2_terms[0] + a2_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a3=a3_terms[0] + a3_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a4=a4_terms[0] + a4_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a5=a5_terms[0] + a5_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a6=a6_terms[0] + a6_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        model_final=optimum_lorentzian_calc_aj(x, model_final, Hl3, fl3, a1, a2, a3, a4, a5,a6, eta0, asym, Wl3, 3, ratios_l3, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, a2, a3, a4, a5,a6, eta0, asym, inclination;// mode_vec;
            Line=Line+1;
        } 
    }       
    //std::cin.ignore();    
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
    //std::cout << "noise_params = " << noise_params.transpose() << std::endl;
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
  if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / a2 / a3 / a4 / a5 / a6 / eta0 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
//    std::cout << "FINAL" << std::endl;
   // exit(EXIT_FAILURE);
    return model_final;
}




VectorXd model_MS_Global_ajAlm_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams) // Added on 31 Mar 2021
    {
    /* Model of the power spectrum of a Main sequence solar-like star
     * Make use of Gizon 2002, AN 323, 251 for describing the perturbation from Active Region on a2
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)    
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = M_PI; //3.141592653589793238462643383279502884L;  
 
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const int decompose_Alm=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2]; // TO VERIFY THE SLOT
    const std::string filter_type="gate";

    double inclination;
    
    VectorXd xfit, rfit, aj(6);
    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());
    VectorXd a1_terms(2), a3_terms(2),a5_terms(2);
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), epsilon_terms(2), thetas(2); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1, a3, a5, eta0, epsilon_nl, asym;
    double rho;

    int Nharvey;

    const int Nrows=1000, Ncols=16; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    aj.setZero();
    outparams=true;
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    
    a1_terms=params.segment(Nmax + lmax + Nf, 2);
    a3_terms=params.segment(Nmax + lmax + Nf + 2,2);
    a5_terms=params.segment(Nmax + lmax + Nf + 4,2);
    epsilon_terms=params.segment(Nmax + lmax + Nf + 6,2); // This is the intensity of active region (in the Sun it is around 10^-3). epsilon terms are after the a3 and asym terms  //three terms: one constant term + one linear in nu + one quadratic in nu
    thetas=params.segment(Nmax + lmax + Nf + 8,2)*M_PI/180.; // Extension of the active regions.
    asym=params[Nmax+lmax + Nf + 11];    
    inclination=params[Nmax + lmax + Nf + Nsplit + Nwidth + Nnoise];
    
    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);
    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }

    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);
    Hl0_all=params.segment(0, Nmax);

    model_final.setZero();

    // --- Large separation and centrifugal force ---
    if (params[Nmax+lmax + Nf + 10] == 1){ // IT IS HIGHLY RECOMMENDED TO ALWAYS KEEP IT TO 1 HERE FOR THIS MODEL AS a2 IS NOT A PARAMETER (contrary to the aj model)
        eta0=eta0_fct(fl0_all); 
    } else{
        eta0=0; // We decide here to include eta0 inside a2 term. eta0 effect can always be removed a posteriori when doing the aj decomposition
    }   
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    for(long n=0; n<Nfl0; n++){        
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);     
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
        } else{
            Hl0=std::abs(params[n]);
        }
        model_final=optimum_lorentzian_calc_aj(x, model_final, Hl0, fl0, 0, 0, 0, 0, 0, 0, 0, asym, Wl0, 0, ratios_l0, step, trunc_c);
       if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, asym, inclination;
            Line=Line+1;
        }
    }   
        
    for(long n=0; n<Nfl1; n++){        
        fl1=params[Nmax+lmax+Nfl0+n];
        Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
        if(do_amp){
            Hl1=std::abs(lin_interpol(fl0_all, Hl0_all, fl1)/(pi*Wl1)*Vl1); // Modified on 1 Mar 2022
        } else{
            Hl1=std::abs(lin_interpol(fl0_all, Hl0_all, fl1)); // Modified on 1 Mar 2022
        }
        a1=a1_terms[0] + a1_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        a3=a3_terms[0] + a3_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a5=a5_terms[0] + a5_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        epsilon_nl=epsilon_terms[0] + epsilon_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz            
        //eta0=0; // FOR TEST PURPOSE ONLY
        switch (decompose_Alm){
            case -1:
                model_final=optimum_lorentzian_calc_ajAlm(x, model_final, Hl1, fl1, a1, a3, a5, eta0, epsilon_nl, thetas, asym, Wl1, 1, ratios_l1, step, trunc_c);
                break;
            // do_a2, do_a4, do_a6
            case 0: // decompose_Alm = 0 ==> do_a2=true, do_a4=true, do_a6=true
                aj=decompose_Alm_fct(1, fl1, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                //std::cout << "l=1   thetas=" << thetas.transpose()  << "    epsilon_nl=" << epsilon_nl << "   odds aj=" << aj.transpose() << std::endl;
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl1, fl1, a1, aj[1], a3, aj[3], a5,aj[5], eta0, asym, Wl1, 1, ratios_l1, step, trunc_c);   
                break; 
            case 1:  // decompose_Alm = 1 ==> do_a2=true, do_a4=true, do_a6=false
                aj=decompose_Alm_fct(1, fl1, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                aj[5]=0;
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl1, fl1, a1, aj[1], a3, aj[3], a5, 0, eta0, asym, Wl1, 1, ratios_l1, step, trunc_c);   
                break; 
            case 2:  // decompose_Alm = 1 ==> do_a2=true, do_a4=false, do_a6=false
                aj=decompose_Alm_fct(1, fl1, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                aj[3]=0; aj[5]=0;
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl1, fl1, a1, aj[1], a3, 0, a5, 0, eta0, asym, Wl1, 1, ratios_l1, step, trunc_c);   
                break; 
            default:
                std::cout << " Error in model_MS_Global_ajAlm_HarveyLike: decompose_Alm must be set to -1 (use Alm directly), 0 (do_a2=true, do_a4=true, do_a6=true), 1 (do_a2=true, do_a4=true, do_a6=fasle) or 2 (do_a2=true, do_a4=false, do_a6=false) " << std::endl;
                exit(EXIT_FAILURE);
        }
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, aj[1], a3, aj[3], a5, aj[5], eta0*1e-6, epsilon_nl, thetas[0], thetas[1], asym, inclination;// mode_vec;
            Line=Line+1;
        } 
   }

    for(long n=0; n<Nfl2; n++){        
        fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
        Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
        if(do_amp){
            Hl2=std::abs(lin_interpol(fl0_all, Hl0_all, fl2)/(pi*Wl2)*Vl2);  // Modified on 1 Mar 2022
        } else{
            Hl2=std::abs(lin_interpol(fl0_all, Hl0_all, fl2)*Vl2); // Modified on 1 Mar 2022
        }   
        a1=a1_terms[0] + a1_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        a3=a3_terms[0] + a3_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a5=a5_terms[0] + a5_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        epsilon_nl=epsilon_terms[0] + epsilon_terms[1]*(fl2*1e-3); //two terms: one constant term + one linear in nu, after a11 
        switch (decompose_Alm){
            case -1:
                model_final=optimum_lorentzian_calc_ajAlm(x, model_final, Hl2, fl2, a1, a3, a5, eta0, epsilon_nl, thetas, asym, Wl2, 2, ratios_l2, step, trunc_c);
                break;
            // do_a2, do_a4, do_a6
            case 0: // decompose_Alm = 0 ==> do_a2=true, do_a4=true, do_a6=true
                aj=decompose_Alm_fct(2, fl2, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl2, fl2, a1, aj[1], a3, aj[3], a5,aj[5], eta0, asym, Wl2, 2, ratios_l2, step, trunc_c);   
                break; 
            case 1:  // decompose_Alm = 1 ==> do_a2=true, do_a4=true, do_a6=false
                aj=decompose_Alm_fct(2, fl2, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                aj[5]=0;
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl2, fl2, a1, aj[1], a3, aj[3], a5, 0, eta0, asym, Wl2, 2, ratios_l2, step, trunc_c);   
                break; 
            case 2:  // decompose_Alm = 1 ==> do_a2=true, do_a4=false, do_a6=false
                aj=decompose_Alm_fct(2, fl2, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                aj[3]=0; aj[5]=0;
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl2, fl2, a1, aj[1], a3, 0, a5, 0, eta0, asym, Wl2, 2, ratios_l2, step, trunc_c);   
                break; 
            default:
                std::cout << " Error in model_MS_Global_ajAlm_HarveyLike: decompose_Alm must be set to -1 (use Alm directly), 0 (do_a2=true, do_a4=true, do_a6=true), 1 (do_a2=true, do_a4=true, do_a6=fasle) or 2 (do_a2=true, do_a4=false, do_a6=false) " << std::endl;
                exit(EXIT_FAILURE);
        }
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, aj[1], a3, aj[3], a5, aj[5], eta0*1e-6, epsilon_nl, thetas[0], thetas[1], asym, inclination;// mode_vec;
            Line=Line+1;
        }      
    }
  
    for(long n=0; n<Nfl3; n++){        
        fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
        Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
        if(do_amp){
            Hl3=std::abs(lin_interpol(fl0_all, Hl0_all, fl3)/(pi*Wl3)*Vl3);  // Modified on 1 Mar 2022
        } else{
            Hl3=std::abs(lin_interpol(fl0_all, Hl0_all, fl3)*Vl3); // Modified on 1 Mar 2022
        }       
        a1=a1_terms[0] + a1_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        a3=a3_terms[0] + a3_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        a5=a5_terms[0] + a5_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        epsilon_nl=epsilon_terms[0] + epsilon_terms[1]*(fl3*1e-3); //two terms: one constant term + one linear in nu, after a11  
        switch (decompose_Alm){
            case -1:
                model_final=optimum_lorentzian_calc_ajAlm(x, model_final, Hl3, fl3, a1, a3, a5, eta0, epsilon_nl, thetas, asym, Wl3, 3, ratios_l3, step, trunc_c);
                break;
            // do_a2, do_a4, do_a6
            case 0: // decompose_Alm = 0 ==> do_a2=true, do_a4=true, do_a6=true
                aj=decompose_Alm_fct(3, fl3, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl3, fl3, a1, aj[1], a3, aj[3], a5,aj[5], eta0, asym, Wl3, 3, ratios_l3, step, trunc_c);   
                break; 
            case 1:  // decompose_Alm = 1 ==> do_a2=true, do_a4=true, do_a6=false
                aj=decompose_Alm_fct(3, fl3, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                aj[5]=0;
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl3, fl3, a1, aj[1], a3, aj[3], a5, 0, eta0, asym, Wl3, 3, ratios_l3, step, trunc_c);   
                break; 
            case 2:  // decompose_Alm = 1 ==> do_a2=true, do_a4=false, do_a6=false
                aj=decompose_Alm_fct(3, fl3, eta0, a1, epsilon_nl, thetas, filter_type); // Generate fl1(m) and decompose it to get aj coefficients. decompose_Alm contains the rule to apply: do we use a2? a4? a6? )
                aj[3]=0; aj[5]=0;
                model_final=optimum_lorentzian_calc_aj(x, model_final, Hl3, fl3, a1, aj[1], a3, 0, a5, 0, eta0, asym, Wl3, 3, ratios_l3, step, trunc_c);   
                break; 
            default:
                std::cout << " Error in model_MS_Global_ajAlm_HarveyLike: decompose_Alm must be set to -1 (use Alm directly), 0 (do_a2=true, do_a4=true, do_a6=true), 1 (do_a2=true, do_a4=true, do_a6=fasle) or 2 (do_a2=true, do_a4=false, do_a6=false) " << std::endl;
                exit(EXIT_FAILURE);
        }
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, aj[1], a3, aj[3], a5, aj[5], eta0*1e-6, epsilon_nl, thetas[0], thetas[1], asym, inclination;// mode_vec;
            Line=Line+1;
        } 
    }       
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
  if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "#Note: decompose_Alm=-1 : Alm directly, 0: do a2, a4, a6, +1: do a2, a4, +2: do a2 only \n# Input mode parameters. degree / freq / H / W / a1  / a2 /  a3  / a4 / a5 / eta0(10^6) / epsilon_nl / theta0 (rad) /  delta (rad) / asymetry / inclination";
        //std::string name_params = "# Input mode parameters. degree / freq / H / W / a1  / a2 /  a3  / a4 / a5 / eta0(10^6) / epsilon_nl / theta0 (rad) /  delta (rad) / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    
    //exit(EXIT_SUCCESS);
    return model_final;
}

VectorXd model_MS_Global_a1etaa3_Harvey1985(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

	double inclination, trunc_c;

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;

	int Nharvey;
	long cpt;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];

	inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);

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

	//a1=std::abs(params[Nmax + lmax + Nf]);
	//eta0=params[Nmax + lmax + Nf + 1];
    eta0=eta0_fct(fl0_all);
	a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	std::cout << "a1=" << a1 << std::endl;
	std::cout << "eta0=" << eta0 << std::endl;
	std::cout << "a3=" << a3 << std::endl;
	std::cout << "asym=" << asym << std::endl;
	exit(EXIT_SUCCESS);
	
	model_final.setZero();
	
	//std::cout << "Here 2 " << std::endl;
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
						
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);
			
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }     
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
    	}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	//noise_params.array().abs();
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey1985(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1/ eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
	return model_final;
}

// Last modifcation: on 05/02/2020: Handling the possibility that the user asked for Heights Hlm to be fitted instead of inclination
// This might be a temporary fix for testing. If successfull, a full derivation into a specific model might be better
// Added on 05/02/2020: Handling the possibility that the user asked for Heights Hlm to be fitted instead of inclination
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

	double inclination, trunc_c;

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;

	int Nharvey;
	long cpt;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];

    inclination=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise]; // This version of the code take inclination as direct argument

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
	//eta=params[Nmax + lmax + Nf + 1];
	eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
				
		fl0=fl0_all[n];
		Wl0=std::abs(Wl0_all[n]);	
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
			//std::cout << "[0] do conversion" << std::endl;
		} else{
			Hl0=std::abs(params[n]);
		}		
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }    		
		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
			if(do_amp){
				Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
				//std::cout << "[1] do conversion" << std::endl;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
 		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
				//std::cout << "[2] do conversion" << std::endl;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    		
        }
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
			if(do_amp){
				Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
		}		
	}

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }	

	return model_final;
}


// Added on 10 Feb 2020
// Alternative version  of the legacy function that is not dealing with the inclination
// but instead fit the mode heights, given a visibility
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Equal to the  number of heights for m components: Sum(l*(l+1))_l={1,2,3}
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double trunc_c;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;

    int Nharvey;
    long cpt;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    double inclination = -1;  // Dummy value here
    if (outparams){
        std::cout << " ERROR: outparams is not implemented for  " << __func__ << std::endl;
        std::cout << "       If you want to get outputs in ascii format, you must review this function and implement outparams = true" << std::endl;
        std::cout << "       The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    } 

    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        //std::cout << "Vl2 " << Vl2 << std::endl;

    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
    }

    // l=1 Heights are defined by two parameters 
    ratios_l0[0]=1;

    ratios_l1[0]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+1]);  // m=-1
    ratios_l1[1]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise]); // m=0
    ratios_l1[2]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+1]); // m=+1 

    // l=2 Heights are defined by three parameters  
    ratios_l2[0]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+4]);  // m=-2
    ratios_l2[1]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+3]); // m=-1
    ratios_l2[2]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+2]); // m=0 
    ratios_l2[3]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+3]); // m=+1
    ratios_l2[4]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+4]); // m=+2 
 
    // l=3 Heights are defined by four parameters 
    ratios_l3[0]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+8]);  // m=-3
    ratios_l3[1]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+7]); // m=-2
    ratios_l3[2]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+6]); // m=-1 
    ratios_l3[3]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+5]);  // m=0
    ratios_l3[4]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+6]); // m=+1
    ratios_l3[5]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+7]); // m=+2 
    ratios_l3[6]=std::abs(params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+8]); // m=+3
 
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

    a1=std::abs(params[Nmax + lmax + Nf]);
    //eta=params[Nmax + lmax + Nf + 1];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    cpt=0;
    for(long n=0; n<Nmax; n++){
                
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);   
        if(do_amp){
            Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
            //std::cout << "[0] do conversion" << std::endl;
        } else{
            Hl0=std::abs(params[n]);
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }      
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
                //std::cout << "[1] do conversion" << std::endl;
            } else{
                Hl1=std::abs(params[n]*Vl1);
            }               
            model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }           
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            if(do_amp){
                Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
                //std::cout << "[2] do conversion" << std::endl;
            } else{
                Hl2=std::abs(params[n]*Vl2);
            }   
            model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            if(do_amp){
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
            } else{
                Hl3=std::abs(params[n]*Vl3);            
            }       
            model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            } 
        }     
    }

    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise backgound
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }

    return model_final;
}


// Added on 10 Feb 2020
// Alternative version  of the legacy function that is not dealing with the inclination
// but instead fit the mode heights directly. Visibilities are not considered here.
// WARNING: This model is not really suitable for a global fit due to the large set of parameters that might
//          Arise from letting all heights free for all m E ([l*(l+1)-1]/2 + 1) * Ntotal ~ 100 free parameters for a global fit with l=3
VectorXd model_MS_Global_a1etaa3_HarveyLike_Classic_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of l=1 frequencies
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Equal to the  number of heights for m components: Sum(l*(l+1))_l={1,2,3}
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];

    double trunc_c;

    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
    VectorXd Hl0(1), Hl1(3), Hl2(5), Hl3(7);
    double fl0, fl1, fl2, fl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;
    int Nharvey;
    long cpt, pos0;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    double inclination=-1;

    if (outparams){
        std::cout << " ERROR: outparams is not implemented for  " << __func__ << std::endl;
        std::cout << "       If you want to get outputs in ascii format, you must review this function and implement outparams = true" << std::endl;
        std::cout << "       The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    } 
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];

    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths
    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nmax);

    a1=std::abs(params[Nmax + lmax + Nf]);
    //eta=params[Nmax + lmax + Nf + 1];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 2];
    asym=params[Nmax+lmax + Nf + 5];
    
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    cpt=0;
    for(long n=0; n<Nmax; n++){
                
        fl0=fl0_all[n];
        Wl0=std::abs(Wl0_all[n]);   
        if(do_amp){
            Hl0[0]=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
            //std::cout << "[0] do conversion" << std::endl;
        } else{
            Hl0[0]=std::abs(params[n]);
        }       
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }       
        if(lmax >=1){
            fl1=params[Nmax+lmax+Nfl0+n];
            Wl1=std::abs(lin_interpol(fl0_all, Wl0_all, fl1));
            pos0=2*n;
            Hl1[0]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            Hl1[1]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0];
            Hl1[2]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            if(do_amp){
                //tmp=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0, 2);
                Hl1=Hl1/(pi*Wl1);
            } 
            Hl1=Hl1.cwiseAbs();
            model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
        }
        if(lmax >=2){
            fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
            Wl2=std::abs(lin_interpol(fl0_all, Wl0_all, fl2));
            pos0=3*n;
            Hl2[0]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            Hl2[1]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];     
            Hl2[2]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0];
            Hl2[3]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            Hl2[4]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            if(do_amp){
                //Hl2=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0, 3)/(pi*Wl2);
                Hl2=Hl2/(pi*Wl2);
           }
            Hl2=Hl2.cwiseAbs();
            model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
        }
        if(lmax >=3){
            fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
            Wl3=std::abs(lin_interpol(fl0_all, Wl0_all, fl3));
            pos0=4*n;
            Hl3[0]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+3];
            Hl3[1]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            Hl3[2]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];     
            Hl3[3]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0];
            Hl3[4]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+1];
            Hl3[5]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+2];
            Hl3[6]=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0+3];
            if(do_amp){
               //Hl3=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise + pos0, 4)/(pi*Wl3);
                Hl3=Hl3/(pi*Wl3);
            }
            Hl3=Hl3.cwiseAbs();
            model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            } 
        }  
 /*       std::cout << "[" << n << "]" << std::endl;
        std::cout << "fl0[" << n << "] = " << fl0 << std::endl;
        std::cout << "Wl0[" << n << "] = " << Wl0 << std::endl;
        std::cout << "Hl0[" << n << "] = " << Hl0 << std::endl;
         std::cout << "---------" << std::endl;
        std::cout << "fl1[" << n << "] = " << fl1 << std::endl;
        std::cout << "Wl1[" << n << "] = " << Wl1 << std::endl;
        std::cout << "Hl1[" << n << "] = " << Hl1 << std::endl;
        std::cout << "---------" << std::endl;
        std::cout << "fl2[" << n << "] = " << fl2 << std::endl;
        std::cout << "Wl2[" << n << "] = " << Wl2 << std::endl;
        std::cout << "Hl2[" << n << "] = " << Hl2 << std::endl;
        std::cout << "---------" << std::endl;
        if(lmax>=3){
            std::cout << "fl3[" << n << "] = " << fl3 << std::endl;
            std::cout << "Wl3[" << n << "] = " << Hl3 << std::endl;
            std::cout << "Hl3[" << n << "] = " << Hl3 << std::endl;
            std::cout << "---------" << std::endl;
        }
        std::cout << "a1 = " << Wl0 << std::endl;
        std::cout << "eta = " << Wl1 << std::endl;
        std::cout << "a3 = " << Wl2 << std::endl;

       for(long indd=0; indd<model_final.size(); indd++){
            if (std::isfinite(model_final[indd]) == false){
                std::cout << indd << " Not finite" << std::endl;
                std::cout << "      " << model_final[indd] << std::endl;
                exit(EXIT_SUCCESS);
            }
        }  
*/
    }
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    return model_final;
}
// ------------------------------
// ------------------------------
// ------------------------------

VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v1(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 * The v1 version has the particularity to not explicitly impose numax as a parameter. Instead, it is 
	 * calculated at each iteration by a weighted average of the frequencies, with the weight being the heights
	 * This model has therefore 5 parameters for the widths
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
	
	double inclination;

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;
	double numax, Htot, lnGamma0, lnLorentz;
	double e;
	
	int Nharvey;
	long cpt;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited

	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

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
	
	// ----- Widths ----
	numax=0.;
	Htot=0.;
	for(long n=0; n<Nmax; n++){
		numax=numax+params[n]*params[Nmax + lmax + n]; // Adding up the l=0: H(l=0, n)*nu(l=0, n)
		Htot=Htot+params[n];
		if(lmax>=1){
			numax=numax+params[n]*Vl1*params[Nmax+lmax+Nfl0+n]; // Adding up the l=1: H(l=1, n)*nu(l=1, n)
			Htot=Htot+params[n]*Vl1;
		}
		if(lmax>=2){
			numax=numax+params[n]*Vl2*params[Nmax+lmax+Nfl0+Nfl1+n]; // Adding up the l=2: H(l=2, n)*nu(l=2, n)
			Htot=Htot+params[n]*Vl2;
		}
		if(lmax>=3){
			numax=numax+params[n]*Vl3*params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]; // Adding up the l=3: H(l=3, n)*nu(l=3, n)
			Htot=Htot+params[n]*Vl3;
		}
	}
	numax=numax/Htot;
	
	std::cout << "numax :" << numax << std::endl;
	// -----------------
	//a1=std::abs(params[Nmax + lmax + Nf]);
	//eta=params[Nmax + lmax + Nf + 1];
	eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	cpt=0;
	for(long n=0; n<Nmax; n++){
		fl0=fl0_all[n];
		lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl0/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
		e=2.*log(fl0/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
		lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
		Wl0=exp(lnGamma0 + lnLorentz);

		//std::cout << "Wl0=" << Wl0 << std::endl;
		 
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		

		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   	
		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl1/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
			e=2.*log(fl1/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
			Wl1=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl1=" << Wl1 << std::endl;

			if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
           if (outparams){
                mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }       
		}
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl2/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
			e=2.*log(fl2/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
			Wl2=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl2=" << Wl2 << std::endl;

			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }    
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+1] * log(fl3/numax) + log(params[Nmax+lmax+Nf+Nsplit+2]);
			e=2.*log(fl3/params[Nmax+lmax+Nf+Nsplit+0]) / log(params[Nmax+lmax+Nf+Nsplit+3]/numax);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+4])/(1. + pow(e,2));		
			Wl3=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl3=" << Wl3 << std::endl;

			if(do_amp){
              Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            if (outparams){
                mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                Line=Line+1;
            }   
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
	return model_final;
}


VectorXd model_MS_Global_a1etaa3_AppWidth_HarveyLike_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
	/* Model of the power spectrum of a Main sequence solar-like star
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
	 * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
	 *          Size MUST be 0 otherwise (this check is not made in this function)
	 * The v2 version has the particularity to explicitly use numax as a parameter. 
	 * This model has therefore 6 parameters for the widths contrary to the 5 parameters of the v1 version
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	const int Nmax=params_length[0]; // Number of Heights
	const int lmax=params_length[1]; // number of degree - 1
	const int Nfl0=params_length[2]; // number of l=0 frequencies
	const int Nfl1=params_length[3]; // number of l=1 frequencies
	const int Nfl2=params_length[4]; // number of l=2 frequencies
	const int Nfl3=params_length[5]; // number of l=3 frequencies
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a ALL global MS model (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7 for a global MS model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
	
	double inclination;
	

	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd fl0_all(Nmax), Wl0_all(Nmax), noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;
	double Htot, lnGamma0, lnLorentz;
	double e;
	
	int Nharvey;
	long cpt;
	

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	inclination=std::atan(params[Nmax + lmax + Nf+4]/params[Nmax + lmax + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + lmax + Nf+3],2)+ pow(params[Nmax + lmax + Nf+4],2);

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

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
	
	// -----------------
	//a1=std::abs(params[Nmax + lmax + Nf]);
	//eta=params[Nmax + lmax + Nf + 1];
	eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 2];
	asym=params[Nmax+lmax + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
   
	cpt=0;
	for(long n=0; n<Nmax; n++){
		fl0=fl0_all[n];
		lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl0/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
		e=2.*log(fl0/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
		lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
		
		Wl0=exp(lnGamma0 + lnLorentz);

		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		

		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
               mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
               Line=Line+1;
        } 		
		if(lmax >=1){
			fl1=params[Nmax+lmax+Nfl0+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl1/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
			e=2.*log(fl1/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
			Wl1=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl1=" << Wl1 << std::endl;

			if(do_amp){
                Hl1=std::abs(params[n]/(pi*Wl1))*Vl1;
			} else{
				Hl1=std::abs(params[n]*Vl1);
			}				
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
		    //debug(model_final, Hl1, fl1, a1, eta, a3, asym, Wl1, 1, step, inclination, ratios_l1, trunc_c, true);
            if (outparams){
               mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
               Line=Line+1;
            }         
        }
		if(lmax >=2){
			fl2=params[Nmax+lmax+Nfl0+Nfl1+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl2/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
			e=2.*log(fl2/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
			Wl2=exp(lnGamma0 + lnLorentz);

			//std::cout << "Wl2=" << Wl2 << std::endl;

			if(do_amp){
				Hl2=std::abs(params[n]/(pi*Wl2))*Vl2;
			} else{
				Hl2=std::abs(params[n]*Vl2);
			}	
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
            //debug(model_final, Hl2, fl2, a1, eta, a3, asym, Wl2, 2, step, inclination, ratios_l2, trunc_c, true);
            if (outparams){
                   mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                   Line=Line+1;
            }             
		}
        if(lmax >=3){
			fl3=params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n];
			lnGamma0=params[Nmax + lmax + Nf + Nsplit+2] * log(fl3/params[Nmax + lmax + Nf + Nsplit + 0]) + log(params[Nmax+lmax+Nf+Nsplit+3]);
			e=2.*log(fl3/params[Nmax+lmax+Nf+Nsplit+1]) / log(params[Nmax+lmax+Nf+Nsplit+4]/params[Nmax + lmax +Nf + Nsplit + 0]);
			lnLorentz=-log(params[Nmax+lmax+Nf+Nsplit+5])/(1. + pow(e,2));		
			Wl3=exp(lnGamma0 + lnLorentz);

			if(do_amp){
                Hl3=std::abs(params[n]/(pi*Wl3))*Vl3;
			} else{
				Hl3=std::abs(params[n]*Vl3);			
			}		
			model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
            //debug(model_final, Hl3, fl3, a1, eta, a3, asym, Wl3, 3, step, inclination, ratios_l3, trunc_c, true);
            if (outparams){
                   mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
                   Line=Line+1;
            }        
		}		
	}
    //std::cin.ignore();

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
	//std::cout << "Before computing noise model" << std::endl;
	noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
		
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	
     if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6/ a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }

	return model_final;
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// //////////////////////////////   Models for local fit \\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

VectorXd model_MS_local_basic(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){

	/* Model for a local fit of the power spectrum of a solar-like star
	 * Assumptions of constant splitting over the fit window make it more suitable for a fit of a MS star
	 * but RGB might be handled also provided that the user focus on single mode fitting (ie, avoids fitting groups of different l)
	 * param is a vector of parameters
	 * param_length defines the structure of the parameters
	 * x is the frequency assumed to be in microHz
	 */

	const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
	const long double pi = 3.141592653589793238462643383279502884L;
	
	const int Nmax=params_length[0]; // Total Number of Heights
	const int Nvis=params_length[1]; // number of visibilities
	const int Nfl0=params_length[2]; // number of l=0 frequencies, heights and widths
	const int Nfl1=params_length[3]; // number of l=1 frequencies, heights and widths
	const int Nfl2=params_length[4]; // number of l=2 frequencies, heights and widths
	const int Nfl3=params_length[5]; // number of l=3 frequencies, heights and widths
	const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a MS_local (a1,eta,a3, magb, magalfa, asym)
	const int Nwidth=params_length[7]; // Total number of parameters for the widths. Should be the same as Nmax for a global MS model
	const int Nnoise=params_length[8]; // number of parameters for the noise. Might be only the white noise for a local model
	const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
	const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

	const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
	const double trunc_c=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc];
	const bool do_amp=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
	
	double inclination;
	
	VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
	VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

	VectorXd noise_params(Nnoise); //Hl0_all[Nmax],
	double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;

	int Nharvey;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
	/*
	   -------------------------------------------------------
	   ------- Gathering information about the modes ---------
	   -------------------------------------------------------
	*/
	inclination=std::atan(params[Nmax + Nvis + Nf+4]/params[Nmax + Nvis + Nf+3]); 
	inclination=inclination*180./pi;
	a1=pow(params[Nmax + Nvis + Nf+3],2)+ pow(params[Nmax + Nvis + Nf+4],2);

    //std::cout << "a1 = " << a1 << std::endl;
    //std::cout << "inclination = " << inclination << std::endl;

	// Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
	ratios_l0.setOnes();
	if(Nfl1 >=1){
//        Vl1=std::abs(params[Nmax]);
		ratios_l1=amplitude_ratio(1, inclination);
        //std::cout << "Vl1 " << Vl1 << std::endl;
    }
	if(Nfl2 >=1){
//		Vl2=std::abs(params[Nmax+1]);
		ratios_l2=amplitude_ratio(2, inclination);
        //std::cout << "Vl2 " << Vl2 << std::endl;

	}
	if(Nfl3 >=1){
//		Vl3=std::abs(params[Nmax+2]);
		ratios_l3=amplitude_ratio(3, inclination);
	}

	//a1=std::abs(params[Nmax + Nvis + Nf]);
	eta0=params[Nmax + Nvis + Nf + 1]; // 7 Dec 2021: Set by io_local to set it using Dnu
	//eta0=eta0_fct(fl0_all);
    a3=params[Nmax + Nvis + Nf + 2];
	asym=params[Nmax+Nvis + Nf + 5];
	
	model_final.setZero();
	
	/* -------------------------------------------------------
	   --------- Computing the models for the modes  ---------
	   -------------------------------------------------------
	*/
	for(long n=0; n<Nfl0; n++){		
		fl0=params[Nmax + Nvis + n];
		Wl0=std::abs(params[Nmax + Nvis + Nf + Nsplit + n ]);		
		if(do_amp){
			Hl0=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
		} else{
			Hl0=std::abs(params[n]);
		}		
		//std::cout << "fl0 = " << fl0 << "      Hl0 = " << Hl0 << "      Wl0 = " << Wl0 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
	}
	for(long n=0; n<Nfl1; n++){		
		fl1=params[Nmax + Nvis + Nfl0 + n];
		Wl1=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + n ]);	
		if(do_amp){
    		Hl1=std::abs(params[Nfl0 + n]/(pi*Wl1));
		} else{
			Hl1=std::abs(params[Nfl0 + n]);
		}				
		//std::cout << "fl1 = " << fl1 << "      Hl1 = " << Hl1 << "      Wl1 = " << Wl1 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
	}
	
	for(long n=0; n<Nfl2; n++){		
		fl2=params[Nmax + Nvis + Nfl0 + Nfl1 + n];
		Wl2=std::abs(params[Nmax+ Nvis + Nf + Nsplit + Nfl0 + Nfl1 + n ]);	
		if(do_amp){
			Hl2=std::abs(params[Nfl0 + Nfl1 + n]/(pi*Wl2));
		} else{
			Hl2=std::abs(params[Nfl0 + Nfl1 + n]);
		}	
		//std::cout << "fl2 = " << fl2 << "      Hl2 = " << Hl2 << "      Wl2 = " << Wl2 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
	}
	for(long n=0; n<Nfl3; n++){		
        //std::cout << "BEFORE l=3" << std::endl;
		fl3=params[Nmax + Nvis + Nfl0 + Nfl1 + Nfl2 + n];
        //std::cout << "fl3 = " << fl3 << std::endl;
		Wl3=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + Nfl1 + Nfl2 + n ]);
        //std::cout << "Wl3 = " << Wl3 << std::endl;
		if(do_amp){
        	Hl3=std::abs(params[Nfl0 + Nfl1 + Nfl2 + n]/(pi*Wl3));
		} else{
			Hl3=std::abs(params[Nfl0 + Nfl1 + Nfl2 + n]);			
		}		
		//std::cout << "fl3 = " << fl3 << "      Hl3 = " << Hl3 << "      Wl3 = " << Wl3 << std::endl;
		model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }  
	}

	/* -------------------------------------------------------
	   ------- Gathering information about the noise ---------
	   -------------------------------------------------------
	*/
    //std::cout << "Before noise_params" << std::endl;
	noise_params=params.segment(Nmax+Nvis+Nf+Nsplit+Nwidth, Nnoise);
	Nharvey=0; //(Nnoise-1)/3;
	//std::cout << "Nharvey = " << Nharvey << std::endl;
	/* -------------------------------------------------------
	   ---------- Computing the mode of the noise ------------
	   -------------------------------------------------------
	*/
	model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
	//std::cout << "End test" << std::endl;
    //exit(EXIT_SUCCESS);
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
	return model_final;
}


VectorXd model_MS_local_Hnlm(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){

    /* Model for a local fit of the power spectrum of a solar-like star
     * Assumptions of constant splitting over the fit window make it more suitable for a fit of a MS star
     * but RGB might be handled also provided that the user focus on single mode fitting (ie, avoids fitting groups of different l)
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * This model differs from model_MS_local_basic() by the fact that it fits the heights Hnlm instead than the inclination.
     * Therefore it compares to the model_MS_Global_a1etaa3_HarveyLike_Classic_v3() but for the local class of models
     * This means that the it preserve the symetry H(n,l,+m) = H(n,l, -m). It has to be used jointly with extra_priors[3]=2
     * Which imposes Sum(Hnlm)_{m=-l, m+l} = 1 in priors_local()
     */

    const double step=x[1]-x[0]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    
    const int Nmax=params_length[0]; // Total Number of Heights
    const int Nvis=params_length[1]; // number of visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies, heights and widths
    const int Nfl1=params_length[3]; // number of l=1 frequencies, heights and widths
    const int Nfl2=params_length[4]; // number of l=2 frequencies, heights and widths
    const int Nfl3=params_length[5]; // number of l=3 frequencies, heights and widths
    const int Nsplit=params_length[6]; // number of splitting parameters. Should be 6 for a MS_local (a1,eta,a3, magb, magalfa, asym)
    const int Nwidth=params_length[7]; // Total number of parameters for the widths. Should be the same as Nmax for a global MS model
    const int Nnoise=params_length[8]; // number of parameters for the noise. Might be only the white noise for a local model
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1 for a global MS model
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)

    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+Nvis+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    
    int pos0;
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd noise_params(Nnoise);
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Wl0, Wl1, Wl2, Wl3, a1,eta0,a3, asym;
    VectorXd Hl0(1), Hl1(3), Hl2(5), Hl3(7);
    int Nharvey;
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    a1=std::abs(params[Nmax + Nvis + Nf]);
    eta0=params[Nmax + Nvis + Nf + 1]; // 7 Dec 2021: Set in the io_local using Dnu
    //eta0=eta0_fct(fl0_all);
    a3=params[Nmax + Nvis + Nf + 2];
    asym=params[Nmax+Nvis + Nf + 5];
    
    model_final.setZero();

    if (outparams){
        std::cout << " ERROR: outparams is not implemented for  " << __func__ << std::endl;
        std::cout << "       If you want to get outputs in ascii format, you must review this function and implement outparams = true" << std::endl;
        std::cout << "       The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    } 
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */
    for(long n=0; n<Nfl0; n++){     
        fl0=params[Nmax + Nvis + n];
        Wl0=std::abs(params[Nmax + Nvis + Nf + Nsplit + n ]);       
        if(do_amp){
            Hl0[0]=std::abs(params[n]/(pi*Wl0)); // A^2/(pi.Gamma)
        } else{
            Hl0[0]=std::abs(params[n]);
        }       
        //std::cout << "fl0 = " << fl0 << "      Hl0 = " << Hl0 << "      Wl0 = " << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl0, fl0, a1, eta0, a3, asym, Wl0, 0, step, trunc_c);
    }
    for(long n=0; n<Nfl1; n++){     
        fl1=params[Nmax + Nvis + Nfl0 + n];
        Wl1=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + n ]);    
        pos0=2*n;
        Hl1[0]=params[Nfl0 + pos0+1];
        Hl1[1]=params[Nfl0 + pos0];
        Hl1[2]=params[Nfl0 + pos0+1];
        if(do_amp){
                Hl1=Hl1/(pi*Wl1);
        }
        Hl1=Hl1.cwiseAbs(); 
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl1, fl1, a1, eta0, a3,asym, Wl1, 1, step, trunc_c);
    }    
    for(long n=0; n<Nfl2; n++){     
        fl2=params[Nmax + Nvis + Nfl0 + Nfl1 + n];
        Wl2=std::abs(params[Nmax+ Nvis + Nf + Nsplit + Nfl0 + Nfl1 + n ]);  
        pos0=3*n;
        Hl2[0]=params[Nfl0 + Nfl1 + pos0+2];
        Hl2[1]=params[Nfl0 + Nfl1 + pos0+1];     
        Hl2[2]=params[Nfl0 + Nfl1 + pos0];
        Hl2[3]=params[Nfl0 + Nfl1 + pos0+1];
        Hl2[4]=params[Nfl0 + Nfl1 + pos0+2];
        if(do_amp){
            Hl2=Hl2/(pi*Wl2);
        }
        Hl2=Hl2.cwiseAbs();
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl2, fl2, a1, eta0, a3,asym, Wl2, 2, step, trunc_c);
    }
    for(long n=0; n<Nfl3; n++){     
        fl3=params[Nmax + Nvis + Nfl0 + Nfl1 + Nfl2 + n];
        Wl3=std::abs(params[Nmax + Nvis + Nf + Nsplit + Nfl0 + Nfl1 + Nfl2 + n ]);
        pos0=4*n;
        Hl3[0]=params[Nfl0 + Nfl1 + Nfl2 + pos0+3];
        Hl3[1]=params[Nfl0 + Nfl1 + Nfl2 + pos0+2];
        Hl3[2]=params[Nfl0 + Nfl1 + Nfl2 + pos0+1];     
        Hl3[3]=params[Nfl0 + Nfl1 + Nfl2 + pos0];
        Hl3[4]=params[Nfl0 + Nfl1 + Nfl2 + pos0+1];
        Hl3[5]=params[Nfl0 + Nfl1 + Nfl2 + pos0+2];
        Hl3[6]=params[Nfl0 + Nfl1 + Nfl2 + pos0+3];
        if(do_amp){
            Hl3=Hl3/(pi*Wl3);
       }
        Hl3=Hl3.cwiseAbs();
        model_final=optimum_lorentzian_calc_a1etaa3_v2(x, model_final, Hl3, fl3, a1, eta0, a3, asym, Wl3, 3, step, trunc_c);
    }

    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    noise_params=params.segment(Nmax+Nvis+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=0; //(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    
    //std::cout << "End test" << std::endl;
    //exit(EXIT_SUCCESS);

    return model_final;
}


////////////////////////////////////
// Models for asymptotic fitting ///
////////////////////////////////////


VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * This version is a pure implementation of the asymptotics, with the possibility of evaluating/accounting for model uncertainties
                using linear model for the inaccuracies for nu(zeta,1) and for nu_s(zeta,1): sigma(nu) = sigma_p + sig*zeta(nu). 
     * This would lead to uncessary complications and increase model computation time
     * This model has therefore 5 parameters for the widths. 
     */
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 5 here as it uses the Appourchaux profile
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const double sigma_limit=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];

    int i_dbg=0;
   
    VectorXi posOK;
    VectorXd gamma_params(6);
    gamma_params << std::abs(params[Nmax + lmax + Nf + Nsplit + 0]) , std::abs(params[Nmax+lmax+Nf+Nsplit+1]) , std::abs(params[Nmax + lmax + Nf + Nsplit+2]),
            std::abs(params[Nmax+lmax+Nf+Nsplit+3]) , std::abs(params[Nmax+lmax+Nf+Nsplit+4]) , std::abs(params[Nmax+lmax+Nf+Nsplit+5]); 
   
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    VectorXd xfit, rfit;//, fmin, fmax;
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym;
    double Dnu_p, epsilon, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    
    int Nharvey;
    long cpt;

    //outparams=true;
    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=std::abs(params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]); 
 
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

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 

    // --- Determining the large separation and the epsilon for the current set of free l=0 mode frequencies ---
    xfit=linspace(0, fl0_all.size()-1, fl0_all.size());
    rfit=linfit(xfit, fl0_all); // fit[0] is the slope ==> Dnu and fit[1] is the ordinate at origin ==> fit[1]/fit[0] = epsilon
    Dnu_p=rfit[0];
    epsilon=rfit[1]/rfit[0];
    epsilon=epsilon - floor(epsilon);
    const double fmin=fl0_all.minCoeff() - Dnu_p; // That will define the computation range for l=1 mixed modes
    const double fmax=fl0_all.maxCoeff() + Dnu_p;

    for (int n=0; n<Nmax;n++)
    {
        lnGamma0=gamma_params[2] * log(fl0_all[n]/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl0_all[n]/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl0_all[n]=exp(lnGamma0 + lnLorentz);
    }

   if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       

    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=std::abs(params[Nmax + lmax + Nfl0 + 1]);
    const double alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2]);
    const double q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
    const double sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]); // Evaluation of the inaccuracy in terms of Heights
    const double sigma_m_0_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]); // Evaluation of the inaccuracy in terms of frequencies: Ordinate at origin ==> sigma_p
    const double sigma_m_1_l1=std::abs(params[Nmax + lmax + Nfl0 + 6]); // Evaluation of the inaccuracy in terms of frequencies: slope ==> zeta=1 ==> sigma_g
    const double sigma_a1_0_l1=std::abs(params[Nmax + lmax + Nfl0 + 7]); // Evaluation of the inaccuracy in terms of splitting a1: Ordinate at origin ==> sigma_p
    const double sigma_a1_1_l1=std::abs(params[Nmax + lmax + Nfl0 + 8]); // Evaluation of the inaccuracy in terms of splitting a1: slope
    const double nmax=std::abs(params[Nmax + lmax + Nfl0 + 9]); // l=1 p modes: nmax ~ numax/Dnu + epsilon
    const double alpha_p=std::abs(params[Nmax + lmax + Nfl0 + 10]); // l=1 p modes: curvature 
    const double rot_env=std::abs(params[Nmax + lmax + Nf]);
    const double rot_core=std::abs(params[Nmax + lmax + Nf+1]); 
 
    VectorXd ksi_pg, h1_h0_ratio,f_interp, h_interp;
 
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen_m(seed);
    std::normal_distribution<double> distrib_fl1m(0.,1);
    std::normal_distribution<double> distrib_Hl1m(0.,sigma_H_l1);
    
    // -------- DEBUG OF l=1 mixed modes -----
    /*
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    std::cout << "sigma_p =" << sigma_p_l1 << std::endl;
    std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    std::cout << "step =" << step << std::endl;
    */
    // ---------------------------------------

    //const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    //const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    //const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    const Data_eigensols freqs_l1=solve_mm_asymptotic_O2p(Dnu_p, epsilon, 1, delta0l, alpha_p, nmax, DPl, alpha_g, q_star, 0, fmin, fmax, step, true, false);
    //fl1_all=freqs_l1.nu_m;

    // Remove modes which may lie beyond the data range
    posOK=where_in_range(freqs_l1.nu_m, x.minCoeff(), x.maxCoeff(), 0); // Remove frequencies out of the requested range
    fl1_all.resize(posOK.size());
    if (posOK[0] != -1){
        for (int k=0;k<posOK.size(); k++){
            fl1_all[k]=freqs_l1.nu_m[posOK[k]];
        }
    }

     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019)

    if (sigma_m_0_l1 != 0 || sigma_m_1_l1 !=0){
       for (int en=0; en<fl1_all.size(); en++)
        {
             r = distrib_fl1m(gen_m);
             r=r * std::abs(sigma_m_0_l1 + sigma_m_1_l1*ksi_pg[en]);
             while (r > sigma_limit){
                r = distrib_fl1m(gen_m);
                r=r * std::abs(sigma_m_0_l1 + sigma_m_1_l1*ksi_pg[en]);
            }
            fl1_all[en]=fl1_all[en] + r;
        }
    } 

   
    // We will add 0-anchors at the edges to avoid negative extrapolation... so we need at least two extra points... we choose to have 4 extra points
    f_interp.resize(fl0_all.size()+4); 
    f_interp[0]=fl0_all.minCoeff()*0.6; // A failsafe beyond the actual range of frequencies
    f_interp[1]=fl0_all.minCoeff()*0.8; 
    f_interp[f_interp.size()-2]=fl0_all.maxCoeff()*1.2;
    f_interp[f_interp.size()-1]=fl0_all.maxCoeff()*1.4;

    h_interp.resize(Hl0_all.size()+4); 
    h_interp[0]=0;
    h_interp[1]=Hl0_all[0]/4;
    h_interp[h_interp.size()-2]=Hl0_all[Hl0_all.size()-1]/4;
    h_interp[h_interp.size()-1]=0;
    for(int j=0; j<fl0_all.size(); j++){
        f_interp[j+2]=fl0_all[j];
        h_interp[j+2]=Hl0_all[j];
    }
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {     
        tmp=lin_interpol(f_interp, h_interp, fl1_all[i]); // interpolate Hl0 to fl1 positions
        if (tmp < 0){ 
            std::cout << "WARNING: WE IMPOSE Hl1p_all = 0 due to Negative interpolation" << std::endl;
            Hl1p_all[i]=0;
        } else{
            Hl1p_all[i]=std::abs(tmp);

        }
    }

    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    if (sigma_H_l1 != 0){
        for (int i=0; i<Hl1_all.size();i++){
             r = distrib_fl1m(gen_m);
             while ((1+r) < 0){ // Avoid negative heights
                r = distrib_fl1m(gen_m);
             }
             Hl1_all[i]=Hl1_all[i]*(1 + r); // Perturbation proportional to the actual value
        }
    }
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1); // generate the mixed modes widths
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
    //eta=params[Nmax + lmax + Nf + 2];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
 
     if (sigma_a1_0_l1 != 0 || sigma_a1_1_l1 !=0){
       for (int en=0; en<fl1_all.size(); en++)
        {
            r = distrib_fl1m(gen_m);
            r=r * std::abs(sigma_a1_0_l1 + sigma_a1_1_l1*ksi_pg[en]);
            a1_l1[en]=std::abs(a1_l1[en] + r);
        }
    } 
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        //std::cout << "Wl0=" << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        //crash=debug(model_final, Hl0, fl0, 0, eta, a3, asym, Wl0, 0, step, inclination, ratios_l0, trunc_c, true);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    for (long n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        //crash=debug(model_final, Hl1, fl1, a1_l1[n], eta, a3, asym, Wl1, 1, step, inclination, ratios_l1, trunc_c, false);
        //if (crash == true){
        //    std::cout << " Full parameter list: " << params.transpose() << std::endl;
        //    std::cout << " Full parameters_length: " << params_length.transpose() << std::endl;
        //    exit(EXIT_FAILURE);
        //}
       if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]);
        lnGamma0=gamma_params[2] * log(fl2/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl2/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl2=exp(lnGamma0 + lnLorentz);
        //std::cout << "Wl2=" << Wl2 << std::endl;
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        //crash=debug(model_final, Hl2, fl2, a1_l2[n], eta, a3, asym, Wl2, 2, step, inclination, ratios_l2, trunc_c, true);
       if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l2[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    for(long n=0; n<Nfl3; n++){ 
        fl3=std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]);
        lnGamma0=gamma_params[2] * log(fl3/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl3/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl3=exp(lnGamma0 + lnLorentz);
        if(do_amp){
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3/(pi*Wl3)*Vl3);
        } else{
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        //crash=debug(model_final, Hl3, fl3, a1_l3[n], eta, a3, asym, Wl3, 3, step, inclination, ratios_l3, trunc_c, true);
       if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l3[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }           
//    std::cout << "--------------- " << std::endl;
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    //debug(model_final, -1, -1, -1, -1, -1, -1, -1, 0, step, -1, ratios_l0, trunc_c, true);
 
     /*std::cout << "fl0    / Wl0     / Hl0" << std::endl;
    for (int i =0; i<fl0_all.size(); i++){
        std::cout << fl0_all[i]  << "    "  << Wl0_all[i]  << "    " << Hl0_all[i] << std::endl;
    }
    std::cout << "fl1    / Wl1     / Hl1    / a1_l1" << std::endl;
    for (int i =0; i<fl1_all.size(); i++){
        std::cout << fl1_all[i]  << "    "  << Wl1_all[i]  << "    " << Hl1_all[i]  << "    " << a1_l1[i]  << std::endl;
    }
    std::cout << "a1_l2 = " << a1_l2.transpose() << std::endl;
    
    std::cout << "a1_l3 = " << a1_l3.transpose() << std::endl;

    std::cout << "End test" << std::endl;
    exit(EXIT_SUCCESS);
    */
    //  ---- DEBUG ----
    /*    std::cout << " DEBUG: Writing Results on-screen... " << std::endl;
        std::cout << "     - Pure p modes:" << std::endl;
        std::cout << freqs_l1.nu_p.transpose() << std::endl;
        std::cout << "     - Derivative of pure p modes dnup/dn:" << std::endl;
        std::cout << freqs_l1.dnup.transpose() << std::endl;
        std::cout << "     - Pure g modes:" << std::endl;
        std::cout << freqs_l1.nu_g.transpose() << std::endl;
        std::cout << "     - Derivative of pure g modes dPg/dn:" << std::endl;
        std::cout << freqs_l1.dPg.transpose() << std::endl;
        std::cout << "     - Mixed modes:" << std::endl;
        std::cout << freqs_l1.nu_m.transpose() << std::endl;
        std::cout << "     - Mixed modes (filtered):" << std::endl;
        std::cout << fl1_all.transpose() << std::endl;
        std::cout << "     - ksi_pg_filtered with original values After Filtering by max(fl0_all) and max(fl0_all):" << std::endl;
        std::cout << ksi_pg.transpose() << std::endl;
        std::cout << "     - h1_h0 with original values After Filtering by max(fl0_all) and max(fl0_all):" << std::endl;
        std::cout << h1_h0_ratio.transpose() << std::endl;
        std::cout << "     - Hl1p_all:" << std::endl;
        std::cout << Hl1p_all.transpose() << std::endl;
        std::cout << "     - Hl1_all:" << std::endl;
        std::cout << Hl1_all.transpose() << std::endl;
        std::cout << "     - Wl1_all:" << std::endl;
        std::cout << Wl1_all.transpose() << std::endl;
    */
    // -------
   
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }

    //exit(EXIT_SUCCESS);
    return model_final;
}



VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v2(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * The v2 version has the particularity to explicitly impose numax as a parameter. This is because we cannot create a v1
     * version unless we include so iterative process to evaluate numax from a weighted average of the frequencies. 
     * This would lead to uncessary complications and increase model computation time
     * This model has therefore 5 parameters for the widths. 
     * NOTE THAT IN THE RGB_asympt model, the v2 version is not implemented
     */
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 5 here as it uses the Appourchaux profile
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const double sigma_limit=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];

    int i_dbg=0;
   
    VectorXi posOK;
    VectorXd gamma_params(6);
    gamma_params << std::abs(params[Nmax + lmax + Nf + Nsplit + 0]) , std::abs(params[Nmax+lmax+Nf+Nsplit+1]) , std::abs(params[Nmax + lmax + Nf + Nsplit+2]),
            std::abs(params[Nmax+lmax+Nf+Nsplit+3]) , std::abs(params[Nmax+lmax+Nf+Nsplit+4]) , std::abs(params[Nmax+lmax+Nf+Nsplit+5]); 
   
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    
    int Nharvey;
    long cpt;

    //outparams=true;
    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=std::abs(params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]); 
 
    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);

    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 

    for (int n=0; n<Nmax;n++)
    {
        lnGamma0=gamma_params[2] * log(fl0_all[n]/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl0_all[n]/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl0_all[n]=exp(lnGamma0 + lnLorentz);
    }

   if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       

    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
  
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=std::abs(params[Nmax + lmax + Nfl0 + 1]);
    const double alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2]);
    const double q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
    const double sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]);
    const double sigma_g_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]);
    const double sigma_m_l1=std::abs(params[Nmax + lmax + Nfl0 + 6]);
    const double rot_env=std::abs(params[Nmax + lmax + Nf]);
    const double rot_core=std::abs(params[Nmax + lmax + Nf+1]);
 
    VectorXd ksi_pg, h1_h0_ratio,f_interp, h_interp;
 
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen_m(seed);    
    std::normal_distribution<double> distrib_fl1m(0.,sigma_m_l1);
    std::normal_distribution<double> distrib_Hl1m(0.,sigma_H_l1);
    
    // -------- DEBUG OF l=1 mixed modes -----
    /*
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    std::cout << "sigma_p =" << sigma_p_l1 << std::endl;
    std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    std::cout << "step =" << step << std::endl;
    */
    // ---------------------------------------

    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    //fl1_all=freqs_l1.nu_m;

    // Remove modes which may lie beyond the data range
    posOK=where_in_range(freqs_l1.nu_m, x.minCoeff(), x.maxCoeff(), 0); // Remove frequencies out of the requested range
    fl1_all.resize(posOK.size());
    if (posOK[0] != -1){
        for (int k=0;k<posOK.size(); k++){
            fl1_all[k]=freqs_l1.nu_m[posOK[k]];
        }
    }

    if (sigma_m_l1 != 0){
       for (int en=0; en<fl1_all.size(); en++)
        {
             r = distrib_fl1m(gen_m);
             while (r > sigma_limit){
                r = distrib_fl1m(gen_m);
            }
            fl1_all[en]=fl1_all[en] + r;
        }
    } 

     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019)
   
    // We will add 0-anchors at the edges to avoid negative extrapolation... so we need at least two extra points... we choose to have 4 extra points
    f_interp.resize(fl0_all.size()+4); 
    f_interp[0]=fl0_all.minCoeff()*0.6; // A failsafe beyond the actual range of frequencies
    f_interp[1]=fl0_all.minCoeff()*0.8; 
    f_interp[f_interp.size()-2]=fl0_all.maxCoeff()*1.2;
    f_interp[f_interp.size()-1]=fl0_all.maxCoeff()*1.4;

    h_interp.resize(Hl0_all.size()+4); 
    h_interp[0]=0;
    h_interp[1]=Hl0_all[0]/4;
    h_interp[h_interp.size()-2]=Hl0_all[Hl0_all.size()-1]/4;
    h_interp[h_interp.size()-1]=0;
    for(int j=0; j<fl0_all.size(); j++){
        f_interp[j+2]=fl0_all[j];
        h_interp[j+2]=Hl0_all[j];
    }
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {     
        tmp=lin_interpol(f_interp, h_interp, fl1_all[i]); // interpolate Hl0 to fl1 positions
        if (tmp < 0){ 
            std::cout << "WARNING: WE IMPOSE Hl1p_all = 0 due to Negative interpolation" << std::endl;
            Hl1p_all[i]=0;
        } else{
            Hl1p_all[i]=std::abs(tmp);

        }
    }

    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    if (sigma_H_l1 != 0){
        for (int i=0; i<Hl1_all.size();i++){
             r = distrib_fl1m(gen_m);
             while ((1+r) < 0){ // Avoid negative heights
                r = distrib_fl1m(gen_m);
             }
             tmp=tmp*(1 + r); // Perturbation proportional to the actual value
        }
    }
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1); // generate the mixed modes widths
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
    //eta=params[Nmax + lmax + Nf + 2];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
 
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        //std::cout << "Wl0=" << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        //crash=debug(model_final, Hl0, fl0, 0, eta, a3, asym, Wl0, 0, step, inclination, ratios_l0, trunc_c, true);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    for (long n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        //crash=debug(model_final, Hl1, fl1, a1_l1[n], eta, a3, asym, Wl1, 1, step, inclination, ratios_l1, trunc_c, false);
        //if (crash == true){
        //    std::cout << " Full parameter list: " << params.transpose() << std::endl;
        //    std::cout << " Full parameters_length: " << params_length.transpose() << std::endl;
        //    exit(EXIT_FAILURE);
        //}
       if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]);
        lnGamma0=gamma_params[2] * log(fl2/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl2/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl2=exp(lnGamma0 + lnLorentz);
        //std::cout << "Wl2=" << Wl2 << std::endl;
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        //crash=debug(model_final, Hl2, fl2, a1_l2[n], eta, a3, asym, Wl2, 2, step, inclination, ratios_l2, trunc_c, true);
       if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l2[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    for(long n=0; n<Nfl3; n++){ 
        fl3=std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]);
        lnGamma0=gamma_params[2] * log(fl3/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl3/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl3=exp(lnGamma0 + lnLorentz);
        if(do_amp){
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3/(pi*Wl3)*Vl3);
        } else{
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        //crash=debug(model_final, Hl3, fl3, a1_l3[n], eta, a3, asym, Wl3, 3, step, inclination, ratios_l3, trunc_c, true);
       if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l3[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }           
//    std::cout << "--------------- " << std::endl;
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    //debug(model_final, -1, -1, -1, -1, -1, -1, -1, 0, step, -1, ratios_l0, trunc_c, true);
 
     /*std::cout << "fl0    / Wl0     / Hl0" << std::endl;
    for (int i =0; i<fl0_all.size(); i++){
        std::cout << fl0_all[i]  << "    "  << Wl0_all[i]  << "    " << Hl0_all[i] << std::endl;
    }
    std::cout << "fl1    / Wl1     / Hl1    / a1_l1" << std::endl;
    for (int i =0; i<fl1_all.size(); i++){
        std::cout << fl1_all[i]  << "    "  << Wl1_all[i]  << "    " << Hl1_all[i]  << "    " << a1_l1[i]  << std::endl;
    }
    std::cout << "a1_l2 = " << a1_l2.transpose() << std::endl;
    
    std::cout << "a1_l3 = " << a1_l3.transpose() << std::endl;

    std::cout << "End test" << std::endl;
    exit(EXIT_SUCCESS);
    */
    //  ---- DEBUG ----
    /*    std::cout << " DEBUG: Writing Results on-screen... " << std::endl;
        std::cout << "     - Pure p modes:" << std::endl;
        std::cout << freqs_l1.nu_p.transpose() << std::endl;
        std::cout << "     - Derivative of pure p modes dnup/dn:" << std::endl;
        std::cout << freqs_l1.dnup.transpose() << std::endl;
        std::cout << "     - Pure g modes:" << std::endl;
        std::cout << freqs_l1.nu_g.transpose() << std::endl;
        std::cout << "     - Derivative of pure g modes dPg/dn:" << std::endl;
        std::cout << freqs_l1.dPg.transpose() << std::endl;
        std::cout << "     - Mixed modes:" << std::endl;
        std::cout << freqs_l1.nu_m.transpose() << std::endl;
        std::cout << "     - Mixed modes (filtered):" << std::endl;
        std::cout << fl1_all.transpose() << std::endl;
        std::cout << "     - ksi_pg_filtered with original values After Filtering by max(fl0_all) and max(fl0_all):" << std::endl;
        std::cout << ksi_pg.transpose() << std::endl;
        std::cout << "     - h1_h0 with original values After Filtering by max(fl0_all) and max(fl0_all):" << std::endl;
        std::cout << h1_h0_ratio.transpose() << std::endl;
        std::cout << "     - Hl1p_all:" << std::endl;
        std::cout << Hl1p_all.transpose() << std::endl;
        std::cout << "     - Hl1_all:" << std::endl;
        std::cout << Hl1_all.transpose() << std::endl;
        std::cout << "     - Wl1_all:" << std::endl;
        std::cout << Wl1_all.transpose() << std::endl;
    */
    // -------
   
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    return model_final;
}


VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * The v3 version is the same as the v2 version (explicitly impose numax as a parameter). 
     * But it requires the l=1 p modes as input to derive the mixed modes
     * This would lead to uncessary complications and increase model computation time
     * This model has therefore 5 parameters for the widths. 
     * NOTE THAT IN THE RGB_asympt model, the v2 version is not implemented
     */
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m + all the fl1p modes
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 5 here as it uses the Appourchaux profile
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const double sigma_limit=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];

    int i_dbg=0;
   
    VectorXd gamma_params(6);
    gamma_params << std::abs(params[Nmax + lmax + Nf + Nsplit + 0]) , std::abs(params[Nmax+lmax+Nf+Nsplit+1]) , std::abs(params[Nmax + lmax + Nf + Nsplit+2]),
            std::abs(params[Nmax+lmax+Nf+Nsplit+3]) , std::abs(params[Nmax+lmax+Nf+Nsplit+4]) , std::abs(params[Nmax+lmax+Nf+Nsplit+5]); 
   
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=std::abs(params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]); 
 
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

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 

    for (int n=0; n<Nmax;n++)
    {
        lnGamma0=gamma_params[2] * log(fl0_all[n]/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl0_all[n]/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl0_all[n]=exp(lnGamma0 + lnLorentz);
    }

   if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       

    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
  
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=std::abs(params[Nmax + lmax + Nfl0 + 1]);
    const double alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2]);
    const double q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
    const double sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]);
    const double sigma_g_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]);
    const double sigma_m_l1=std::abs(params[Nmax + lmax + Nfl0 + 6]);
    const double Hfactor=std::abs(params[Nmax + lmax + Nfl0 + 7]);
    const VectorXd fl1p_all=params.segment(Nmax+ lmax + Nfl0 + 8, Nfl1 - 8); // The total number fl1p modes is the total number of params Nfl1 - 8
    const double rot_env=std::abs(params[Nmax + lmax + Nf]);
    const double rot_core=std::abs(params[Nmax + lmax + Nf+1]);
    VectorXd ksi_pg, h1_h0_ratio,f_interp, h_interp;
 
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen_m(seed);    
    std::normal_distribution<double> distrib_fl1m(0.,sigma_m_l1);
    std::normal_distribution<double> distrib_Hl1m(0.,sigma_H_l1);
    
    // -------- DEBUG OF l=1 mixed modes -----
    /*
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    std::cout << "sigma_p =" << sigma_p_l1 << std::endl;
    std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    std::cout << "step =" << step << std::endl;
    */
    // ---------------------------------------

    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_nupl(fl1p_all, 1,  DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    fl1_all=freqs_l1.nu_m;

    if (sigma_m_l1 != 0){
       for (int en=0; en<fl1_all.size(); en++)
        {
             r = distrib_fl1m(gen_m);
             while (r > sigma_limit){
                r = distrib_fl1m(gen_m);
            }
            fl1_all[en]=fl1_all[en] + r;
        }
    } 

     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg, Hfactor); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019).
   
    // We will add 0-anchors at the edges to avoid negative extrapolation... so we need at least two extra points... we choose to have 4 extra points
    f_interp.resize(fl0_all.size()+4); 
    f_interp[0]=fl0_all.minCoeff()*0.6; // A failsafe beyond the actual range of frequencies
    f_interp[1]=fl0_all.minCoeff()*0.8; 
    f_interp[f_interp.size()-2]=fl0_all.maxCoeff()*1.2;
    f_interp[f_interp.size()-1]=fl0_all.maxCoeff()*1.4;

    h_interp.resize(Hl0_all.size()+4); 
    h_interp[0]=0;
    h_interp[1]=Hl0_all[0]/4;
    h_interp[h_interp.size()-2]=Hl0_all[Hl0_all.size()-1]/4;
    h_interp[h_interp.size()-1]=0;
    for(int j=0; j<fl0_all.size(); j++){
        f_interp[j+2]=fl0_all[j];
        h_interp[j+2]=Hl0_all[j];
    }
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {     
        tmp=lin_interpol(f_interp, h_interp, fl1_all[i]); // interpolate Hl0 to fl1 positions
        if (tmp < 0){ 
            std::cout << "WARNING: WE IMPOSE Hl1p_all = 0 due to Negative interpolation" << std::endl;
            Hl1p_all[i]=0;
        } else{
            Hl1p_all[i]=std::abs(tmp);

        }
    }

    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    if (sigma_H_l1 != 0){
        for (int i=0; i<Hl1_all.size();i++){
             r = distrib_fl1m(gen_m);
             while ((1+r) < 0){ // Avoid negative heights
                r = distrib_fl1m(gen_m);
             }
             tmp=tmp*(1 + r); // Perturbation proportional to the actual value
        }
    }
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1); // generate the mixed modes widths
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
    //eta=params[Nmax + lmax + Nf + 2];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
   
    //bool error, crash;
    //error=debug_solver(x, fl1_all, fl0_all, 1, delta0l, DPl, alpha_g, q_star, sigma_p_l1);
    //if (error == true){
    //    std::cout << " Full parameter list: " << params.transpose() << std::endl;
    //    std::cout << " Full parameters_length: " << params_length.transpose() << std::endl;
    //    crash=debug(model_final, Hl1, fl1, -1, eta, a3, asym, Wl1, 1, step, inclination, ratios_l1, trunc_c, true);
    //}

    // --------------
    // --------------
    // --------------
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        //std::cout << "Wl0=" << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        //crash=debug(model_final, Hl0, fl0, 0, eta, a3, asym, Wl0, 0, step, inclination, ratios_l0, trunc_c, true);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    for (long n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
        //crash=debug(model_final, Hl1, fl1, a1_l1[n], eta, a3, asym, Wl1, 1, step, inclination, ratios_l1, trunc_c, false);
        //if (crash == true){
        //    std::cout << " Full parameter list: " << params.transpose() << std::endl;
        //    std::cout << " Full parameters_length: " << params_length.transpose() << std::endl;
        //    exit(EXIT_FAILURE);
        //}
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]);
        lnGamma0=gamma_params[2] * log(fl2/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl2/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl2=exp(lnGamma0 + lnLorentz);
        //std::cout << "Wl2=" << Wl2 << std::endl;
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        //crash=debug(model_final, Hl2, fl2, a1_l2[n], eta, a3, asym, Wl2, 2, step, inclination, ratios_l2, trunc_c, true);
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l2[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }     
    }

//    std::cout << "--------------- " << std::endl;
    for(long n=0; n<Nfl3; n++){ 
        fl3=std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]);
        lnGamma0=gamma_params[2] * log(fl3/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl3/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl3=exp(lnGamma0 + lnLorentz);
        if(do_amp){
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3/(pi*Wl3)*Vl3);
        } else{
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l3[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }  
        //crash=debug(model_final, Hl3, fl3, a1_l3[n], eta, a3, asym, Wl3, 3, step, inclination, ratios_l3, trunc_c, true);
    }           
//    std::cout << "--------------- " << std::endl;
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    //debug(model_final, -1, -1, -1, -1, -1, -1, -1, 0, step, -1, ratios_l0, trunc_c, true);
 
     /*std::cout << "fl0    / Wl0     / Hl0" << std::endl;
    for (int i =0; i<fl0_all.size(); i++){
        std::cout << fl0_all[i]  << "    "  << Wl0_all[i]  << "    " << Hl0_all[i] << std::endl;
    }
    std::cout << "fl1    / Wl1     / Hl1    / a1_l1" << std::endl;
    for (int i =0; i<fl1_all.size(); i++){
        std::cout << fl1_all[i]  << "    "  << Wl1_all[i]  << "    " << Hl1_all[i]  << "    " << a1_l1[i]  << std::endl;
    }
    std::cout << "a1_l2 = " << a1_l2.transpose() << std::endl;
    
    std::cout << "a1_l3 = " << a1_l3.transpose() << std::endl;

    std::cout << "End test" << std::endl;
    exit(EXIT_SUCCESS);
    */
    //  ---- DEBUG ----
    /*    std::cout << " DEBUG: Writing Results on-screen... " << std::endl;
        std::cout << "     - Pure p modes:" << std::endl;
        std::cout << freqs_l1.nu_p.transpose() << std::endl;
        std::cout << "     - Derivative of pure p modes dnup/dn:" << std::endl;
        std::cout << freqs_l1.dnup.transpose() << std::endl;
        std::cout << "     - Pure g modes:" << std::endl;
        std::cout << freqs_l1.nu_g.transpose() << std::endl;
        std::cout << "     - Derivative of pure g modes dPg/dn:" << std::endl;
        std::cout << freqs_l1.dPg.transpose() << std::endl;
        std::cout << "     - Mixed modes:" << std::endl;
        std::cout << freqs_l1.nu_m.transpose() << std::endl;
        std::cout << "     - Mixed modes (filtered):" << std::endl;
        std::cout << fl1_all.transpose() << std::endl;
        std::cout << "     - ksi_pg_filtered with original values After Filtering by max(fl0_all) and max(fl0_all):" << std::endl;
        std::cout << ksi_pg.transpose() << std::endl;
        std::cout << "     - h1_h0 with original values After Filtering by max(fl0_all) and max(fl0_all):" << std::endl;
        std::cout << h1_h0_ratio.transpose() << std::endl;
        std::cout << "     - Hl1p_all:" << std::endl;
        std::cout << Hl1p_all.transpose() << std::endl;
        std::cout << "     - Hl1_all:" << std::endl;
        std::cout << Hl1_all.transpose() << std::endl;
        std::cout << "     - Wl1_all:" << std::endl;
        std::cout << Wl1_all.transpose() << std::endl;
    */
    // -------
    //exit(EXIT_SUCCESS);
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    return model_final;
}

VectorXd model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v4(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Widths are constant: It is usefull only when fitting over a few Dnu
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * Like v2 and v3, the v4 version explicitly imposes numax as a parameter. 
     * Contrary to v3, it has the following capabilities:
     *      - model_type = 1 : The user can choose to use the l=0, shifted with d01 in order to generate l=1 p modes.
     *             In this case, we make use of:
             solve_mm_asymptotic_O2from_l0(nu_l0_in, el, delta0l, DPl, alpha_g, q, resol, bool returns_pg_freqs=true, verbose=false, freq_min=fmin, freq_max=fmax)
     *      - OR model_type=0: it can use an O2 polynomial:   
               solve_mm_asymptotic_O2p(Dnu_p, epsilon, el, delta0l, alpha_p, nmax, DPl, alpha_g, q, fmin, fmax, resol, true, false); //returns_pg_freqs=true, verbose=false
              Here, Dnu_p, epsilon are calculated using a linear fit of l=0 frequencies.
              alpha_p, nmax are determined using a 2nd order fit of l=0 frequencies.
     *      The switch between these two option is made using a control parameter at the end of vector of parameters (replacing sigma_limit of v2 and v3 and named model_type)
    *   The fit incorporate a list of Nerr parameters situated between ferr=[Nmax+ lmax + Nfl0 + 8] and [Nmax+ lmax + Nfl0 + 8 +Nerr] that are used to correct from the bias of the exact asymptotic model
    *   These are associated to a list of Nerr list of fixed frequencies between fref=[Nmax+ lmax + Nfl0 + 8 + Nerr] and [Nmax+ lmax + Nfl0 + 8 + 2 Nerr] that provided frequencies
    *   of (a) case of bias_type = 0 (cte):
               The reference frequencies fref are used as 'windows' for determining region of 'constant bias': If a mode lies within a window [fref(j),fref(j+1)] of the list, then we apply the correction ferr(j)
               This means that multiple modes may have the same correction factor ferr(j) and that some window may have 0. Be aware of this when post processing: Each window must be cleaned of solutions with 0 modes
           (b) case of bias_type = 1 (Cubic spline, ie strongly correlated as it is twice continuously differentiable):
               The reference frequencies fref are used as 'anchors' for a spline defined by the values ferr.
               The spline fitting is performed here in this function using the spline module from https://github.com/ttk592/spline/ also forked to my repo
           (c) case of bias_type = 2 (Hermite spline, ie less correlated variable as it is once continuously differentiable):
               The reference frequencies fref are used as 'anchors' for a spline defined by the values ferr.
               The spline fitting is performed here in this function using the spline module from https://github.com/ttk592/spline/ also forked to my repo
     * This model has therefore 6 parameters for the widths. 
     */
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = M_PI;//3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m + all the fl1p modes
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 1
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const double sigma_limit=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];
    const double model_type=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+3];
    const double bias_type=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+4];
    const int Nferr=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+5];
    int i_dbg=0;
        
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    Data_eigensols freqs_l1;
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    outparams=true;
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=std::abs(params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]); 
 
    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);

    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 
    Wl0_all.setConstant(std::abs(params[Nmax + lmax + Nf + Nsplit + 0])); 
   
    if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       
    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=std::abs(params[Nmax + lmax + Nfl0 + 1]);
    const double alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2]);
    const double q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
    const double sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]);
    const double sigma_g_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]);
    const double Wfactor=std::abs(params[Nmax + lmax + Nfl0 + 6]);
    const double Hfactor=std::abs(params[Nmax + lmax + Nfl0 + 7]);
    std::vector<double> ferr_all(Nferr); // The total number fl1p modes is the total number of params Nfl1 - 8 // THESE MUST BE VARIABLE
    std::vector<double> fref_all(Nferr); // The total number fl1p modes is the total number of params Nfl1 - 8 // THESE MUST BE CONSTANT
    const double rot_env=std::abs(params[Nmax + lmax + Nf]);
    const double rot_core=std::abs(params[Nmax + lmax + Nf+1]);
    VectorXd ksi_pg, h1_h0_ratio,f_interp, h_interp;
 
    for (int i=0;i<Nferr; i++){
        fref_all[i]=params[Nmax+ lmax + Nfl0 + 8 + i];        
        ferr_all[i]=params[Nmax+ lmax + Nfl0 + 8 + Nferr + i];
    }
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //std::default_random_engine gen_m(seed);    
    //std::normal_distribution<double> distrib_fl1m(0.,sigma_m_l1);
    //std::normal_distribution<double> distrib_Hl1m(0.,sigma_H_l1);
    
    // bias type setup
    tk::spline bias; // declare a bias
    //tk::spline bias(fref_all, ferr_all);
    if (bias_type == 0){ // WINDOW CASE
        std::cout << " Windows case not implemented... Use cubic or hermite spline instead" << std::endl;
        exit(EXIT_SUCCESS);
    } else{
        bias.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0); // imposes second derivative to 0 at the edge (extrapolation if happens will be linear)
        if (bias_type == 1){ // CUBIC SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline);     // this calculates all spline coefficients   
        }
        if (bias_type == 2){ // HERMITE SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline_hermite);     // this calculates all spline coefficients
        }
    }
       // -------- DEBUG OF l=1 mixed modes -----
    /*
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    //std::cout << "sigma_p =" << sigma_p_l1 << std::endl;
    //std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    //std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "Hfactor =" << Hfactor << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    std::cout << "step =" << step << std::endl;
    */
    // ---------------------------------------

    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    //const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_nupl(fl1p_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    if (model_type == 0){
        VectorXd xfit, rfit;
        xfit=linspace(0, fl0_all.size()-1, fl0_all.size());
        rfit=linfit(xfit, fl0_all); // rfit[0] = Dnu
        const double Dnu_p=rfit[0];
        const long double n0=floor(rfit[0]/Dnu_p);
        const double epsilon_p = rfit[0]/Dnu_p -n0;
        //std::cout << " ---- 1st Order Fit ----" << std::endl;
        //std::cout << " Dnu_p =" << rfit[0] << std::endl;
        //std::cout << " epsilon_p + n0 =" << rfit[1]/rfit[0] << std::endl;
        //std::cout << " ---- 2nd Order Fit ----" << std::endl;
        //rfit=polyfitXd(xfit, fl0_all, 2); // rfit[0] = A, rfit[1]=B, rfit[2]=C : alpha_p=2 C / Dnu_p   ; nmax = (1 - B/Dnu_p)/alpha_p ; epsilon = A/Dnu_p - alpha_p * nmax^2 / 2
        //const long double alpha_p=2 *rfit[2]/Dnu_p;
        //const long double nmax = (1-rfit[1]/Dnu_p)/alpha_p;
        //const long double n0=floor(rfit[0]/Dnu_p - alpha_p * pow(nmax,2) / 2);
        //const long double epsilon_p = rfit[0]/Dnu_p - alpha_p * pow(nmax,2) / 2 - n0; // remove n0
        /*
        std::cout << " A=rfit[0] =" << rfit[0] << std::endl;
        std::cout << " B=rfit[1] =" << rfit[1] << std::endl;
        std::cout << " C=rfit[2] =" << rfit[2] << std::endl;
        std::cout << " alpha_p=" << alpha_p << std::endl;    
        std::cout << " nmax =" << nmax << std::endl;
        std::cout << " n0 = " << n0 << std::endl;
        std::cout << " epsilon_p =" << epsilon_p  << std::endl; // epsilon = (n0 + epsilon) - n0
        std::cout << " fl0_all=  " << fl0_all.transpose() << std::endl;
        std::cout << " fitted frequencies: " << std::endl;
        for (int i=0; i<xfit.size();i++){
            std::cout <<   rfit[0] +  rfit[1]*xfit[i] + rfit[2]*pow(xfit[i],2) << std::endl;
        }
        std::cout << " --- " << std::endl;
        //exit(EXIT_SUCCESS);
        std::cout << " asymptotic..." << std::endl;
        std::cout << "delta0l =" << delta0l << std::endl;
        std::cout << "DPl =" << DPl << std::endl;
        std::cout << "alpha_g =" << alpha_g << std::endl;
        std::cout << "q_star =" << q_star << std::endl;
        std::cout << "Hfactor =" << Hfactor << std::endl;
        std::cout << "rot_env =" << rot_env << std::endl;
        std::cout << "rot_core =" << rot_core << std::endl;
        std::cout << "fmin =" << fmin << "   fmax =" << fmax << std::endl;
        std::cout << " Getting inside solver:" << std::endl;
        */
        //freqs_l1=solve_mm_asymptotic_O2p(Dnu_p, epsilon_p, 1, delta0l, alpha_p, nmax, DPl, alpha_g, q_star, 0, fmin - Dnu_p, fmax + Dnu_p, step, true, false);
        freqs_l1=solve_mm_asymptotic_O2p(Dnu_p, epsilon_p, 1, delta0l, 0, 0., DPl, alpha_g, q_star, 0, fmin - Dnu_p, fmax + Dnu_p, step, true, false);
        //std::cout << " Getting outside solver " << std::endl;
        /*
        std::cout << "Dnu_p = " << Dnu_p << std::endl;
        std::cout << " freqs_l1.nu_p = " << freqs_l1.nu_p.transpose() << std::endl;
        std::cout << " freqs_l1.nu_g = " << freqs_l1.nu_g.transpose() << std::endl;
        std::cout << " freqs_l1.nu_m = " << freqs_l1.nu_m.transpose() << std::endl;
        std::cout << " --- " << std::endl;
        */
        //exit(EXIT_SUCCESS);
        } else{
        freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax);
    }
    fl1_all=freqs_l1.nu_m;

    // Adding the bias function to the vector fl1_all
     for (int i=0; i< fl1_all.size(); i++){
        fl1_all[i] = fl1_all[i] + bias(fl1_all[i]);
    }
   /*
    std::cout << " fl1     /   bias     /   Total " << std::endl;
    for (int i=0; i< fl1_all.size(); i++){
        std::cout << fl1_all[i] << "  ";
        fl1_all[i] = fl1_all[i] + bias(fl1_all[i]);
        std::cout << bias(fl1_all[i]) << "  "  << fl1_all[i] << std::endl;
    }
    */
    //exit(EXIT_SUCCESS);
    //if (sigma_m_l1 != 0){
    //   for (int en=0; en<fl1_all.size(); en++)
    //    {
    //         r = distrib_fl1m(gen_m);
    //         while (r > sigma_limit){
    //            r = distrib_fl1m(gen_m);
    //        }
    //        fl1_all[en]=fl1_all[en] + r;
    //    }
    //} 
     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg, Hfactor); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019).
   
    // We will add 0-anchors at the edges to avoid negative extrapolation... so we need at least two extra points... we choose to have 4 extra points
    f_interp.resize(fl0_all.size()+4); 
    f_interp[0]=fl0_all.minCoeff()*0.6; // A failsafe beyond the actual range of frequencies
    f_interp[1]=fl0_all.minCoeff()*0.8; 
    f_interp[f_interp.size()-2]=fl0_all.maxCoeff()*1.2;
    f_interp[f_interp.size()-1]=fl0_all.maxCoeff()*1.4;

    h_interp.resize(Hl0_all.size()+4); 
    h_interp[0]=0;
    h_interp[1]=Hl0_all[0]/4;
    h_interp[h_interp.size()-2]=Hl0_all[Hl0_all.size()-1]/4;
    h_interp[h_interp.size()-1]=0;
    for(int j=0; j<fl0_all.size(); j++){
        f_interp[j+2]=fl0_all[j];
        h_interp[j+2]=Hl0_all[j];
    }
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {     
        tmp=lin_interpol(f_interp, h_interp, fl1_all[i]); // interpolate Hl0 to fl1 positions
        if (tmp < 0){ 
            std::cout << "WARNING: WE IMPOSE Hl1p_all = 0 due to Negative interpolation" << std::endl;
            Hl1p_all[i]=0;
        } else{
            Hl1p_all[i]=std::abs(tmp);

        }
    }

    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    /*if (sigma_H_l1 != 0){
        for (int i=0; i<Hl1_all.size();i++){
             r = distrib_fl1m(gen_m);
             while ((1+r) < 0){ // Avoid negative heights
                r = distrib_fl1m(gen_m);
             }
             tmp=tmp*(1 + r); // Perturbation proportional to the actual value
        }
    }
    */
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1, Wfactor); // generate the mixed modes widths // Added on 28 Apr: Wfactor: in front of the ksi fonction: recommended to have a alpha ~ N(1,sigma=0.1)
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
    //eta=params[Nmax + lmax + Nf + 2];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
   
    // --------------
    // --------------
    // --------------
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    for (long n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]); 
        Wl2=Wl0_all[n];
        //std::cout << "Wl2=" << Wl2 << std::endl;
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l2[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }     
    }

//    std::cout << "--------------- " << std::endl;
    for(long n=0; n<Nfl3; n++){ 
        fl3=std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]);
        //std::cout << "fl3= " << fl3 << std::endl;
        Wl3=Wl0_all[n];
        if(do_amp){
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3/(pi*Wl3)*Vl3);
        } else{
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l3[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }  
    }           
//    std::cout << "--------------- " << std::endl;
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
        
        file_out="params.raw";
        std::ofstream outfile;
        outfile.open(file_out.c_str());
            outfile << "# File generated using " << modelname << std::endl;
            outfile << "# This is the raw vector of (1) params_length and (2) params as provided in input of the model" << std::endl;
            outfile << params_length.transpose() << std::endl;
            outfile << params.transpose() << std::endl;
        outfile.close();
    }
//    exit(EXIT_SUCCESS);
    return model_final;
}


VectorXd model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Width a following the Appourchaux et al. 2014, 566, 20 and Appourchaux et al. 2016, 595, C2 (Corrigendum) relation.
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * Like v2 and v3, the v4 version explicitly imposes numax as a parameter. 
     * Contrary to v3, it has the following capabilities:
     *      - model_type = 1 : The user can choose to use the l=0, shifted with d01 in order to generate l=1 p modes.
     *             In this case, we make use of:
             solve_mm_asymptotic_O2from_l0(nu_l0_in, el, delta0l, DPl, alpha_g, q, resol, bool returns_pg_freqs=true, verbose=false, freq_min=fmin, freq_max=fmax)
     *      - OR model_type=0: it can use an O2 polynomial:   
               solve_mm_asymptotic_O2p(Dnu_p, epsilon, el, delta0l, alpha_p, nmax, DPl, alpha_g, q, fmin, fmax, resol, true, false); //returns_pg_freqs=true, verbose=false
              Here, Dnu_p, epsilon are calculated using a linear fit of l=0 frequencies.
              alpha_p, nmax are determined using a 2nd order fit of l=0 frequencies.
     *      The switch between these two option is made using a control parameter at the end of vector of parameters (replacing sigma_limit of v2 and v3 and named model_type)
    *   The fit incorporate a list of Nerr parameters situated between ferr=[Nmax+ lmax + Nfl0 + 8] and [Nmax+ lmax + Nfl0 + 8 +Nerr] that are used to correct from the bias of the exact asymptotic model
    *   These are associated to a list of Nerr list of fixed frequencies between fref=[Nmax+ lmax + Nfl0 + 8 + Nerr] and [Nmax+ lmax + Nfl0 + 8 + 2 Nerr] that provided frequencies
    *   of (a) case of bias_type = 0 (cte):
               The reference frequencies fref are used as 'windows' for determining region of 'constant bias': If a mode lies within a window [fref(j),fref(j+1)] of the list, then we apply the correction ferr(j)
               This means that multiple modes may have the same correction factor ferr(j) and that some window may have 0. Be aware of this when post processing: Each window must be cleaned of solutions with 0 modes
           (b) case of bias_type = 1 (Cubic spline, ie strongly correlated as it is twice continuously differentiable):
               The reference frequencies fref are used as 'anchors' for a spline defined by the values ferr.
               The spline fitting is performed here in this function using the spline module from https://github.com/ttk592/spline/ also forked to my repo
           (c) case of bias_type = 2 (Hermite spline, ie less correlated variable as it is once continuously differentiable):
               The reference frequencies fref are used as 'anchors' for a spline defined by the values ferr.
               The spline fitting is performed here in this function using the spline module from https://github.com/ttk592/spline/ also forked to my repo
     * This model has therefore 6 parameters for the widths. 
     */
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = M_PI;//3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m + all the fl1p modes
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 5 here as it uses the Appourchaux profile
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const double sigma_limit=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];
    const double model_type=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+3];
    const double bias_type=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+4];
    const int Nferr=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+5];
    int i_dbg=0;
     
    VectorXd gamma_params(6);
    gamma_params << std::abs(params[Nmax + lmax + Nf + Nsplit + 0]) , std::abs(params[Nmax+lmax+Nf+Nsplit+1]) , std::abs(params[Nmax + lmax + Nf + Nsplit+2]),
            std::abs(params[Nmax+lmax+Nf+Nsplit+3]) , std::abs(params[Nmax+lmax+Nf+Nsplit+4]) , std::abs(params[Nmax+lmax+Nf+Nsplit+5]); 
   
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    Data_eigensols freqs_l1;
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    int Nharvey;
    long cpt;
    
    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    //outparams=true;
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=std::abs(params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]); 
 
    // Forcing values of visibilities to be greater than 0... priors will be in charge of the penalisation
    ratios_l0.setOnes();
    if(lmax >=1){
        Vl1=std::abs(params[Nmax]);
        ratios_l1=amplitude_ratio(1, inclination);
    }
    if(lmax >=2){
        Vl2=std::abs(params[Nmax+1]);
        ratios_l2=amplitude_ratio(2, inclination);

    }
    if(lmax >=3){
        Vl3=std::abs(params[Nmax+2]);
        ratios_l3=amplitude_ratio(3, inclination);
    }

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 

    for (int n=0; n<Nmax;n++)
    {
        lnGamma0=gamma_params[2] * log(fl0_all[n]/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl0_all[n]/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl0_all[n]=exp(lnGamma0 + lnLorentz);
    }

   if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       
    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=std::abs(params[Nmax + lmax + Nfl0 + 1]);
    const double alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2]);
    const double q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
    const double sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]);
    const double sigma_g_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]);
    const double Wfactor=std::abs(params[Nmax + lmax + Nfl0 + 6]);
    const double Hfactor=std::abs(params[Nmax + lmax + Nfl0 + 7]);
    std::vector<double> ferr_all(Nferr); // The total number fl1p modes is the total number of params Nfl1 - 8 // THESE MUST BE VARIABLE
    std::vector<double> fref_all(Nferr); // The total number fl1p modes is the total number of params Nfl1 - 8 // THESE MUST BE CONSTANT
    const double rot_env=std::abs(params[Nmax + lmax + Nf]);
    const double rot_core=std::abs(params[Nmax + lmax + Nf+1]);
    VectorXd ksi_pg, h1_h0_ratio,f_interp, h_interp;
 //   std::cout << params << std::endl;
 //   std::cout << " " << std::endl;
 //   std::cout << " i0=" << Nmax+ lmax + Nfl0 + 8 << std::endl;
 //   std::cout << " Nferr = " << Nferr << std::endl;
    for (int i=0;i<Nferr; i++){
        fref_all[i]=params[Nmax+ lmax + Nfl0 + 8 + i];        
        ferr_all[i]=params[Nmax+ lmax + Nfl0 + 8 + Nferr + i];
 //       std::cout << fref_all[i]  << "   " << ferr_all[i] << std::endl;
    }
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //std::default_random_engine gen_m(seed);    
    //std::normal_distribution<double> distrib_fl1m(0.,sigma_m_l1);
    //std::normal_distribution<double> distrib_Hl1m(0.,sigma_H_l1);
    
 //   exit(EXIT_SUCCESS);
    // bias type setup
    tk::spline bias; // declare a bias
    //tk::spline bias(fref_all, ferr_all);
    if (bias_type == 0){ // WINDOW CASE
        std::cout << " Windows case not implemented... Use cubic or hermite spline instead" << std::endl;
        exit(EXIT_SUCCESS);
    } else{
        bias.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0); // imposes second derivative to 0 at the edge (extrapolation if happens will be linear)
        if (bias_type == 1){ // CUBIC SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline);     // this calculates all spline coefficients   
        }
        if (bias_type == 2){ // HERMITE SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline_hermite);     // this calculates all spline coefficients
        }
    }
       // -------- DEBUG OF l=1 mixed modes -----
    /*
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    //std::cout << "sigma_p =" << sigma_p_l1 << std::endl;
    //std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    //std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "Hfactor =" << Hfactor << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    std::cout << "step =" << step << std::endl;
    */
    // ---------------------------------------

    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    //const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_nupl(fl1p_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    if (model_type == 0){
        VectorXd xfit, rfit;
        xfit=linspace(0, fl0_all.size()-1, fl0_all.size());
        rfit=linfit(xfit, fl0_all); // rfit[0] = Dnu
        const double Dnu_p=rfit[0];
        const long double n0=floor(rfit[0]/Dnu_p);
        const double epsilon_p = rfit[0]/Dnu_p -n0;
        //std::cout << " ---- 1st Order Fit ----" << std::endl;
        //std::cout << " Dnu_p =" << rfit[0] << std::endl;
        //std::cout << " epsilon_p + n0 =" << rfit[1]/rfit[0] << std::endl;
        //std::cout << " ---- 2nd Order Fit ----" << std::endl;
        //rfit=polyfitXd(xfit, fl0_all, 2); // rfit[0] = A, rfit[1]=B, rfit[2]=C : alpha_p=2 C / Dnu_p   ; nmax = (1 - B/Dnu_p)/alpha_p ; epsilon = A/Dnu_p - alpha_p * nmax^2 / 2
        //const long double alpha_p=2 *rfit[2]/Dnu_p;
        //const long double nmax = (1-rfit[1]/Dnu_p)/alpha_p;
        //const long double n0=floor(rfit[0]/Dnu_p - alpha_p * pow(nmax,2) / 2);
        //const long double epsilon_p = rfit[0]/Dnu_p - alpha_p * pow(nmax,2) / 2 - n0; // remove n0
        /*
        std::cout << " A=rfit[0] =" << rfit[0] << std::endl;
        std::cout << " B=rfit[1] =" << rfit[1] << std::endl;
        std::cout << " C=rfit[2] =" << rfit[2] << std::endl;
        std::cout << " alpha_p=" << alpha_p << std::endl;    
        std::cout << " nmax =" << nmax << std::endl;
        std::cout << " n0 = " << n0 << std::endl;
        std::cout << " epsilon_p =" << epsilon_p  << std::endl; // epsilon = (n0 + epsilon) - n0
        std::cout << " fl0_all=  " << fl0_all.transpose() << std::endl;
        std::cout << " fitted frequencies: " << std::endl;
        for (int i=0; i<xfit.size();i++){
            std::cout <<   rfit[0] +  rfit[1]*xfit[i] + rfit[2]*pow(xfit[i],2) << std::endl;
        }
        std::cout << " --- " << std::endl;
        //exit(EXIT_SUCCESS);
        std::cout << " asymptotic..." << std::endl;
        std::cout << "delta0l =" << delta0l << std::endl;
        std::cout << "DPl =" << DPl << std::endl;
        std::cout << "alpha_g =" << alpha_g << std::endl;
        std::cout << "q_star =" << q_star << std::endl;
        std::cout << "Hfactor =" << Hfactor << std::endl;
        std::cout << "rot_env =" << rot_env << std::endl;
        std::cout << "rot_core =" << rot_core << std::endl;
        std::cout << "fmin =" << fmin << "   fmax =" << fmax << std::endl;
        std::cout << " Getting inside solver:" << std::endl;
        */
        //freqs_l1=solve_mm_asymptotic_O2p(Dnu_p, epsilon_p, 1, delta0l, alpha_p, nmax, DPl, alpha_g, q_star, 0, fmin - Dnu_p, fmax + Dnu_p, step, true, false);
        freqs_l1=solve_mm_asymptotic_O2p(Dnu_p, epsilon_p, 1, delta0l, 0, 0., DPl, alpha_g, q_star, 0, fmin - Dnu_p, fmax + Dnu_p, step, true, false);
        //std::cout << " Getting outside solver " << std::endl;
        /*
        std::cout << "Dnu_p = " << Dnu_p << std::endl;
        std::cout << " freqs_l1.nu_p = " << freqs_l1.nu_p.transpose() << std::endl;
        std::cout << " freqs_l1.nu_g = " << freqs_l1.nu_g.transpose() << std::endl;
        std::cout << " freqs_l1.nu_m = " << freqs_l1.nu_m.transpose() << std::endl;
        std::cout << " --- " << std::endl;
        */
        //exit(EXIT_SUCCESS);
        } else{
        freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax);
    }
    fl1_all=freqs_l1.nu_m;

    // Adding the bias function to the vector fl1_all

//  std::cout << " fl1     /   bias     /   Total " << std::endl;
    for (int i=0; i< fl1_all.size(); i++){
        //std::cout << fl1_all[i] << "  ";
        fl1_all[i] = fl1_all[i] + bias(fl1_all[i]);
        //std::cout << bias(fl1_all[i]) << "  "  << fl1_all[i] << std::endl;
    }
    
    //exit(EXIT_SUCCESS);
    //if (sigma_m_l1 != 0){
    //   for (int en=0; en<fl1_all.size(); en++)
    //    {
    //         r = distrib_fl1m(gen_m);
    //         while (r > sigma_limit){
    //            r = distrib_fl1m(gen_m);
    //        }
    //        fl1_all[en]=fl1_all[en] + r;
    //    }
    //} 
     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg, Hfactor); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019).
   
    // We will add 0-anchors at the edges to avoid negative extrapolation... so we need at least two extra points... we choose to have 4 extra points
    f_interp.resize(fl0_all.size()+4); 
    f_interp[0]=fl0_all.minCoeff()*0.6; // A failsafe beyond the actual range of frequencies
    f_interp[1]=fl0_all.minCoeff()*0.8; 
    f_interp[f_interp.size()-2]=fl0_all.maxCoeff()*1.2;
    f_interp[f_interp.size()-1]=fl0_all.maxCoeff()*1.4;

    h_interp.resize(Hl0_all.size()+4); 
    h_interp[0]=0;
    h_interp[1]=Hl0_all[0]/4;
    h_interp[h_interp.size()-2]=Hl0_all[Hl0_all.size()-1]/4;
    h_interp[h_interp.size()-1]=0;
    for(int j=0; j<fl0_all.size(); j++){
        f_interp[j+2]=fl0_all[j];
        h_interp[j+2]=Hl0_all[j];
    }
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {     
        tmp=lin_interpol(f_interp, h_interp, fl1_all[i]); // interpolate Hl0 to fl1 positions
        if (tmp < 0){ 
            std::cout << "WARNING: WE IMPOSE Hl1p_all = 0 due to Negative interpolation" << std::endl;
            Hl1p_all[i]=0;
        } else{
            Hl1p_all[i]=std::abs(tmp);

        }
    }

    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    /*if (sigma_H_l1 != 0){
        for (int i=0; i<Hl1_all.size();i++){
             r = distrib_fl1m(gen_m);
             while ((1+r) < 0){ // Avoid negative heights
                r = distrib_fl1m(gen_m);
             }
             tmp=tmp*(1 + r); // Perturbation proportional to the actual value
        }
    }
    */
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1, Wfactor); // generate the mixed modes widths // Added on 28 Apr: Wfactor: in front of the ksi fonction: recommended to have a alpha ~ N(1,sigma=0.1)
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
    //eta=params[Nmax + lmax + Nf + 2];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
   
    // --------------
    // --------------
    // --------------
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    for (long n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]);
        lnGamma0=gamma_params[2] * log(fl2/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl2/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl2=exp(lnGamma0 + lnLorentz);
        //std::cout << "Wl2=" << Wl2 << std::endl;
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l2[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }     
    }

//    std::cout << "--------------- " << std::endl;
    for(long n=0; n<Nfl3; n++){ 
        fl3=std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]);
        lnGamma0=gamma_params[2] * log(fl3/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl3/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl3=exp(lnGamma0 + lnLorentz);
        if(do_amp){
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3/(pi*Wl3)*Vl3);
        } else{
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l3[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }  
    }           
//    std::cout << "--------------- " << std::endl;
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background

    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
        
        file_out="params.raw";
        std::ofstream outfile;
        outfile.open(file_out.c_str());
            outfile << "# File generated using " << modelname << std::endl;
            outfile << "# This is the raw vector of (1) params_length and (2) params as provided in input of the model" << std::endl;
            outfile << params_length.transpose() << std::endl;
            outfile << params.transpose() << std::endl;
        outfile.close();
    }
    //exit(EXIT_SUCCESS);
    return model_final;
}


VectorXd model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Width IS A SINGLE CONSTANT VALUE. THIS IS SUITABLE FOR FIT OVER A FEW Dnu
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * The v3 version is the same as the v2 version (explicitly impose numax as a parameter). 
     * But it requires the l=1 p modes as input to derive the mixed modes
     * This would lead to uncessary complications and increase model computation time
     * NOTE THAT IN THE RGB_asympt model, the v2 version is not implemented
     */
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m + all the fl1p modes
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 1
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const double sigma_limit=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];

    //int i_dbg=1; 
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    
    int Nharvey;
    long cpt;
    
    outparams=0;
    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    
    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=std::abs(params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]); 
 
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

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 

    Wl0_all.setConstant(std::abs(params[Nmax + lmax + Nf + Nsplit + 0]));
  
   if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       

    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
  
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=std::abs(params[Nmax + lmax + Nfl0 + 1]);
    const double alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2]);
    const double q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
    const double sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]);
    const double sigma_g_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]);
    const double sigma_m_l1=std::abs(params[Nmax + lmax + Nfl0 + 6]);
    const double Hfactor=std::abs(params[Nmax + lmax + Nfl0 + 7]);
    const VectorXd fl1p_all=params.segment(Nmax+ lmax + Nfl0 + 8, Nfl1 - 8); // The total number fl1p modes is the total number of params Nfl1 - 7
    const double rot_env=std::abs(params[Nmax + lmax + Nf]);
    const double rot_core=std::abs(params[Nmax + lmax + Nf+1]);
    VectorXd ksi_pg, h1_h0_ratio,f_interp, h_interp;
 
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen_m(seed);    
    std::normal_distribution<double> distrib_fl1m(0.,sigma_m_l1);
    std::normal_distribution<double> distrib_Hl1m(0.,sigma_H_l1);
    
    // -------- DEBUG OF l=1 mixed modes -----
    /*
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    std::cout << "sigma_p =" << sigma_p_l1 << std::endl;
    std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    std::cout << "step =" << step << std::endl;
    */
    // ---------------------------------------

    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_nupl(fl1p_all, 1, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    fl1_all=freqs_l1.nu_m;

    if (sigma_m_l1 != 0){
       for (int en=0; en<fl1_all.size(); en++)
        {
             r = distrib_fl1m(gen_m);
             while (r > sigma_limit){
                r = distrib_fl1m(gen_m);
            }
            fl1_all[en]=fl1_all[en] + r;
        }
    } 

     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg, Hfactor); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019).
   
    // We will add 0-anchors at the edges to avoid negative extrapolation... so we need at least two extra points... we choose to have 4 extra points
    f_interp.resize(fl0_all.size()+4); 
    f_interp[0]=fl0_all.minCoeff()*0.6; // A failsafe beyond the actual range of frequencies
    f_interp[1]=fl0_all.minCoeff()*0.8; 
    f_interp[f_interp.size()-2]=fl0_all.maxCoeff()*1.2;
    f_interp[f_interp.size()-1]=fl0_all.maxCoeff()*1.4;

    h_interp.resize(Hl0_all.size()+4); 
    h_interp[0]=0;
    h_interp[1]=Hl0_all[0]/4;
    h_interp[h_interp.size()-2]=Hl0_all[Hl0_all.size()-1]/4;
    h_interp[h_interp.size()-1]=0;
    for(int j=0; j<fl0_all.size(); j++){
        f_interp[j+2]=fl0_all[j];
        h_interp[j+2]=Hl0_all[j];
    }
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {     
        tmp=lin_interpol(f_interp, h_interp, fl1_all[i]); // interpolate Hl0 to fl1 positions
        if (tmp < 0){ 
            std::cout << "WARNING: WE IMPOSE Hl1p_all = 0 due to Negative interpolation" << std::endl;
            Hl1p_all[i]=0;
        } else{
            Hl1p_all[i]=std::abs(tmp);

        }
    }

    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    if (sigma_H_l1 != 0){
        for (int i=0; i<Hl1_all.size();i++){
             r = distrib_fl1m(gen_m);
             while ((1+r) < 0){ // Avoid negative heights
                r = distrib_fl1m(gen_m);
             }
             tmp=tmp*(1 + r); // Perturbation proportional to the actual value
        }
    }
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1); // generate the mixed modes widths
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
    //eta=params[Nmax + lmax + Nf + 2];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
   
    // --------------
    // --------------
    // --------------
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        //std::cout << "Wl0=" << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    for (long n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]);
        Wl2=params[Nmax + lmax + Nf + Nsplit + 0]; //Wl0_all[0];
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l2[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }
    }
//    std::cout << "--------------- " << std::endl;
    for(long n=0; n<Nfl3; n++){ 
        fl3=std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]);
        Wl3=params[Nmax + lmax + Nf + Nsplit + 0]; //Wl0_all[0];
        if(do_amp){
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3/(pi*Wl3)*Vl3); 
        } else{
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l3[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }  
    }           
//    std::cout << "--------------- " << std::endl;
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    //exit(EXIT_SUCCESS);
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    return model_final;
}




VectorXd model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    /* Model of the power spectrum of a Main sequence solar-like star
     * param is a vector of parameters
     * param_length defines the structure of the parameters
     * x is the frequency assumed to be in microHz
     * Width ARE FREE PARAMETERS
     * Warning: Although we have Nfli terms, all these MUST have same size, provided that 0<i<lmax. 
     *          Size MUST be 0 otherwise (this check is not made in this function)
     * The v3 version is the same as the v2 version (explicitly impose numax as a parameter). 
     * But it requires the l=1 p modes as input to derive the mixed modes
     * This would lead to uncessary complications and increase model computation time
     * This model has therefore 5 parameters for the widths. 
     * NOTE THAT IN THE RGB_asympt model, the v2 version is not implemented
     */
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const long double pi = 3.141592653589793238462643383279502884L;
    const int Nmax=params_length[0]; // Number of Heights
    const int lmax=params_length[1]; // number of degree - 1, ie, visibilities
    const int Nfl0=params_length[2]; // number of l=0 frequencies
    const int Nfl1=params_length[3]; // number of parameters to describe the l=1 mixed modes: delta0l, DPl, alpha_g, q, sigma_p, sigma_g, sigma_m + all the fl1p modes
    const int Nfl2=params_length[4]; // number of l=2 frequencies
    const int Nfl3=params_length[5]; // number of l=3 frequencies
    const int Nsplit=params_length[6]; // number of zones to describe the rotation profile (2 for rot_core and rot_renv + 1 eta + 1 a3 + 1 asym )
    const int Nwidth=params_length[7]; // number of parameters for the l=0 widths. Should be 5 here as it uses the Appourchaux profile
    const int Nnoise=params_length[8]; // number of parameters for the noise. Should be 7
    const int Ninc=params_length[9]; // number of parameters for the stellar inclination. Should be 1
    const int Ncfg=params_length[10]; // number of extra configuration parameters (e.g. truncation parameter)
    const int Nf=Nfl0+Nfl1+Nfl2+Nfl3;
    const double trunc_c=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc];
    const bool do_amp=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+1];
    const double sigma_limit=params[Nmax+lmax+Nf+Nsplit+Nwidth+Nnoise+Ninc+2];

    std::cout << "THIS MODEL HAS A PROBLEM: IT DOES NOT CONVERGE PROPERLY. DEBUG REQUIRED HERE AND ON PRIOR_APPLICATION" << std::endl;
    exit(EXIT_SUCCESS);
    int i_dbg=0;
     
    double inclination;

    VectorXd ratios_l0(1), ratios_l1(3), ratios_l2(5), ratios_l3(7);
    VectorXd model_l0(x.size()), model_l1(x.size()), model_l2(x.size()), model_l3(x.size()), model_noise(x.size()), model_final(x.size());

    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, Wl2_all(Nfl2), Wl3_all(Nfl3), a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3,eta0,a3, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    
    int Nharvey;
    long cpt;

    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited

    /*
       -------------------------------------------------------
       ------- Gathering information about the modes ---------
       -------------------------------------------------------
    */
    inclination=std::abs(params[Nmax + lmax + Nf+Nsplit + Nwidth + Nnoise]); 
 
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

    // --- Preparing profiles for l=0 modes ---
    fl0_all=params.segment(Nmax + lmax, Nfl0); // required for the interpolation of the widths and mixed modes determination 

    Wl0_all=params.segment(Nmax + lmax + Nf + Nsplit, Nwidth);

    //std::cout << "Wl0_all =" << Wl0_all << std::endl;

   if(do_amp){
        Hl0_all=params.segment(0, Nmax).cwiseProduct(Wl0_all.cwiseInverse() / pi); // A^2/(pi.Gamma)
        Hl0_all=Hl0_all.cwiseAbs();
    } else{
        Hl0_all=params.segment(0, Nmax).cwiseAbs();
    }       

    // --------------
    // --------------
    // --------------
    // --------------
    // ---- Mixed modes handling ----
    // --------------
  
    const double delta0l=params[Nmax + lmax + Nfl0];
    const double DPl=std::abs(params[Nmax + lmax + Nfl0 + 1]);
    const double alpha_g=std::abs(params[Nmax + lmax + Nfl0 + 2]);
    const double q_star=std::abs(params[Nmax + lmax + Nfl0 + 3]);
    const double sigma_H_l1=std::abs(params[Nmax + lmax + Nfl0 + 4]);
    const double sigma_g_l1=std::abs(params[Nmax + lmax + Nfl0 + 5]);
    const double sigma_m_l1=std::abs(params[Nmax + lmax + Nfl0 + 6]);
    const double Hfactor=std::abs(params[Nmax + lmax + Nfl0 + 7]);
    const VectorXd fl1p_all=params.segment(Nmax+ lmax + Nfl0 + 8, Nfl1 - 8); // The total number fl1p modes is the total number of params Nfl1 - 7
    const double rot_env=std::abs(params[Nmax + lmax + Nf]);
    const double rot_core=std::abs(params[Nmax + lmax + Nf+1]);
    
    VectorXd ksi_pg, h1_h0_ratio,f_interp, h_interp;
 
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen_m(seed);    
    std::normal_distribution<double> distrib_fl1m(0.,sigma_m_l1);
    std::normal_distribution<double> distrib_Hl1m(0.,sigma_H_l1);
    
    // -------- DEBUG OF l=1 mixed modes -----
    /*
    std::cout << "delta0l =" << delta0l << std::endl;
    std::cout << "DPl =" << DPl << std::endl;
    std::cout << "alpha_g =" << alpha_g << std::endl;
    std::cout << "q_star =" << q_star << std::endl;
    std::cout << "sigma_p =" << sigma_H_l1 << std::endl;
    std::cout << "sigma_g =" << sigma_g_l1 << std::endl;
    std::cout << "sigma_m =" << sigma_m_l1 << std::endl;
    std::cout << "rot_env =" << rot_env << std::endl;
    std::cout << "rot_core =" << rot_core << std::endl;
    std::cout << "Hfactor =" << Hfactor << std::endl;
    std::cout << "step =" << step << std::endl;
    */
      // ---------------------------------------    
    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_nupl(fl1p_all, 1, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    fl1_all=freqs_l1.nu_m;

    if (sigma_m_l1 != 0){
       for (int en=0; en<fl1_all.size(); en++)
        {
             r = distrib_fl1m(gen_m);
             while (r > sigma_limit){
                r = distrib_fl1m(gen_m);
            }
            fl1_all[en]=fl1_all[en] + r;
        }
    }     
     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg, Hfactor); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019). Note that Factor must be in [0,1]
    
    // We will add 0-anchors at the edges to avoid negative extrapolation... so we need at least two extra points... we choose to have 4 extra points
    f_interp.resize(fl0_all.size()+4); 
    f_interp[0]=fl0_all.minCoeff()*0.6; // A failsafe beyond the actual range of frequencies
    f_interp[1]=fl0_all.minCoeff()*0.8; 
    f_interp[f_interp.size()-2]=fl0_all.maxCoeff()*1.2;
    f_interp[f_interp.size()-1]=fl0_all.maxCoeff()*1.4;

    h_interp.resize(Hl0_all.size()+4); 
    h_interp[0]=0;
    h_interp[1]=Hl0_all[0]/4;
    h_interp[h_interp.size()-2]=Hl0_all[Hl0_all.size()-1]/4;
    h_interp[h_interp.size()-1]=0;
    for(int j=0; j<fl0_all.size(); j++){
        f_interp[j+2]=fl0_all[j];
        h_interp[j+2]=Hl0_all[j];
    }
    Hl1p_all.resize(fl1_all.size());
    for (int i=0; i<fl1_all.size();i++)
    {     
        tmp=lin_interpol(f_interp, h_interp, fl1_all[i]); // interpolate Hl0 to fl1 positions
        if (tmp < 0){ 
            std::cout << "WARNING: WE IMPOSE Hl1p_all = 0 due to Negative interpolation" << std::endl;
            Hl1p_all[i]=0;
        } else{
            Hl1p_all[i]=std::abs(tmp);

        }
    }
    
    Hl1_all=h1_h0_ratio.cwiseProduct(Hl1p_all*Vl1);
    if (sigma_H_l1 != 0){
        for (int i=0; i<Hl1_all.size();i++){
             r = distrib_fl1m(gen_m);
             while ((1+r) < 0){ // Avoid negative heights
                r = distrib_fl1m(gen_m);
             }
             tmp=tmp*(1 + r); // Perturbation proportional to the actual value
        }
    }
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1); // generate the mixed modes widths
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
    //eta=params[Nmax + lmax + Nf + 2];
    eta0=eta0_fct(fl0_all);
    a3=params[Nmax + lmax + Nf + 3];
    asym=params[Nmax+lmax + Nf + 4];
   
    //bool error, crash;
    //error=debug_solver(x, fl1_all, fl0_all, 1, delta0l, DPl, alpha_g, q_star, sigma_p_l1);
    //if (error == true){
    //    std::cout << " Full parameter list: " << params.transpose() << std::endl;
    //    std::cout << " Full parameters_length: " << params_length.transpose() << std::endl;
    //    crash=debug(model_final, Hl1, fl1, -1, eta, a3, asym, Wl1, 1, step, inclination, ratios_l1, trunc_c, true);
    //}

    // --------------
    // --------------
    // --------------
    // --------------
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 
    cpt=0;
    for(long n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        //std::cout << "Wl0=" << Wl0 << std::endl;
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, a3, asym, Wl0, 0, ratios_l0, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
        //crash=debug(model_final, Hl0, fl0, 0, eta, a3, asym, Wl0, 0, step, inclination, ratios_l0, trunc_c, true);
    }
    for (long n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, a3,asym, Wl1, 1, ratios_l1, step, trunc_c);
        if (outparams){
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1_l1[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    for(long n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]);         
        Wl2=lin_interpol(fl0_all, Wl0_all, fl2); // interpolate Wl0 to fl2 positions
        Wl2=std::abs(Wl2);
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        //crash=debug(model_final, Hl2, fl2, a1_l2[n], eta, a3, asym, Wl2, 2, step, inclination, ratios_l2, trunc_c, true);
        if (outparams){
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1_l2[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }

//    std::cout << "--------------- " << std::endl;
    for(long n=0; n<Nfl3; n++){ 
        fl3=std::abs(params[Nmax+lmax+Nfl0+Nfl1+Nfl2+n]);
        Wl3=lin_interpol(fl0_all, Wl0_all, fl3); // interpolate Wl0 to fl2 positions
        Wl3=std::abs(Wl3);
          if(do_amp){
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3/(pi*Wl3)*Vl3);
        } else{
            Hl3=lin_interpol(fl0_all, Hl0_all, fl3);
            Hl3=std::abs(Hl3*Vl3);            
        }       
        model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        //crash=debug(model_final, Hl3, fl3, a1_l3[n], eta, a3, asym, Wl3, 3, step, inclination, ratios_l3, trunc_c, true);
        if (outparams){
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1_l3[n], eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }  
    }           
//    std::cout << "--------------- " << std::endl;
    /* -------------------------------------------------------
       ------- Gathering information about the noise ---------
       -------------------------------------------------------
    */
    //std::cout << "Before computing noise model" << std::endl;
    noise_params=params.segment(Nmax+lmax+Nf+Nsplit+Nwidth, Nnoise);
    Nharvey=(Nnoise-1)/3;
    //std::cout << "Nharvey = " << Nharvey << std::endl;
        
    /* -------------------------------------------------------
       ---------- Computing the mode of the noise ------------
       -------------------------------------------------------
    */
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey); // this function increment the model_final with the noise background
    //debug(model_final, -1, -1, -1, -1, -1, -1, -1, 0, step, -1, ratios_l0, trunc_c, true);
 
    if(outparams){
        int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input mode parameters. degree / freq / H / W / a1 / eta0*1e-6 / a3 / asymetry / inclination";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    return model_final;
}


VectorXd model_Harvey_Gaussian(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
/*
 * A model in which we fit a Gaussian + Harvey-Like profile.
 * Parameters are assumed to be in that order: Maximum Height, variance/width, central frequency, White noise
*/
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const int Nrows=1, Ncols=3; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested

	VectorXd nu0(x.size()), model_final(x.size()), noise_params;
	int c=0;
    int Nharvey;

    /*if (outparams){
        std::cout << " ERROR: outparams is not implemented for  " << __func__ << std::endl;
        std::cout << "       If you want to get outputs in ascii format, you must review this function and implement outparams = true" << std::endl;
        std::cout << "       The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    } 
    */
	// ------ Setting the Gaussian -------
	model_final= -0.5 * (x - nu0.setConstant(params[8])).array().square() /pow(std::abs(params[9]),2);
	model_final= std::abs(params[7])*model_final.array().exp();
	// ----------------------------------

    // ---- Setting the Noise model -----
    Nharvey=2;
    noise_params=params.segment(0, 7); // pick the first 7 elements, begining from the index 3: [H1, tc1, p1, H2, tc2, p2, B0]
    model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey);
    // ----------------------------------

    if(outparams){
        //int c=0;
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input Gaussian envelope parameters. Amax / numax / sigma";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                std::cout << "i,j : " << i << ", " << j << std::endl;
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        std::cout << noise << std::endl;
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, 0, mode_params.cols()), noise, file_out, modelname, name_params);
    }
	return model_final;
}

VectorXd model_Harvey1985_Gaussian(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
/*
 * A model in which we fit a Gaussian + Harvey-Like profile.
 * Parameters are assumed to be in that order: Maximum Height, variance/width, central frequency, White noise
*/
    const double step=x[2]-x[1]; // used by the function that optimise the lorentzian calculation
    const int Nrows=1, Ncols=3; // Number of parameters for each mode
    int c=0;
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested

	VectorXd nu0(x.size()), model_final(x.size()), noise_params;
	int Nharvey;

	// ------ Setting the Gaussian -------
	model_final= -0.5 * (x - nu0.setConstant(params[2])).array().square() /pow(std::abs(params[1]),2);
	model_final= std::abs(params[0])*model_final.array().exp();
	// ----------------------------------

	// ---- Setting the Noise model -----
	Nharvey=1;
	noise_params=params.segment(3, 4); // pick the 3 elements, begining from the index 3
	model_final=harvey1985(noise_params.array().abs(), x, model_final, Nharvey);
	//model_final=harvey_like(noise_params.array().abs(), x, model_final, Nharvey);
	// ----------------------------------

    if(outparams){
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input Gaussian envelope parameters. Amax / numax / sigma";
        VectorXd spec_params(3);
        spec_params << x.minCoeff() , x.maxCoeff(), step;
        MatrixXd noise(Nharvey+1, 3);
        noise.setConstant(-2);
        for(int i=0;  i<noise.cols(); i++){
            for(int j=0;j<Nharvey; j++){
                noise(i,j)=noise_params(c);
                c=c+1;
           }
        }
        noise(Nharvey, 0) = noise_params(c); // White noise 
        std::cout << "here" << std::endl;
        write_star_params(spec_params, params, params_length, mode_params.block(0,0, 1, mode_params.cols()), noise, file_out, modelname, name_params);
    }
	return model_final;
}


VectorXd model_Test_Gaussian(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
/*
 * A model in which we fit a Gaussian + White noise.
 * Parameters are assumed to be in that order: Maximum Height, variance/width, central frequency, White noise
*/
	VectorXd nu0(x.size()), model_final(x.size()), tmp(x.size());

	tmp.setConstant(params[3]);

    if (outparams){
        std::cout << " ERROR: outparams is not implemented for  " << __func__ << std::endl;
        std::cout << "       If you want to get outputs in ascii format, you must review this function and implement outparams = true" << std::endl;
        std::cout << "       The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    } 

	model_final= -0.5 * (x - nu0.setConstant(params[2])).array().square() /pow(params[1],2);
	model_final= params[0]*model_final.array().exp();
	model_final= model_final + tmp;
	return model_final;
}


VectorXd model_ajfit(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
    // Model to fit the aj terms. It is a simplified version of the python code within fit_a2sig.py::do_stats_ongrid_for_observations()
    // It handles mean values aj coefficients
    outparams=true;
    //
    const int Nrows=1000, Ncols=10; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int Line=0; // will be used to trim the mode_params table where suited
    //
    const int lmax=3;
    const std::string data_type="mean_nu_l"; 
    const std::string filter_type="gate";
    const bool do_a2=params[params_length.segment(0, 4).sum()];
    const bool do_a4=params[params_length.segment(0, 4).sum()+1];
    const bool do_a6=params[params_length.segment(0, 4).sum()+2];
    const bool do_CFonly=params[params_length.segment(0,4).sum()+3];
    int i;
    long double a2_mod_data=-9999;
    long double a4_mod_data=-9999;
    long double a6_mod_data=-9999;
    VectorXi posl, pos;
    VectorXd model_final(x.size());
    VectorXd acoefs;
    VectorXd a2_nl_mod(params_length[2]), a4_nl_mod(params_length[2]), a6_nl_mod(params_length[2]);
    VectorXd thetas(2), els, nu_nl_obs;
    long double epsilon_nl0, eta0, a1_obs, Dnu_obs;
    
    // Deploying the parameters
    epsilon_nl0 = params[0];
    thetas << params[1]*M_PI/180., params[2]*M_PI/180.;
    Dnu_obs= params[3];
    a1_obs=params[4];
    els = params.segment(5, params_length[2]);
    nu_nl_obs = params.segment(5+params_length[2], params_length[3]);
    eta0=eta0_fct(Dnu_obs);
    //std::cout << "els[i]  / nu_nl_obs    / eta0   /  a1_obs    / epsilon_nl0   / thetas[0]    / thetas[1]    / filter_type" << std::endl;
    for (int i=0; i<params_length[2];i++){
        //std::cout << els[i]  << "  "  << nu_nl_obs[i] << "   " <<  eta0  <<  "   " << a1_obs   << "  " << epsilon_nl0  <<  "  "  << thetas[0] << "  " << thetas[1]  << "  " << filter_type << std::endl;
        if(do_CFonly == false){
            acoefs=decompose_Alm_fct(els[i], nu_nl_obs[i], eta0, a1_obs, epsilon_nl0, thetas, filter_type); // eta0 will be used only if a1_obs != 0 (case where CF was not removed from data)
        } else{
            acoefs=decompose_CFonly(els[i],  nu_nl_obs[i], eta0, a1_obs); // basically returns only a2_CF. The rest is set to by definition
        }
        a2_nl_mod[i]=acoefs[1]*1e3; //  convert a2 in nHz, because we assume here that nu_nl is in microHz
        a4_nl_mod[i]=acoefs[3]*1e3; //  convert a4 in nHz, because we assume here that nu_nl is in microHz
        a6_nl_mod[i]=acoefs[5]*1e3; //  convert a6 in nHz, because we assume here that nu_nl is in microHz
        if (outparams){
            mode_params.row(Line) << els[i], nu_nl_obs[i], epsilon_nl0 , thetas[0], thetas[1], eta0*1e-6,  a1_obs, a2_nl_mod[i], a4_nl_mod[i], a6_nl_mod[i];// mode_vec;
            Line=Line+1;
        }      }
        // Calculate the mean over nu and l of the a-coefficients
    if (data_type == "mean_nu_l"){ // Here, all of the data we use for the posterior are the results of the mean of <aj>_ln = cte
        if (do_a2 == true){
            a2_mod_data=a2_nl_mod.mean();
            pos=where_dbl(x, 2, 1e-3); // Locate where the a2 coefficient is in the x-vector
            model_final[pos[0]]=a2_mod_data;    
        }
        if (do_a4 == true){
            posl=where_in_range(els, 2, lmax, 0); // a4 makes sense only for l>=2, we need to filter out the values to select only l>=2 before the fit
            if (posl[0] != -1){
                a4_mod_data=0;
                for (int k=0;k<posl.size(); k++){
                    a4_mod_data=a4_mod_data + a4_nl_mod[posl[k]]/posl.size();
                }
            }
            //std::cout << "a4_mod_data =" << a4_mod_data << std::endl;
            pos=where_dbl(x, 4, 1e-3); // Locate where the a4 coefficient is in the x-vector
            model_final[pos[0]]=a4_mod_data;
        }
        if (do_a6 == true){
            posl=where_in_range(els, 3, lmax, 0); // a4 makes sense only for l>=2, we need to filter out the values to select only l>=2 before the fit
            if (posl[0] != -1){
                a6_mod_data=0;
                for (int k=0;k<posl.size(); k++){
                    a6_mod_data=a6_mod_data + a6_nl_mod[posl[k]]/posl.size();
                }
            }
            pos=where_dbl(x, 6, 1e-3); // Locate where the a2 coefficient is in the x-vector
            model_final[pos[0]]=a6_mod_data;
        }
    }
    if(outparams){
        std::string file_out="params.model";
        std::string modelname = __func__;
        std::string name_params = "# Input aj model parameters. l / nu_nl_obs  / epsilon_nl0 / theta0 / delta /  eta0   / a1_obs  / a2_nl   / a4_nl  /  a6_nl";
        VectorXd do_aj(3);
        do_aj << do_a2 , do_a4, do_a6;
        MatrixXd noise(1, 3);
        noise(0,0)=a2_mod_data; noise(0,1)=a4_mod_data; noise(0,2)=a6_mod_data;
        write_star_params(do_aj, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, file_out, modelname, name_params);
    }
    return model_final;
}



///// ----------- FOR ASCII OUTPUTS  ------------ //////
void write_star_params(const VectorXd& spec_params, const VectorXd& raw_params, const VectorXi& plength, const MatrixXd& mode_params, const MatrixXd& noise_params, 
    const std::string file_out, const std::string modelname, const std::string name_params){

    VectorXi Nchars_spec(3), Nchars_params(17), Nchars_noise(3), precision_spec(3), precision_params(17), precision_noise(3);

    std::ofstream outfile;
    outfile.open(file_out.c_str());

    if(outfile.is_open()){
        // ---------------------
        Nchars_spec << 20, 20 , 20;
        precision_spec << 10, 10, 10;

        //outfile << "ID= " << identifier << std::endl;
        outfile << "Model = " << modelname << std::endl;
        outfile << "plength = ";
        for(int i=0;i<plength.size(); i++){
            outfile << std::setw(8) << std::setprecision(0) << plength(i);
        }
        outfile << std::endl;
        outfile << "raw_params = ";
        for(int i=0;i<raw_params.size(); i++){
            outfile << std::setw(20) << std::setprecision(8) << raw_params(i);
        }
        outfile << std::endl;
        outfile << "# Spectrum parameters. Frequency range min/max (microHz) / Resolution (microHz)" << std::endl;
        outfile << std::setw(Nchars_spec[0]) << std::setprecision(precision_spec[0]) << spec_params[0];
        outfile << std::setw(Nchars_spec[1]) << std::setprecision(precision_spec[1]) << spec_params[1];
        outfile << std::setw(Nchars_spec[2]) << std::setprecision(precision_spec[2]) << spec_params[2] << std::endl;


        // ---------------------
        outfile << "# Configuration of mode parameters. This file was generated by models.cpp (TAMCMC code version >1.61)" << std::endl;
        Nchars_params    << 5, 20, 20, 20, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 , 16; // Handle max 15 input here. The only constrain is to have deg/freq/H/W first
        precision_params << 1, 10, 10, 10, 8, 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8 , 8, 8 , 8;
        outfile << name_params << std::endl; //"# Input mode parameters. degree / freq / H / W / splitting a1 / a2_0  / a2_1  / a2_2 / a3 / asymetry / inclination" << std::endl;   
        
        for(int i=0; i<mode_params.rows(); i++){
            for(int j=0;j<mode_params.cols(); j++){
                outfile << std::setw(Nchars_params[j]) << std::setprecision(precision_params[j]) << mode_params(i,j);
            }
            outfile << std::endl;
        }

        // ---------------------
        Nchars_noise << 16, 16, 16;
        precision_noise << 6, 6, 6;

        outfile << "# Configuration of noise parameters. This file was generated by write_star_mode_params (write_star_params.cpp)" << std::endl;
        outfile << "# Input mode parameters. H0 , tau_0 , p0 / H1, tau_1, p1 / N0. Set at -1 if not used. -2 means that the parameter is not even written on the file (because irrelevant)." << std::endl;

        for(int i=0; i<noise_params.rows(); i++){
            for(int j=0;j<noise_params.cols(); j++){
                outfile << std::setw(Nchars_noise[j]) << std::setprecision(precision_noise[j]) << noise_params(i,j);
            }
            outfile << std::endl;
        }
        outfile.close();
    }  
    else {
        std::cout << " Unable to open file " << file_out << std::endl;  
        std::cout << " Check that the full path exists" << std::endl;
        std::cout << " The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    }
}


///// ----------- FOR DEBUG ------------/////

bool debug(const VectorXd& model, const long double Hl, const long double fl, const long double a1, const long double eta, const long double a3,
               const long double asym, const long double Wl, const long double el, const long double step, const double inclination, const VectorXd& ratios,
               const long double trunc_c, const bool exit_c){
// This function must be put after an 'optimum_lorentzian_Calc_a1etaa3()' function to verify that the set of paraneters given to that same function
// Is not responsible of NaN in the likelihood, which can occur when models have negative entries of if the model itself contains NaN.
// All the input parameters have the usual definition for a mode fit, execpt 'exit' that if set to true (default) force program exit

    const long double mini=1e10;
    const long double maxi=0;

    bool leadtoNaN;
    VectorXi posNotOK;

    for (int i=0; i<model.size();i++){
        if (model.array().isFinite()[i] != true){
            leadtoNaN=true;
            if (model.array().isNaN()[i] == true){
                std::cout << "NaN detected into the model: " << std::endl;
            } else{
                if (model.array().isInf()[i] == true){
                    std::cout << "Inf detected into the model:" << std::endl;
                } else{
                    std::cout << "Something wrong in model leading isFinite() to say that it is not finite, but it is neither a NaN or an Inf" << std::endl;
                }
            }
            goto summary;
        }
    }

    posNotOK=where_in_range(model, mini, maxi, 1);
    if (posNotOK[0] != -1){
        leadtoNaN=true;
        std::cout << "Negative values in the model detected:" << std::endl;
    }

    summary:
    if(leadtoNaN == true){
        std::cout << " The error happens for the following parameter condition of the last calculated mode :" << std::endl;
        std::cout << " l= " << el << std::endl;
        std::cout << " fl= " << fl << std::endl;
        std::cout << " Hl= " << Hl << std::endl;
        std::cout << " Wl= " << Wl << std::endl;
        std::cout << " a1= " << a1 << std::endl;
        std::cout << " eta= " << eta << std::endl;
        std::cout << " a3= " << a3 << std::endl;
        std::cout << " asym= " << asym << std::endl;
        std::cout << " inclination= " << inclination << std::endl;
        std::cout << " ratios= " << ratios.transpose() << std::endl;
        std::cout << " step= " << step << std::endl;
        std::cout << " trunc_c= " << trunc_c << std::endl;   
    }
    if (exit_c == true && leadtoNaN == true){
        std::cout << " Exit requested on debug() due to the detection of a NaN condition for the likelihood " << std::endl;
        exit(EXIT_FAILURE);
    }
    return leadtoNaN;
}

bool debug_solver(const VectorXd& x, const VectorXd& fl1_all, const VectorXd& fl0_all, const int el, const long double delta0l, 
                  const long double DPl, const long double alpha_g, const long double q_star, const long double sigma_p_l1){

    bool error;

    if ( fl1_all.minCoeff() < fl0_all.minCoeff() || fl1_all.minCoeff() < x.minCoeff() ||
         fl1_all.maxCoeff() > fl0_all.maxCoeff() || fl1_all.maxCoeff() > x.maxCoeff()){
        error=true;
        std::cout << " Detected out of bound solution that may cause issues: " << std::endl;
        std::cout << " fmin =" << x.minCoeff() << std::endl;
        std::cout << " fmacx =" << x.maxCoeff() << std::endl;
        std::cout << " el = " <<  el << std::endl;
        std::cout << " fl0_all = " << fl0_all.transpose() << std::endl;
        std::cout << " fl1_all = " << fl1_all.transpose() << std::endl;
        std::cout << " delta0l = " << delta0l << std::endl;
        std::cout << " DPl = " << DPl << std::endl;
        std::cout << " alpha_g = " << alpha_g << std::endl;
        std::cout << " q_star = " << q_star << std::endl;
        std::cout << " sigma_p_l1 = " << sigma_p_l1 << std::endl;
    }
    return error;
}


// ------ Common ------
/*
double eta0_fct(const VectorXd& fl0_all){
    const double G=6.667e-8;
    const double Dnu_sun=135.1;
    const double R_sun=6.96342e5; //in km
    const double M_sun=1.98855e30; //in kg
    const double rho_sun=M_sun*1e3/(4*M_PI*std::pow(R_sun*1e5,3)/3); //in g.cm-3
    double rho, eta0;
    VectorXd xfit, rfit;
    xfit=linspace(0, fl0_all.size()-1, fl0_all.size());
    rfit=linfit(xfit, fl0_all); // rfit[0] = Dnu 
    rho=pow(rfit[0]/Dnu_sun,2.) * rho_sun;
    //eta0=3./(4.*M_PI*rho*G); WRONG DUE TO WRONG ASSUMPTION OMEGA ~ a1. It SHOULD BE OMEGA ~ 2.pi.a1
    eta0=3.*M_PI/(rho*G);
    return eta0;
}
*/
double eta0_fct(const VectorXd& fl0_all){
    VectorXd xfit, rfit;
    xfit=linspace(0, fl0_all.size()-1, fl0_all.size());
    rfit=linfit(xfit, fl0_all); // rfit[0] = Dnu 
    return eta0_fct(rfit[0]);
}

double eta0_fct(const double Dnu_obs){
    const double G=6.667e-8;
    const double Dnu_sun=135.1;
    const double R_sun=6.96342e5; //in km
    const double M_sun=1.98855e30; //in kg
    const double rho_sun=M_sun*1e3/(4*M_PI*std::pow(R_sun*1e5,3)/3); //in g.cm-3
    double rho, eta0;
    rho=pow(Dnu_obs/Dnu_sun,2.) * rho_sun;
    //eta0=3./(4.*M_PI*rho*G); WRONG DUE TO WRONG ASSUMPTION OMEGA ~ a1. It SHOULD BE OMEGA ~ 2.pi.a1
    eta0=3.*M_PI/(rho*G);
    return eta0;
}
/*
    This function takes for argument elements necessary for:
        - accounting of the centrifugal effects: eta0, a1
        - accounting of the activity using Gizon prescription for activity: epsilon_nl, thetas=[theta0, delta], filter_type
    And it perturbates accordingly the central frequency fc_l of a given mode of degree l
    Note: This will only generate valid even a-coefficients (up to 6). Because there is no elements
           of perturbation on even a-coefficients, a1, a3, a5 will be returned as 0. DO NOT USE THEM IN THE MODELING
*/
VectorXd decompose_Alm_fct(const int l, const long double fc_l, const long double eta0, const long double a1, const long double epsilon_nl, const VectorXd& thetas, const std::string filter_type){
    VectorXd aj, nu_nlm(2*l+1);
    for (int m=-l; m<=l; m++){
        nu_nlm[m+l]=fc_l;
        if (eta0 > 0){
            nu_nlm[m+l] = nu_nlm[m+l] + fc_l*eta0*Qlm(l,m)*pow(a1*1e-6,2);
        } 
        nu_nlm[m+l]=nu_nlm[m+l] + fc_l*epsilon_nl*Alm(l, m, thetas[0], thetas[1], filter_type);
    }
    //std::cout << " Alm(l, m, thetas[0], thetas[1], filter_type) = " << Alm(l, m, thetas[0], thetas[1], filter_type) << std::endl;
    //std::cout << " epsilon_nl = " << epsilon_nl << std::endl;
    //std::cout << "fc_l= " << fc_l << "  nu_nlm =" << nu_nlm.transpose() << std::endl;
    aj=eval_acoefs(l, nu_nlm);
    return aj;
}

VectorXd decompose_CFonly(const int l, const long double fc_l, const long double eta0, const long double a1){
    /*
        This is for models that seek to evaluate only the Centrifugal force effect on a2
        Basically, this is mostly usefull if one wants to evaluate the global likelihood of 
        a result being consistent with a centrifugal distorsion only
    */
    VectorXd aj, nu_nlm(2*l+1);
    for (int m=-l; m<=l; m++){
        nu_nlm[m+l]=fc_l;
        if (eta0 > 0){
            nu_nlm[m+l] = nu_nlm[m+l] + fc_l*eta0*Qlm(l,m)*pow(a1*1e-6,2);
        }
    }
    aj=eval_acoefs(l, nu_nlm);
    return aj;
}