# include <iostream>
# include <iomanip>
#include <fstream>
# include <Eigen/Dense>
#include "models.h"
#include "noise_models.h"
#include "data.h"

#include "../../external/ARMM/solver_mm.h"
#include "../../external/ARMM/bump_DP.h"
#include "../../external/Alm/Alm_cpp/activity.h"
#include "../../external/Alm/Alm_cpp/Alm_interpol.h"
#include "../../external/Alm/Alm_cpp/data.h"
#include "acoefs.h"
#include "interpol.h"
#include "linfit.h"
#include "polyfit.h"
#include "../../external/spline/src/spline.h"

using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXd amplitude_ratio(int l, double beta);


////////////////////////////////////
// Models for asymptotic fitting ///
////////////////////////////////////

VectorXd model_RGB_asympt_aj_CteWidth_HarveyLike_v4_parallel(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
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
               The reference frequencies fref are NOT used. You basically fit the asymptotic. Note that the setups should have imposed ferr to fix in that scenario to avoid loss of efficiency while fitting
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
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3;
    double a1,eta0,a3, asym;
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

    // bias type setup. Note that you need at least 3 points to define a spline. Otherwise it will crash here.
    tk::spline bias; // declare a bias
    if (bias_type != 0){ // 0 is No bias case. Defined when priors are set to fix in the hyper priors
        bias.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0); // imposes second derivative to 0 at the edge (extrapolation if happens will be linear)
        if (bias_type == 1){ // CUBIC SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline);     // this calculates all spline coefficients   
        }
        if (bias_type == 2){ // HERMITE SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline_hermite);     // this calculates all spline coefficients
        }
    }

    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    //const Data_eigensols freqs_l1=solve_mm_asymptotic_O2from_nupl(fl1p_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax); // note that we use the true data resolution (step) for optimising computation
    if (model_type == 0){
        VectorXd xfit, rfit;
        xfit = Eigen::VectorXd::LinSpaced(fl0_all.size(), 0, fl0_all.size()-1);
        rfit=linfit(xfit, fl0_all); // rfit[0] = Dnu
        const double Dnu_p=rfit[0];
        const int n0=floor(rfit[1]/Dnu_p);
        const double epsilon_p = rfit[1]/Dnu_p -n0;
        if(fmin - Dnu_p <0){
            std::cerr << " ERROR: THE ARMM WILL NOT CONVERGE AS YOU HAVE fmin - Dnu_p < 0, leading the an INFINITE g modes density!" << std::endl;
            std::cerr << "           You need to consider a higher fmin and/or smaller Dnu_p..." << std::endl;
            std::cerr << "        fmin  : " << fmin << std::endl;
            std::cerr << "        Dnu_p : " << Dnu_p << std::endl;
            exit(EXIT_FAILURE);
        }    
        freqs_l1=solve_mm_asymptotic_O2p(Dnu_p, epsilon_p, 1, delta0l, 0, 0., DPl, alpha_g, q_star, 0, fmin - Dnu_p, fmax + Dnu_p, step, true, false);
    } else{
        freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax);
    }
    fl1_all=freqs_l1.nu_m;

    // Adding the bias function to the vector fl1_all
    if (bias_type !=0){
        VectorXd b;
        b.resize(fl1_all.size());
        b.setZero();
        for (int i=0; i< fl1_all.size(); i++){
            b[i]=bias(fl1_all[i]);
            fl1_all[i] = fl1_all[i] + b[i];
            //std::cout << bias(fl1_all[i]) << "  "  << fl1_all[i] << std::endl;
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
    model_final.setZero();
    
    /* -------------------------------------------------------
       --------- Computing the models for the modes  ---------
       -------------------------------------------------------
    */ 

    cpt=0;
    int n;
    #pragma omp parallel for default(shared) private(n, fl0, Wl0, Hl0)  
    for(n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        //model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl0, fl0, 0, eta0, 0, asym, Wl0, 0, ratios_l0, step, trunc_c);
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl0, fl0, 0, 0, 0,0,0,0,0,asym, Wl0, 0, ratios_l0, step, trunc_c);
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }

    #pragma omp parallel for default(shared) private(n, fl1, Wl1, Hl1, a1)  
    for (n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        //model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl1, fl1, a1_l1[n], eta0, 0,asym, Wl1, 1, ratios_l1, step, trunc_c);
        double a1=a1_l1[n];
        double a2=0;
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl1, fl1, a1, a2, 0,0, 0,0,eta0,asym, Wl1, 1, ratios_l1, step, trunc_c);
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }

    #pragma omp parallel for default(shared) private(n, fl2, Wl2, Hl2, a1)  
    for(n=0; n<Nfl2; n++){ 
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
        //model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl2, fl2, a1_l2[n], eta0, a3,asym, Wl2, 2, ratios_l2, step, trunc_c);
        double a1=a1_l2[n]; //a1_terms[0] + a1_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        double a2=0;//a2_terms[0] + a2_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        //a3=a3_terms[0] + a3_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a4=0;//a4_terms[0] + a4_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl2, fl2, a1, a2, a3,a4, 0,0,eta0,asym, Wl2, 2, ratios_l2, step, trunc_c);        
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }     
    }

    #pragma omp parallel for default(shared) private(n, fl3, Wl3, Hl3, a1)  
    for(n=0; n<Nfl3; n++){ 
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
        //model_final=optimum_lorentzian_calc_a1etaa3(x, model_final, Hl3, fl3, a1_l3[n], eta0, a3, asym, Wl3, 3, ratios_l3, step, trunc_c);
        double a1=a1_l3[n]; //a1_terms[0] + a1_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        double a2=0;//a2_terms[0] + a2_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        //a3=a3_terms[0] + a3_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a4=0;//a4_terms[0] + a4_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a5=0;//a5_terms[0] + a5_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a6=0;//a6_terms[0] + a6_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl3, fl3, a1, a2, a3,a4, a5,a6,eta0,asym, Wl3, 3, ratios_l3, step, trunc_c);        
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
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
        
        file_out="params.raw";
        std::ofstream outfile;
        outfile.open(file_out.c_str());
            outfile << "# File generated using " << modelname << std::endl;
            outfile << "# This is the raw vector of (1) params_length and (2) params as provided in input of the model" << std::endl;
            outfile << params_length.transpose() << std::endl;
            outfile << params.transpose() << std::endl;
        outfile.close();
    }
    return model_final;
}


VectorXd model_RGB_asympt_aj_AppWidth_HarveyLike_v4_parallel(const VectorXd& params, const VectorXi& params_length, const VectorXd& x, bool outparams){
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
               The reference frequencies fref are NOT used. You basically fit the asymptotic. Note that the setups should have imposed ferr to fix in that scenario to avoid loss of efficiency while fitting
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
    VectorXd model_final(x.size());
    
    VectorXd fl0_all(Nmax), Wl0_all(Nmax), Hl0_all(Nmax), noise_params(Nnoise), fl1_all, Wl1_all, Hl1p_all, Hl1_all,a1_l1, a1_l2(Nfl2), a1_l3(Nfl3); //Hl0_all[Nmax],
    Data_eigensols freqs_l1;
    double fl0, fl1, fl2, fl3, Vl1, Vl2, Vl3, Hl0, Hl1, Hl2, Hl3, Wl0, Wl1, Wl2, Wl3;
    double a1,a2,a3,a4,a5,a6, eta0, asym;
    double numax, Htot, lnGamma0, lnLorentz;
    double e, tmp, r;
    int Nharvey;
    long cpt;
    // For outparams table
    const int Nrows=1000, Ncols=9; // Number of parameters for each mode
    MatrixXd mode_params(Nrows, Ncols); // For ascii outputs, if requested
    int k, L, Line=0; // will be used to trim the mode_params table where suited
    // For the mixed modes table
    VectorXd b, global_mixed_modes_params(10); // b is a bias vector
    bool included;
    MatrixXd mixed_modes_params;
    std::string mixed_modes_name_params;
    std::vector<std::string> global_mixed_modes_name_params(10);
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
    for (int i=0;i<Nferr; i++){
        fref_all[i]=params[Nmax+ lmax + Nfl0 + 8 + i];        
        ferr_all[i]=params[Nmax+ lmax + Nfl0 + 8 + Nferr + i];
    }
    // bias type setup. Note that you need at least 3 points to define a spline. Otherwise it will crash here.
    tk::spline bias; // declare a bias
    if (bias_type != 0){ // 0 is No bias case. Defined when priors are set to fix in the hyper priors
        bias.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0); // imposes second derivative to 0 at the edge (extrapolation if happens will be linear)
        if (bias_type == 1){ // CUBIC SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline);     // this calculates all spline coefficients   
        }
        if (bias_type == 2){ // HERMITE SPLINE
            bias.set_points(fref_all,ferr_all,tk::spline::cspline_hermite);     // this calculates all spline coefficients
        }
    }

    const double fmin=fl0_all.minCoeff(); // This is to avoid a change in the number of modes
    const double fmax=fl0_all.maxCoeff();  // This is to avoid a change in the number of modes 
    VectorXd xfit, rfit;
    xfit = Eigen::VectorXd::LinSpaced(fl0_all.size(), 0, fl0_all.size()-1);
    rfit=linfit(xfit, fl0_all); // rfit[0] = Dnu
    const double Dnu_p=rfit[0];
    const int n0=floor(rfit[1]/Dnu_p);
    const double epsilon_p = rfit[1]/Dnu_p -n0;
    if(fmin - Dnu_p <0){
        std::cerr << " ERROR: THE ARMM WILL NOT CONVERGE AS YOU HAVE fmin - Dnu_p < 0, leading the an INFINITE g modes density!" << std::endl;
        std::cerr << "           You need to consider a higher fmin and/or smaller Dnu_p..." << std::endl;
        std::cerr << "        fmin  : " << fmin << std::endl;
        std::cerr << "        Dnu_p : " << Dnu_p << std::endl;
        exit(EXIT_FAILURE);
    }
    if (model_type == 0){
        freqs_l1=solve_mm_asymptotic_O2p(Dnu_p, epsilon_p, 1, delta0l, 0, 0., DPl, alpha_g, q_star, 0, fmin - Dnu_p, fmax + Dnu_p, step, true, false);
    } else{
        freqs_l1=solve_mm_asymptotic_O2from_l0(fl0_all, 1, delta0l, DPl, alpha_g, q_star, 0, step, true, false, fmin, fmax);
    }
    fl1_all=freqs_l1.nu_m;

    // Adding the bias function to the vector fl1_all
    if (bias_type !=0){
        b.resize(fl1_all.size());
        b.setZero();
        for (int i=0; i< fl1_all.size(); i++){
            b[i]=bias(fl1_all[i]);
            fl1_all[i] = fl1_all[i] + b[i];
            //std::cout << bias(fl1_all[i]) << "  "  << fl1_all[i] << std::endl;
        }
    }
     // Generating widths profiles for l=1 modes using the ksi function
    ksi_pg=ksi_fct2(fl1_all, freqs_l1.nu_p, freqs_l1.nu_g, freqs_l1.dnup, freqs_l1.dPg, q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
    h1_h0_ratio=h_l_rgb(ksi_pg, Hfactor); //  Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019). Hence Hfactor to depart from asymptotic
   
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
    Wl1_all=gamma_l_fct2(ksi_pg, fl1_all, fl0_all, Wl0_all, h1_h0_ratio, 1, Wfactor); // generate the mixed modes widths // Added on 28 Apr: Wfactor: in front of the ksi fonction: recommended to have a alpha ~ N(1,sigma=0.1)
    // Generating splittings with a two-zone averaged rotation rates
    a1_l1=dnu_rot_2zones(ksi_pg, rot_env, rot_core);
    a1_l1=a1_l1.array().abs();
    a1_l2.setConstant(std::abs(rot_env));
    a1_l3.setConstant(std::abs(rot_env));
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
    int n;
    #pragma omp parallel for default(shared) private(n, fl0, Wl0, Hl0)  
    for(n=0; n<Nfl0; n++){
        fl0=fl0_all[n];
        Wl0=Wl0_all[n];
        Hl0=Hl0_all[n];
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl0, fl0, 0, 0, 0,0,0,0,0,asym, Wl0, 0, ratios_l0, step, trunc_c);
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 0, fl0, Hl0 , Wl0, 0, 0, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }    
    }
    #pragma omp parallel for default(shared) private(n, fl1, Wl1, Hl1, a1)   
    for (n=0; n<fl1_all.size(); n++){
        fl1=fl1_all[n];
        Wl1=Wl1_all[n];
        Hl1=std::abs(Hl1_all[n]);  
        double a1=a1_l1[n]; //a1_terms[0] + a1_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        double a2=0;//a2_terms[0] + a2_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl1, fl1, a1, a2, 0,0, 0,0,eta0,asym, Wl1, 1, ratios_l1, step, trunc_c);
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 1, fl1, Hl1 , Wl1, a1, eta0*1e-6, 0, asym, inclination;// mode_vec;
            Line=Line+1;
        }   
    }
    #pragma omp parallel for default(shared) private(n, fl2, Wl2, Hl2, a1)  
    for(n=0; n<Nfl2; n++){ 
        fl2=std::abs(params[Nmax+lmax+Nfl0+Nfl1+n]);
        lnGamma0=gamma_params[2] * log(fl2/gamma_params[0]) + log(gamma_params[3]);
        e=2.*log(fl2/gamma_params[1]) / log(gamma_params[4]/gamma_params[0]);
        lnLorentz=-log(gamma_params[5])/(1. + pow(e,2));     
        Wl2=exp(lnGamma0 + lnLorentz);
        if(do_amp){
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2/(pi*Wl2)*Vl2);
        } else{
            Hl2=lin_interpol(fl0_all, Hl0_all, fl2);
            Hl2=std::abs(Hl2*Vl2);
        }   
        double a1=a1_l2[n]; //a1_terms[0] + a1_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        double a2=0;//a2_terms[0] + a2_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        //a3=a3_terms[0] + a3_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a4=0;//a4_terms[0] + a4_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl2, fl2, a1, a2, a3,a4, 0,0,eta0,asym, Wl2, 2, ratios_l2, step, trunc_c);        
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 2, fl2, Hl2 , Wl2, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
            Line=Line+1;
        }     
    }
    #pragma omp parallel for default(shared) private(n, fl3, Wl3, Hl3, a1)  
    for(n=0; n<Nfl3; n++){ 
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
        double a1=a1_l3[n]; //a1_terms[0] + a1_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz  
        double a2=0;//a2_terms[0] + a2_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        //a3=a3_terms[0] + a3_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a4=0;//a4_terms[0] + a4_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a5=0;//a5_terms[0] + a5_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        double a6=0;//a6_terms[0] + a6_terms[1]*(fl1*1e-3); //two terms: one constant term + one linear in nu, after a11, set in mHz
        Optim_L model_tmp=optimum_lorentzian_calc_aj(x, Hl3, fl3, a1, a2, a3,a4, a5,a6,eta0,asym, Wl3, 3, ratios_l3, step, trunc_c);        
        #pragma omp critical
        model_final.segment(model_tmp.i0, model_tmp.N)=model_final.segment(model_tmp.i0, model_tmp.N) + model_tmp.y;
        if (outparams){
            #pragma omp critical
            mode_params.row(Line) << 3, fl3, Hl3 , Wl3, a1, eta0*1e-6, a3, asym, inclination;// mode_vec;
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
        // Prepare Mixed modes data arrays
        mixed_modes_params.resize(freqs_l1.nu_m.size(), 9);
        global_mixed_modes_name_params[0] = "model_type"; global_mixed_modes_name_params[1] = "Dnu_p"  ; global_mixed_modes_name_params[2] = "epsilon_p"; 
        global_mixed_modes_name_params[3] = "delta0l"   ; global_mixed_modes_name_params[4] = "alpha_p"; global_mixed_modes_name_params[5] = "nmax"; 
        global_mixed_modes_name_params[6] = "DPl"       ; global_mixed_modes_name_params[7] = "alpha_g"; global_mixed_modes_name_params[8] = "q_star"; 
        global_mixed_modes_name_params[9] = "sigma_p";
        if (model_type == 0){
            global_mixed_modes_params << model_type, Dnu_p , epsilon_p , delta0l , 0 , 0 , DPl , alpha_g , q_star , 0;
        } else{
            global_mixed_modes_params << model_type, -1 , -1 , delta0l , 0 , 0 , DPl , alpha_g , q_star , 0;            
        }
        mixed_modes_name_params = "# el / fl1_asymptotic / included? / spline_corr / ksi_pg / h1_h0 / Hl1p_all / Hl1_all / Wl1_all ";
        for(L=0;L<freqs_l1.nu_m.size();L++){
            k=0;
            included=0;
            while((fl1_all[k]- b[k] != freqs_l1.nu_m[L]) && k<fl1_all.size()-1){
                k=k+1;
            }
            if(fl1_all[k] - b[k] == freqs_l1.nu_m[L]){
                included=1;
            }
            if(included == 0){
                mixed_modes_params.row(L) << 1, freqs_l1.nu_m[L], included, -1 , -1 , -1 , -1, -1 , -1;
            } else{
                mixed_modes_params.row(L) << 1, freqs_l1.nu_m[L], included, b[k], ksi_pg[k], h1_h0_ratio[k], Hl1p_all[k], Hl1_all[k], Wl1_all[k];
            }   
        }
        // Writing on files all of the data 
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
        write_star_params_mixed(spec_params, params, params_length, mode_params.block(0,0, Line, mode_params.cols()), noise, 
                        file_out, modelname, name_params, mixed_modes_params, mixed_modes_name_params,
                        global_mixed_modes_params,  global_mixed_modes_name_params, freqs_l1.nu_p, freqs_l1.nu_g);
        file_out="params.raw";
        std::ofstream outfile;
        outfile.open(file_out.c_str());
            outfile << "# File generated using " << modelname << std::endl;
            outfile << "# This is the raw vector of (1) params_length and (2) params as provided in input of the model" << std::endl;
            outfile << params_length.transpose() << std::endl;
            outfile << params.transpose() << std::endl;
        outfile.close();
    }
    return model_final;
}


