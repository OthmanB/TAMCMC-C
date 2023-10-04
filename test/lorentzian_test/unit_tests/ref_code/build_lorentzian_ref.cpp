/*
 * build_lorentzian.cpp
 *
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
#include "build_lorentzian_ref.h"
#include <iostream>
#include <iomanip>
#include "../../external/Alm/Alm_cpp/activity.h"
#include "../../external/Alm/Alm_cpp/Alm_interpol.h"
#include "../../tamcmc/headers/acoefs.h"
using Eigen::VectorXd;
using Eigen::VectorXi;


VectorXd build_l_mode_a1l_etaa3_ref(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double f_s;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            if(l == 1){
                f_s=f_s1;
            }
            if(l == 2){
                f_s=f_s2;
            }
            if(l == 3){
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
            }
            profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + m*f_s + Pslm(3,l,m)*a3)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }
    
    return result;
}

VectorXd build_l_mode_a1l_a2a3_ref(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double f_s, a2_terms;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            a2_terms=Pslm(2,l,m)*a2;
            if(l == 1){
                //clm=0; // a3=0 for l=1 BY DEFINITION
                f_s=f_s1;
                //a2_terms=(3*m*m - 2)*a2;  // From Takashi note and Pnl decomposition: c2(n,l) = [3m*m - l(l+1)] / (2l-1)
            }
            if(l == 2){
                //clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                f_s=f_s2;
                //a2_terms=(m*m -2)*a2;
            }
            if(l == 3){
                //clm=(pow(m,3)-7*m)/2; // a3 implemented on 30/04/2021
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
                //a2_terms=(3*m*m - 12)*a2/5;
            }
            profile=(x_l - tmp.setConstant(fc_l + m*f_s + a2_terms + Pslm(3,l,m)*a3)).array().square(); // a1=f_s , a2 and a3 coefficients
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }
    
    return result;
}

VectorXd build_l_mode_a1etaa3_ref(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + m*f_s + Pslm(3,l,m)*a3)).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		if(asym == 0){ //Model with no asymetry
			result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
		} else{
			tmp.setConstant(1);
			asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
			result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
		}
	}

return result;
}

VectorXd build_l_mode_ajAlm_ref(const VectorXd& x_l, const double H_l, const double fc_l, const double a1, const double a3, const double a5, 
    const double eta0, const double epsilon_nl, const VectorXd& thetas, const double asym, const double gamma_l, const int l, 
    const VectorXd& V, const std::string filter_type, const gsl_funcs interp_funcs){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - acoefficients: a1, a3, a5
 *      - an Asphericity term eta (Centrifugal effect) + Alm term due to Active region following Gizon 2002, AN, 323, 251.
 *             Currently we use a hard-coded filter type "gate" which is rough but match the Gizon paper. "gauss" is also available and might be our final choice.
 *             Once we could compare the method adapted from Gizon works on global fits
 */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double nu_nlm;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            nu_nlm=fc_l + a1*Pslm(1,l,m) + a3*Pslm(3,l,m) + a5*Pslm(5,l,m); // Only even terms as odds terms are accounted by Alm
            if (eta0 > 0){
                nu_nlm = nu_nlm + fc_l*eta0*Qlm(l,m)*pow(a1*1e-6,2);
            } 
            if(filter_type == "gauss"){
                nu_nlm=nu_nlm + fc_l*epsilon_nl*Alm(l, m, thetas[0], thetas[1], filter_type); // Adding the Activity terms
            } else{
                nu_nlm=nu_nlm + fc_l*epsilon_nl*Alm_interp_iter_preinitialised(l, m, thetas[0], thetas[1], filter_type, interp_funcs); // Adding the Activity terms
            }
            profile=(x_l - tmp.setConstant(nu_nlm)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }
return result;
}

VectorXd build_l_mode_aj_ref(const VectorXd& x_l, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splittings in the form of a-coefficients aj with j={1,6}
 *      - eta: If eta0 > 0, account for the centrifugal distorsion in the a2 term (a2_CF) ==> This means the measured a2 = a2_AR : It WILL NOT include centrifugal effects
 *                         but the model DO account for it. In that case, eta0 =3./(4.*pi*rho*G) Should be the value given to that function: NOTE THAT Benomar2018 HAS A MISTAKE IN THAT FORMULATION (Eq. S6)
 *             If eta0 <=0 0, set a2_CF(eta) = 0 such that the measured a2 = a2_CF + a2_AR
*/
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double nu_nlm;

    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            nu_nlm=fc_l + a1*Pslm(1,l,m) + a2*Pslm(2,l,m) + a3*Pslm(3,l,m) + a4*Pslm(4,l,m)+ a5*Pslm(5,l,m) + a6*Pslm(6,l,m);
            if (eta0 > 0){
                nu_nlm = nu_nlm + fc_l*eta0*Qlm(l,m)*pow(a1*1e-6,2);
            } 
            profile=(x_l - tmp.setConstant(nu_nlm)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_l*V(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }

return result;
}

VectorXd build_l_mode_a1l_etaa3_v2_ref(const VectorXd& x_l, const VectorXd& H_lm, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     *		- Inclination IS NOT IMPOSED contrary to build_l_mode_a1l_etaa3. INSTEAD H_lm should be of dimension l(l+1) and provide all the heights 
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double f_s;
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            switch (l){
                case 1:
                    f_s=f_s1;
                    break;
                case 2:
                    f_s=f_s2;
                    break;
                case 3:
                    f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
                    break;
            }
            profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + Pslm(1,l,m)*f_s + Pslm(3,l,m)*a3)).array().square();
            //profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + Pslm(1,l,m)*a1 + clm*a3)).array().square();
            profile=4*profile/pow(gamma_l,2);
        } else{
            profile=(x_l - tmp.setConstant(fc_l)).array().square();
            profile=4*profile/pow(gamma_l,2);
        }
        if(asym == 0){ //Model with no asymetry
            result=result+ H_lm(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
        } else{
            tmp.setConstant(1);
            asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
            result=result+ H_lm(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
        }
    }
    
    std::cout << "NEED CHECKS in build_l_mode_a1l_etaa3_v2: The function was never verified" << std::endl;
    exit(EXIT_SUCCESS);

    return result;
}

VectorXd build_l_mode_a1etaa3_v2_ref(const VectorXd& x_l, const VectorXd& H_lm, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double clm;

	/*std::cout << " ---------- " << std::endl;
	std::cout << " l = " << l << std::endl;
	std::cout << "H_lm =" << H_lm << std::endl;
	std::cout << "fc_l =" << fc_l << std::endl;
	std::cout << "eta =" << eta << std::endl;
	std::cout << "a3 =" << a3 << std::endl;
	std::cout << "asym =" << asym << std::endl;
	std::cout << "gamma_l =" << gamma_l << std::endl;
	std::cout << " ---------- " << std::endl;
*/
	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta0*pow(f_s*1e-6,2)*Qlm(l,m)) + m*f_s + Pslm(3,l,m)*a3)).array().square();
			profile=4*profile/pow(gamma_l,2);
		} else{
			profile=(x_l - tmp.setConstant(fc_l)).array().square();
			profile=4*profile/pow(gamma_l,2);
		}
		if(asym == 0){ //Model with no asymetry
			result=result+ H_lm(m+l)* ((tmp.setConstant(1) + profile)).cwiseInverse();
		} else{
			tmp.setConstant(1);
			asymetry=(tmp + asym*(x_l/fc_l - tmp)).array().square() + (tmp2.setConstant(0.5*gamma_l*asym/fc_l)).array().square();
			result=result+ H_lm(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
		}
	}

return result;
}





VectorXd optimum_lorentzian_calc_a1l_a2a3_ref(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
    BEWARE: USES build_l_mode_a1a2a3() ==> LINEAR DEPENDENCE OF Asphericity IN NU IS NOT ACCOUNTED FOR... THIS DEPENDENCE MAY BE IMPLEMENTED AT HIGHER LEVEL when calling this function
     */
    double pmin, pmax, f_s;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
    switch(l){
        case 0:
            f_s=0.;
            break;
        case 1:
            f_s=f_s1;
            break;
        case 2:
            f_s=f_s2;
            break;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
            break;
    }
    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);

    m0=build_l_mode_a1l_a2a3_ref(x_l, H_l, fc_l, f_s1, f_s2, a2, a3, asym, gamma_l, l, V);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
    return y_out;
}


VectorXd optimum_lorentzian_calc_a1etaa3_ref(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
	double pmin, pmax;
	VectorXi ivals;
	VectorXd m0, x_l, y_out(y.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
 
	m0=build_l_mode_a1etaa3_ref(x_l, H_l, fc_l, f_s, eta0, a3, asym, gamma_l, l, V);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}

VectorXd optimum_lorentzian_calc_ajAlm_ref(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double a1, 
    const double a3, const double a5, const double eta0, const double epsilon_nl, const VectorXd& thetas, const double asym, 
    const double gamma_l, const int l, const VectorXd& V, const double step, const double c, const std::string filter_type, gsl_funcs interp_funcs){
/*
    function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
    that contains the lorentzian model.
*/
    double pmin, pmax;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, a1, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
 
    m0=build_l_mode_ajAlm_ref(x_l, H_l, fc_l, a1, a3, a5, eta0, epsilon_nl, thetas, asym, gamma_l, l, V, filter_type, interp_funcs);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}


VectorXd optimum_lorentzian_calc_aj_ref(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, 
        const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, 
        const double eta0, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
    function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
    that contains the lorentzian model.
    BEWARE: USES build_l_mode_aj() ==> LINEAR DEPENDENCE OF Asphericity IN NU IS NOT ACCOUNTED FOR... THIS DEPENDENCE MAY BE IMPLEMENTED AT HIGHER LEVEL when calling this function
*/
    double pmin, pmax;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(x.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, a1, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);

    m0=build_l_mode_aj_ref(x_l, H_l, fc_l, a1, a2, a3, a4, a5, a6, eta0, asym, gamma_l, l, V);
    y_out.segment(ivals[0], ivals[1]-ivals[0]) = y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}


VectorXd optimum_lorentzian_calc_a1etaa3_v2_ref(const VectorXd& x, const VectorXd& y, const VectorXd& H_lm, const double fc_l, const double f_s, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
	double pmin, pmax;
    VectorXi ivals;
	VectorXd m0, x_l, y_out(y.size());
    y_out=y;

    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);

	m0=build_l_mode_a1etaa3_v2_ref(x_l, H_lm, fc_l, f_s, eta0, a3, asym, gamma_l, l);
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
return y_out;
}

VectorXd optimum_lorentzian_calc_a1l_etaa3_ref(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
     BEWARE: USES build_l_mode_a1l_etaa3() ==> Asphericity is a linear term in nu
     */
    //const double c=20.;
    double pmin, pmax, f_s;
    VectorXi ivals;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
    switch(l){
        case 0:
            f_s=0.;
            break;
        case 1:
            f_s=f_s1;
            break;
        case 2:
            f_s=f_s2;
            break;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
            break;
    }
    
    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
    
    m0=build_l_mode_a1l_etaa3_ref(x_l, H_l, fc_l, f_s1, f_s2, eta0, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
    return y_out;
}

VectorXd optimum_lorentzian_calc_a1l_etaa3_v2_ref(const VectorXd& x, const VectorXd& y, const VectorXd& H_lm, const double fc_l, const double f_s1, const double f_s2, const double eta0, const double a3, const double asym, const double gamma_l, const int l, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
     BEWARE: USES build_l_mode_a1l_etaa3() ==> Asphericity is a linear term in nu

     This function differs from optimum_lorentzian_calc_a1l_etaa3 by the fact that it fits directly the (l,m) heights instead of considering H_l and V==ratios
     Thus an l=1 will have H_1m = [ H(m=-1), H(m=0), H(m=1)] components, etc...
     */
    double f_s;
    VectorXd m0, x_l, y_out(y.size());
    VectorXi ivals;
    y_out=y;;
    switch(l){
        case 0:
            f_s=0.;
            break;
        case 1:
            f_s=f_s1;
            break;
        case 2:
            f_s=f_s2;
            break;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
            break;
    }
    ivals=set_imin_imax(x, l, fc_l, gamma_l, f_s, c, step);    
    x_l=x.segment(ivals[0], ivals[1]-ivals[0]);
    m0=build_l_mode_a1l_etaa3_v2_ref(x_l, H_lm, fc_l, f_s1, f_s2, eta0, a3, asym, gamma_l, l);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(ivals[0], ivals[1]-ivals[0])= y_out.segment(ivals[0], ivals[1]-ivals[0]) + m0;
    return y_out;
}
