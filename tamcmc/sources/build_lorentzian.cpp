/*
 * build_lorentzian.cpp
 *
 *  Created on: 22 Feb 2016
 *      Author: obenomar
 */
#include <math.h>
#include <Eigen/Dense>
#include "build_lorentzian.h"
# include <iostream>
# include <iomanip>

using Eigen::VectorXd;

VectorXd build_l_mode_a1l_etaa3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double eta, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double Qlm, clm, f_s;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
            if(l == 1){
                clm=0; // a3=0 for l=1 BY DEFINITION
                f_s=f_s1;
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                f_s=f_s2;
            }
            if(l == 3){
                clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
            }
            profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
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

VectorXd build_l_mode_a1l_a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
    /*
     * This model IS WITHOUT THE ASSUMPTION S11=S22. It includes:
     *      - Asymetry of Lorentzian asym
     *      - splitting a1(l=1) and a1(l=2). ASSUMES a1(l=3) =  (a1(1) + a1(2))/2.
     *      - an Asphericity parameter eta
     *      - latitudinal effect a3(l=2) only. We consider a3(l=3)=0
     */
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double Qlm, clm, f_s, a2_terms;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            //Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
            if(l == 1){
                clm=0; // a3=0 for l=1 BY DEFINITION
                f_s=f_s1;
                a2_terms=(3*m*m - 2)*a2;  // From Takashi note and Pnl decomposition: c2(n,l) = [3m*m - l(l+1)] / (2l-1)
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                f_s=f_s2;
                a2_terms=(m*m -2)*a2;
            }
            if(l == 3){
                clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
                a2_terms=(3*m*m - 12)*a2/5;
            }
            //profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
            profile=(x_l - tmp.setConstant(fc_l + m*f_s + a2_terms + clm*a3)).array().square(); // a1=f_s , a2 and a3 coefficients
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

VectorXd build_l_mode_a1etaa3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double Qlm, clm;

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
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

VectorXd build_l_mode_a1etaGlma3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, 
    const double eta0, const double epsilon_nl, const VectorXd& thetas, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity term eta (Centrifugal effect) and due to Active region following Gizon 2002, AN, 323, 251.
 *             Additional assumption are made in order to have a model that contain both a pole and equatorial activity band
 *             using shift(Equator) = epsilon*alpha_nl(Equator) and shift(Pole) = epsilon*(1-alpha_nl). Mutual exclusions conditions
 *             might need to be used (see prior_application.cpp) in order to avoid degeneracies
 *      - latitudinal effect a3
 */
    const long Nxl=x_l.size();
    const long double Dnl=0.75;
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double Qlm, clm, Glm;

    result.setZero();
    for(int m=-l; m<=l; m++){
        Glm=Glm(l, m, thetas[0], thetas[1]); // WARNING WARNING WARNING: HERE WE APPLY GLM ONLY TO l>0 BUT MAY BE NOT CORRECT. NOTE: l=0,m=0 might only be shifted
        if(l != 0){
            Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
            if(l == 1){
                clm=m; // a3 for l=1
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
            }
            if(l == 3){
                clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
            }
            CF_term=eta0*Dnl*f_s*f_s*Qlm;
            AR_term=epsilon_nl*Glm;
            profile=(x_l - tmp.setConstant(fc_l*(1. + CF_term + AR_term) + m*f_s + clm*a3)).array().square();
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

VectorXd build_l_mode_a1a2a3(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
    const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
    double Qlm, clm, a2_terms;

    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            //Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
            if(l == 1){
                clm=m; // a3 for l=1
                a2_terms=(3*m*m - 2)*a2;  // From Takashi note and Pnl decomposition: c2(n,l) = [3m*m - l(l+1)] / (2l-1)
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                a2_terms=(m*m -2)*a2;
            }
            if(l == 3){
                clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
                a2_terms=(3*m*m - 12)*a2/5;
            }
            profile=(x_l - tmp.setConstant(fc_l + m*f_s + a2_terms + clm*a3)).array().square();
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


VectorXd build_l_mode_a1l_etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm, const double fc_l, const double f_s1, const double f_s2, const double eta, const double a3, const double asym, const double gamma_l, const int l){
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
    double Qlm, clm, f_s;
    
    result.setZero();
    for(int m=-l; m<=l; m++){
        if(l != 0){
            Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
            if(l == 1){
                clm=0; // a3=0 for l=1 BY DEFINITION
                f_s=f_s1;
            }
            if(l == 2){
                clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
                f_s=f_s2;
            }
            if(l == 3){
                clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
                f_s=(f_s1 + f_s2)/2.; // APPROXIMATION
            }
            profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
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

VectorXd build_l_mode_a1etaa3_v2(const VectorXd& x_l, const VectorXd& H_lm, const double fc_l, const double f_s, const double eta, const double a3, const double asym, const double gamma_l, const int l){
/*
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - an Asphericity parameter eta
 *      - latitudinal effect a3
*/
	const long Nxl=x_l.size();
    VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double Qlm, clm;

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
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
			profile=(x_l - tmp.setConstant(fc_l*(1. + eta*Qlm) + m*f_s + clm*a3)).array().square();
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

VectorXd build_l_mode_asym_act(const VectorXd& x_l, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, const double b, const double alpha, const double asym, const double gamma_l, const int l, const VectorXd& V){
/*
 * ---- OBSELETE FUNCTION ----
 * This model includes:
 *      - Asymetry of Lorentzian asym
 *      - splitting a1
 *      - centrifugal force effect eta.nu.a1^2 (rotation-induced oblateness)
 *      - latitudinal effect a3
 *      - effect of magnetic field of the form b.nu^alpha
 * ---------------------------
*/
	const long Nxl=x_l.size();
        VectorXd profile(Nxl), tmp(Nxl), tmp2(Nxl), result(Nxl), asymetry(Nxl);
	double Qlm, clm;

	result.setZero();
	for(int m=-l; m<=l; m++){
		if(l != 0){
			Qlm=(l*(l+1) - 3*pow(m,2))/((2*l - 1)*(2*l + 3)); // accounting for eta
			if(l == 1){
				clm=m; // a3 for l=1
			}
			if(l == 2){
				clm=(5*pow(m,3) - 17*m)/3.; // a3 for l=2
			}
			if(l == 3){
				clm=0; // a3 NOT YET IMPLEMENTED FOR l=3
			}
		
			profile=(x_l - tmp.setConstant(fc_l + Qlm*(eta*fc_l*pow(f_s,2) + b*pow(fc_l*1e-3,alpha)) + m*f_s + clm*a3) ).array().square();
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
			//std::cout << "asymetry=" << asymetry.transpose() <<std::endl;
			//tmp2=((tmp.setConstant(1) + profile)).cwiseInverse();
			result=result+ H_l*V(m+l)*asymetry.cwiseProduct(((tmp.setConstant(1) + profile)).cwiseInverse());
			//std::cout << "in build_lorentzian" <<std::endl;
		}
	}

//std::cout << "result=" << result.transpose() <<std::endl;
return result;
}

VectorXd optimum_lorentzian_calc_a1l_etaa3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double eta, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
     BEWARE: USES build_l_mode_a1l_etaa3() ==> Asphericity is a linear term in nu
     */
    //const double c=20.;
    double pmin, pmax, f_s;
    long imin, imax;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
    //VectorXd mall(x.size());
    switch(l){
        case 0:
            f_s=0.;
        case 1:
            f_s=f_s1;
        case 2:
            f_s=f_s2;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
    }
    if(gamma_l >= 1 && f_s >= 1){
        if(l != 0){
            pmin=fc_l - c*(l*f_s + gamma_l);
            pmax=fc_l + c*(l*f_s + gamma_l);
        } else{
            pmin=fc_l - c*gamma_l*2.2;
            pmax=fc_l + c*gamma_l*2.2;
        }
    }
    if(gamma_l <= 1 && f_s >= 1){
        if(l !=0){
            pmin=fc_l -c*(l*f_s + 1);
            pmax=fc_l + c*(l*f_s + 1);
        } else{
            pmin=fc_l - c*2.2;
            pmax=fc_l + c*2.2;
        }
    }
    if(gamma_l >= 1 && f_s <= 1){
        if(l != 0){
            pmin=fc_l - c*(l + gamma_l);
            pmax=fc_l + c*(l + gamma_l);
        } else{
            pmin=fc_l - c*2.2*gamma_l;
            pmax=fc_l + c*2.2*gamma_l;
        }
    }
    if(gamma_l <= 1 && f_s <= 1){
        if(l !=0){
            pmin=fc_l - c*(l+1);
            pmax=fc_l + c*(l+1);
        } else{
            pmin=fc_l -c*2.2;
            pmax=fc_l +c*2.2;
        }
    }
    // ---------- Handling boundaries ----------
    //     case when the proposed values lead to (pmax - step) < min(x)....
    if( (pmax - step) < x.head(1)(0)){
        pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
    }
    //     case when the proposed values lead to (pmin + step) >= max(x)...
    if( (pmin + step) >= x.tail(1)(0)){
        pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
    }
    
    imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
    imax=ceil((pmax-x.head(1)(0))/step);
    
    if(imin < 0){imin=0;} // ensure that we always stay within the vector x
    if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
    if(imax-imin <= 0){
        std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
        std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
        std::cout << " - step=" << step << std::endl;
        std::cout << " --------" << std::endl;
        std::cout << " - l=" << l << std::endl;
        std::cout << " - fc_l=" << fc_l << std::endl;
        std::cout << " - H_l=" << H_l << std::endl;
        std::cout << " - gamma_l=" << gamma_l << std::endl;
        std::cout << " - f_s=" << f_s << std::endl;
        std::cout << " - eta=" << eta << std::endl;
        std::cout << " - a3=" << a3 << std::endl;
        std::cout << " - asym=" << asym << std::endl;
        std::cout << " --------" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
    //std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
    
    //std::cout << "step=" << step << std::endl;
    //std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
    //std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();
    
    x_l=x.segment(imin, imax-imin);
    
    m0=build_l_mode_a1l_etaa3(x_l, H_l, fc_l, f_s1, f_s2, eta, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
    return y_out;
}


VectorXd optimum_lorentzian_calc_a1l_a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s1, const double f_s2, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
    BEWARE: USES build_l_mode_a1a2a3() ==> LINEAR DEPENDENCE OF Asphericity IN NU IS NOT ACCOUNTED FOR... THIS DEPENDENCE MAY BE IMPLEMENTED AT HIGHER LEVEL when calling this function
     */
    //const double c=20.;
    double pmin, pmax, f_s;
    long imin, imax;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
    //VectorXd mall(x.size());
    switch(l){
        case 0:
            f_s=0.;
        case 1:
            f_s=f_s1;
        case 2:
            f_s=f_s2;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
    }
    if(gamma_l >= 1 && f_s >= 1){
        if(l != 0){
            pmin=fc_l - c*(l*f_s + gamma_l);
            pmax=fc_l + c*(l*f_s + gamma_l);
        } else{
            pmin=fc_l - c*gamma_l*2.2;
            pmax=fc_l + c*gamma_l*2.2;
        }
    }
    if(gamma_l <= 1 && f_s >= 1){
        if(l !=0){
            pmin=fc_l -c*(l*f_s + 1);
            pmax=fc_l + c*(l*f_s + 1);
        } else{
            pmin=fc_l - c*2.2;
            pmax=fc_l + c*2.2;
        }
    }
    if(gamma_l >= 1 && f_s <= 1){
        if(l != 0){
            pmin=fc_l - c*(l + gamma_l);
            pmax=fc_l + c*(l + gamma_l);
        } else{
            pmin=fc_l - c*2.2*gamma_l;
            pmax=fc_l + c*2.2*gamma_l;
        }
    }
    if(gamma_l <= 1 && f_s <= 1){
        if(l !=0){
            pmin=fc_l - c*(l+1);
            pmax=fc_l + c*(l+1);
        } else{
            pmin=fc_l -c*2.2;
            pmax=fc_l +c*2.2;
        }
    }
    // ---------- Handling boundaries ----------
    //     case when the proposed values lead to (pmax - step) < min(x)....
    if( (pmax - step) < x.head(1)(0)){
        pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
    }
    //     case when the proposed values lead to (pmin + step) >= max(x)...
    if( (pmin + step) >= x.tail(1)(0)){
        pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
    }
    
    imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
    imax=ceil((pmax-x.head(1)(0))/step);
    
    if(imin < 0){imin=0;} // ensure that we always stay within the vector x
    if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
    if(imax-imin <= 0){
        std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
        std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
        std::cout << " - step=" << step << std::endl;
        std::cout << " --------" << std::endl;
        std::cout << " - l=" << l << std::endl;
        std::cout << " - fc_l=" << fc_l << std::endl;
        std::cout << " - H_l=" << H_l << std::endl;
        std::cout << " - gamma_l=" << gamma_l << std::endl;
        std::cout << " - f_s=" << f_s << std::endl;
        std::cout << " - a2=" << a2 << std::endl;
        std::cout << " - a3=" << a3 << std::endl;
        std::cout << " - asym=" << asym << std::endl;
        std::cout << " --------" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
    //std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
    
    //std::cout << "step=" << step << std::endl;
    //std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
    //std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();
    
    x_l=x.segment(imin, imax-imin);
    
    m0=build_l_mode_a1l_a2a3(x_l, H_l, fc_l, f_s1, f_s2, a2, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
    return y_out;
}


VectorXd optimum_lorentzian_calc_a1etaa3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
	//const double c=20.;
	double pmin, pmax;
	long imin, imax;
	VectorXd m0, x_l, y_out(y.size());
    y_out=y;
	//VectorXd mall(x.size());

	if(gamma_l >= 1 && f_s >= 1){
		if(l != 0){
			pmin=fc_l - c*(l*f_s + gamma_l);
			pmax=fc_l + c*(l*f_s + gamma_l);
		} else{
			pmin=fc_l - c*gamma_l*2.2;
			pmax=fc_l + c*gamma_l*2.2;
		}
	}
	if(gamma_l <= 1 && f_s >= 1){
		if(l !=0){
			pmin=fc_l -c*(l*f_s + 1);
			pmax=fc_l + c*(l*f_s + 1);
		} else{
			pmin=fc_l - c*2.2;
			pmax=fc_l + c*2.2;
		}
	}
	if(gamma_l >= 1 && f_s <= 1){
		if(l != 0){
			pmin=fc_l - c*(l + gamma_l);
			pmax=fc_l + c*(l + gamma_l);
		} else{
			pmin=fc_l - c*2.2*gamma_l;
			pmax=fc_l + c*2.2*gamma_l;
		}
	}
	if(gamma_l <= 1 && f_s <= 1){
		if(l !=0){
			pmin=fc_l - c*(l+1);
			pmax=fc_l + c*(l+1);
		} else{
			pmin=fc_l -c*2.2;
			pmax=fc_l +c*2.2;
		}
	}
	// ---------- Handling boundaries ----------
	//     case when the proposed values lead to (pmax - step) < min(x)....
	if( (pmax - step) < x.head(1)(0)){ 
		pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
	}
	//     case when the proposed values lead to (pmin + step) >= max(x)...
	if( (pmin + step) >= x.tail(1)(0)){
		pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
	}

	imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
	imax=ceil((pmax-x.head(1)(0))/step);

	if(imin < 0){imin=0;} // ensure that we always stay within the vector x
	if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
	if(imax-imin <= 0){
		std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
		std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
		std::cout << " - step=" << step << std::endl;
		std::cout << " --------" << std::endl;
		std::cout << " - l=" << l << std::endl;
		std::cout << " - fc_l=" << fc_l << std::endl;
		std::cout << " - H_l=" << H_l << std::endl;
		std::cout << " - gamma_l=" << gamma_l << std::endl;
		std::cout << " - f_s=" << f_s << std::endl;
		std::cout << " - eta=" << eta << std::endl;
		std::cout << " - a3=" << a3 << std::endl;
		std::cout << " - asym=" << asym << std::endl;
		std::cout << " --------" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
	//std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
	
	//std::cout << "step=" << step << std::endl;
	//std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
	//std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

	x_l=x.segment(imin, imax-imin);
 
	m0=build_l_mode_a1etaa3(x_l, H_l, fc_l, f_s, eta, a3, asym, gamma_l, l, V);
	//mall.setZero();
	//mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
return y_out;
}

VectorXd optimum_lorentzian_calc_a1etaGlma3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double eta0, const VectorXd& thetas, const double epsilon_nl, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
    function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
    that contains the lorentzian model.
    BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
    //const double c=20.;
    double pmin, pmax;
    long imin, imax;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;
    //VectorXd mall(x.size());

    if(gamma_l >= 1 && f_s >= 1){
        if(l != 0){
            pmin=fc_l - c*(l*f_s + gamma_l);
            pmax=fc_l + c*(l*f_s + gamma_l);
        } else{
            pmin=fc_l - c*gamma_l*2.2;
            pmax=fc_l + c*gamma_l*2.2;
        }
    }
    if(gamma_l <= 1 && f_s >= 1){
        if(l !=0){
            pmin=fc_l -c*(l*f_s + 1);
            pmax=fc_l + c*(l*f_s + 1);
        } else{
            pmin=fc_l - c*2.2;
            pmax=fc_l + c*2.2;
        }
    }
    if(gamma_l >= 1 && f_s <= 1){
        if(l != 0){
            pmin=fc_l - c*(l + gamma_l);
            pmax=fc_l + c*(l + gamma_l);
        } else{
            pmin=fc_l - c*2.2*gamma_l;
            pmax=fc_l + c*2.2*gamma_l;
        }
    }
    if(gamma_l <= 1 && f_s <= 1){
        if(l !=0){
            pmin=fc_l - c*(l+1);
            pmax=fc_l + c*(l+1);
        } else{
            pmin=fc_l -c*2.2;
            pmax=fc_l +c*2.2;
        }
    }
    // ---------- Handling boundaries ----------
    //     case when the proposed values lead to (pmax - step) < min(x)....
    if( (pmax - step) < x.head(1)(0)){ 
        pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
    }
    //     case when the proposed values lead to (pmin + step) >= max(x)...
    if( (pmin + step) >= x.tail(1)(0)){
        pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
    }

    imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
    imax=ceil((pmax-x.head(1)(0))/step);

    if(imin < 0){imin=0;} // ensure that we always stay within the vector x
    if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
    if(imax-imin <= 0){
        std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
        std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
        std::cout << " - step=" << step << std::endl;
        std::cout << " --------" << std::endl;
        std::cout << " - l=" << l << std::endl;
        std::cout << " - fc_l=" << fc_l << std::endl;
        std::cout << " - H_l=" << H_l << std::endl;
        std::cout << " - gamma_l=" << gamma_l << std::endl;
        std::cout << " - f_s=" << f_s << std::endl;
        std::cout << " - eta=" << eta << std::endl;
        std::cout << " - a3=" << a3 << std::endl;
        std::cout << " - asym=" << asym << std::endl;
        std::cout << " --------" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
    //std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
    
    //std::cout << "step=" << step << std::endl;
    //std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
    //std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

    x_l=x.segment(imin, imax-imin);
 
    m0=build_l_mode_a1etaGlma3(x_l, H_l, fc_l, f_s, eta0, epsilon_nl, thetas, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
return y_out;
}


VectorXd optimum_lorentzian_calc_a1a2a3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double a2, const double a3, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
    function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
    that contains the lorentzian model.
    BEWARE: USES build_l_mode_a1a2a3() ==> LINEAR DEPENDENCE OF Asphericity IN NU IS NOT ACCOUNTED FOR... THIS DEPENDENCE MAY BE IMPLEMENTED AT HIGHER LEVEL when calling this function
*/
    //const double c=20.;
    double pmin, pmax;
    long imin, imax;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;
    //VectorXd mall(x.size());

    if(gamma_l >= 1 && f_s >= 1){
        if(l != 0){
            pmin=fc_l - c*(l*f_s + gamma_l);
            pmax=fc_l + c*(l*f_s + gamma_l);
        } else{
            pmin=fc_l - c*gamma_l*2.2;
            pmax=fc_l + c*gamma_l*2.2;
        }
    }
    if(gamma_l <= 1 && f_s >= 1){
        if(l !=0){
            pmin=fc_l -c*(l*f_s + 1);
            pmax=fc_l + c*(l*f_s + 1);
        } else{
            pmin=fc_l - c*2.2;
            pmax=fc_l + c*2.2;
        }
    }
    if(gamma_l >= 1 && f_s <= 1){
        if(l != 0){
            pmin=fc_l - c*(l + gamma_l);
            pmax=fc_l + c*(l + gamma_l);
        } else{
            pmin=fc_l - c*2.2*gamma_l;
            pmax=fc_l + c*2.2*gamma_l;
        }
    }
    if(gamma_l <= 1 && f_s <= 1){
        if(l !=0){
            pmin=fc_l - c*(l+1);
            pmax=fc_l + c*(l+1);
        } else{
            pmin=fc_l -c*2.2;
            pmax=fc_l +c*2.2;
        }
    }
    // ---------- Handling boundaries ----------
    //     case when the proposed values lead to (pmax - step) < min(x)....
    if( (pmax - step) < x.head(1)(0)){ 
        pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
    }
    //     case when the proposed values lead to (pmin + step) >= max(x)...
    if( (pmin + step) >= x.tail(1)(0)){
        pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
    }

    imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
    imax=ceil((pmax-x.head(1)(0))/step);

    if(imin < 0){imin=0;} // ensure that we always stay within the vector x
    if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
    if(imax-imin <= 0){
        std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
        std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
        std::cout << " - step=" << step << std::endl;
        std::cout << " --------" << std::endl;
        std::cout << " - l=" << l << std::endl;
        std::cout << " - fc_l=" << fc_l << std::endl;
        std::cout << " - H_l=" << H_l << std::endl;
        std::cout << " - gamma_l=" << gamma_l << std::endl;
        std::cout << " - f_s=" << f_s << std::endl;
        std::cout << " - a2=" << a2 << std::endl;
        std::cout << " - a3=" << a3 << std::endl;
        std::cout << " - asym=" << asym << std::endl;
        std::cout << " --------" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
    //std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
    
    //std::cout << "step=" << step << std::endl;
    //std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
    //std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

    x_l=x.segment(imin, imax-imin);
 
    m0=build_l_mode_a1a2a3(x_l, H_l, fc_l, f_s, a2, a3, asym, gamma_l, l, V);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
return y_out;
}

VectorXd optimum_lorentzian_calc_a1etaa3_v2(const VectorXd& x, const VectorXd& y, const VectorXd& H_lm, const double fc_l, const double f_s, const double eta, const double a3, const double asym, const double gamma_l, const int l, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_a1etaa3() ==> Asphericity is a linear term in nu
*/
	//const double c=20.;
	double pmin, pmax;
	long imin, imax;
	VectorXd m0, x_l, y_out(y.size());
    y_out=y;
	//VectorXd mall(x.size());

	if(gamma_l >= 1 && f_s >= 1){
		if(l != 0){
			pmin=fc_l - c*(l*f_s + gamma_l);
			pmax=fc_l + c*(l*f_s + gamma_l);
		} else{
			pmin=fc_l - c*gamma_l*2.2;
			pmax=fc_l + c*gamma_l*2.2;
		}
	}
	if(gamma_l <= 1 && f_s >= 1){
		if(l !=0){
			pmin=fc_l -c*(l*f_s + 1);
			pmax=fc_l + c*(l*f_s + 1);
		} else{
			pmin=fc_l - c*2.2;
			pmax=fc_l + c*2.2;
		}
	}
	if(gamma_l >= 1 && f_s <= 1){
		if(l != 0){
			pmin=fc_l - c*(l + gamma_l);
			pmax=fc_l + c*(l + gamma_l);
		} else{
			pmin=fc_l - c*2.2*gamma_l;
			pmax=fc_l + c*2.2*gamma_l;
		}
	}
	if(gamma_l <= 1 && f_s <= 1){
		if(l !=0){
			pmin=fc_l - c*(l+1);
			pmax=fc_l + c*(l+1);
		} else{
			pmin=fc_l -c*2.2;
			pmax=fc_l +c*2.2;
		}
	}
	// ---------- Handling boundaries ----------
	//     case when the proposed values lead to (pmax - step) < min(x)....
	if( (pmax - step) < x.head(1)(0)){ 
		pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
	}
	//     case when the proposed values lead to (pmin + step) >= max(x)...
	if( (pmin + step) >= x.tail(1)(0)){
		pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
	}

	imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
	imax=ceil((pmax-x.head(1)(0))/step);

	if(imin < 0){imin=0;} // ensure that we always stay within the vector x
	if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
	if(imax-imin <= 0){
		std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
		std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
		std::cout << " - step=" << step << std::endl;
		std::cout << " --------" << std::endl;
		std::cout << " - l=" << l << std::endl;
		std::cout << " - fc_l=" << fc_l << std::endl;
		std::cout << " - H_lm=" << H_lm << std::endl;
		std::cout << " - gamma_l=" << gamma_l << std::endl;
		std::cout << " - f_s=" << f_s << std::endl;
		std::cout << " - eta=" << eta << std::endl;
		std::cout << " - a3=" << a3 << std::endl;
		std::cout << " - asym=" << asym << std::endl;
		std::cout << " --------" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
	//std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
	
	//std::cout << "step=" << step << std::endl;
	//std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
	//std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

	x_l=x.segment(imin, imax-imin);
 
	m0=build_l_mode_a1etaa3_v2(x_l, H_lm, fc_l, f_s, eta, a3, asym, gamma_l, l);
	//mall.setZero();
	//mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
return y_out;
}


VectorXd optimum_lorentzian_calc_a1acta3(const VectorXd& x, const VectorXd& y, const double H_l, const double fc_l, const double f_s, const double eta, const double a3, 
		const double b, const double alpha, const double asym, const double gamma_l, const int l, const VectorXd& V, const double step, const double c){
/*
	function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
	that contains the lorentzian model.
	BEWARE: USES build_l_mode_asym_act() ==> Include mode asymetry and effect of activity

*/
	//const double c=20.;
	double pmin, pmax;
	long imin, imax;
	VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
	//VectorXd mall(x.size());

	if(gamma_l >= 1 && f_s >= 1){
		if(l != 0){
			pmin=fc_l - c*(l*f_s + gamma_l);
			pmax=fc_l + c*(l*f_s + gamma_l);
		} else{
			pmin=fc_l - c*gamma_l*2.2;
			pmax=fc_l + c*gamma_l*2.2;
		}
	}
	if(gamma_l <= 1 && f_s >= 1){
		if(l !=0){
			pmin=fc_l -c*(l*f_s + 1);
			pmax=fc_l + c*(l*f_s + 1);
		} else{
			pmin=fc_l - c*2.2;
			pmax=fc_l + c*2.2;
		}
	}
	if(gamma_l >= 1 && f_s <= 1){
		if(l != 0){
			pmin=fc_l - c*(l + gamma_l);
			pmax=fc_l + c*(l + gamma_l);
		} else{
			pmin=fc_l - c*2.2*gamma_l;
			pmax=fc_l + c*2.2*gamma_l;
		}
	}
	if(gamma_l <= 1 && f_s <= 1){
		if(l !=0){
			pmin=fc_l - c*(l+1);
			pmax=fc_l + c*(l+1);
		} else{
			pmin=fc_l -c*2.2;
			pmax=fc_l +c*2.2;
		}
	}
	// ---------- Handling boundaries ----------
	//     case when the proposed values lead to (pmax - step) < min(x)....
	if( (pmax - step) < x.head(1)(0)){ 
		pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
	}
	//     case when the proposed values lead to (pmin + step) >= max(x)...
	if( (pmin + step) >= x.tail(1)(0)){
		pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
	}

	imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
	imax=ceil((pmax-x.head(1)(0))/step);

	if(imin < 0){imin=0;} // ensure that we always stay within the vector x
	if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
	if(imax-imin <= 0){
		std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
		std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
		std::cout << " - step=" << step << std::endl;
		std::cout << " --------" << std::endl;
		std::cout << " - l=" << l << std::endl;
		std::cout << " - fc_l=" << fc_l << std::endl;
		std::cout << " - H_l=" << H_l << std::endl;
		std::cout << " - gamma_l=" << gamma_l << std::endl;
		std::cout << " - f_s=" << f_s << std::endl;
		std::cout << " - eta=" << eta << std::endl;
		std::cout << " - a3=" << a3 << std::endl;
		std::cout << " - b=" << b << std::endl;
		std::cout << " - alpha=" << alpha << std::endl;
		std::cout << " - asym=" << asym << std::endl;
		std::cout << " --------" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
	//std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
	
	//std::cout << "step=" << step << std::endl;
	//std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
	//std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();

	x_l=x.segment(imin, imax-imin);
 
	//m0=build_l_mode(x_l, H_l, fc_l, f_s, a2, a3, gamma_l, l, V);
	m0=build_l_mode_asym_act(x_l, H_l, fc_l, f_s, eta, a3, b, alpha, asym, gamma_l, l, V);
	//mall.setZero();
	//mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
return y_out;
}

VectorXd optimum_lorentzian_calc_a1l_etaa3_v2(const VectorXd& x, const VectorXd& y, const VectorXd& H_lm, const double fc_l, const double f_s1, const double f_s2, const double eta, const double a3, const double asym, const double gamma_l, const int l, const double step, const double c){
    /*
     function that calculates the lorentzian on a optimized range of frequency. It returns a Vector of same size as the original vector x
     that contains the lorentzian model.
     BEWARE: USES build_l_mode_a1l_etaa3() ==> Asphericity is a linear term in nu

     This function differs from optimum_lorentzian_calc_a1l_etaa3 by the fact that it fits directly the (l,m) heights instead of considering H_l and V==ratios
     Thus an l=1 will have H_1m = [ H(m=-1), H(m=0), H(m=1)] components, etc...
     */
    //const double c=20.;
    double pmin, pmax, f_s;
    long imin, imax;
    VectorXd m0, x_l, y_out(y.size());
    y_out=y;;
    //VectorXd mall(x.size());
    switch(l){
        case 0:
            f_s=0.;
        case 1:
            f_s=f_s1;
        case 2:
            f_s=f_s2;
        case 3:
            f_s=(f_s1 + f_s2)/2.;
    }
    if(gamma_l >= 1 && f_s >= 1){
        if(l != 0){
            pmin=fc_l - c*(l*f_s + gamma_l);
            pmax=fc_l + c*(l*f_s + gamma_l);
        } else{
            pmin=fc_l - c*gamma_l*2.2;
            pmax=fc_l + c*gamma_l*2.2;
        }
    }
    if(gamma_l <= 1 && f_s >= 1){
        if(l !=0){
            pmin=fc_l -c*(l*f_s + 1);
            pmax=fc_l + c*(l*f_s + 1);
        } else{
            pmin=fc_l - c*2.2;
            pmax=fc_l + c*2.2;
        }
    }
    if(gamma_l >= 1 && f_s <= 1){
        if(l != 0){
            pmin=fc_l - c*(l + gamma_l);
            pmax=fc_l + c*(l + gamma_l);
        } else{
            pmin=fc_l - c*2.2*gamma_l;
            pmax=fc_l + c*2.2*gamma_l;
        }
    }
    if(gamma_l <= 1 && f_s <= 1){
        if(l !=0){
            pmin=fc_l - c*(l+1);
            pmax=fc_l + c*(l+1);
        } else{
            pmin=fc_l -c*2.2;
            pmax=fc_l +c*2.2;
        }
    }
    // ---------- Handling boundaries ----------
    //     case when the proposed values lead to (pmax - step) < min(x)....
    if( (pmax - step) < x.head(1)(0)){
        pmax=x.head(1)(0) +c;  // bundary = first value of x plus c
    }
    //     case when the proposed values lead to (pmin + step) >= max(x)...
    if( (pmin + step) >= x.tail(1)(0)){
        pmin=x.tail(1)(0) - c; // bundary = last value of x minus c
    }
    
    imin=floor((pmin-x.head(1)(0))/step); // Here it is assumed that there is a regular grid WITHOUT WHOLES (CANNOT USE a .FILT file)
    imax=ceil((pmax-x.head(1)(0))/step);
    
    if(imin < 0){imin=0;} // ensure that we always stay within the vector x
    if(imax > x.size()){imax=x.size();} // ensure that we always stay within the vector x
    if(imax-imin <= 0){
        std::cout << "Warning imax -imin <= 0 : imin=" << imin << "   imax=" << imax << std::endl;
        std::cout << " - pmin=" << pmin << "   pmax=" << pmax << std::endl;
        std::cout << " - step=" << step << std::endl;
        std::cout << " --------" << std::endl;
        std::cout << " - l=" << l << std::endl;
        std::cout << " - fc_l=" << fc_l << std::endl;
        std::cout << " - H_lm=" << H_lm << std::endl;
        std::cout << " - gamma_l=" << gamma_l << std::endl;
        std::cout << " - f_s=" << f_s << std::endl;
        std::cout << " - eta=" << eta << std::endl;
        std::cout << " - a3=" << a3 << std::endl;
        std::cout << " - asym=" << asym << std::endl;
        std::cout << " --------" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //std::cout << "xmin=" << x.head(1) << "(microHz)" << std::endl;
    //std::cout << "xmax=" << x.tail(1) << "(microHz)" << std::endl;
    
    //std::cout << "step=" << step << std::endl;
    //std::cout << "pmin=" << pmin << "(microHz)  pmax=" << pmax << "(microHz)" << std::endl;
    //std::cout << "imin=" << imin << "  imax=" << imax << std::endl;
    //std::cin.ignore();
    
    x_l=x.segment(imin, imax-imin);
    
    m0=build_l_mode_a1l_etaa3_v2(x_l, H_lm, fc_l, f_s1, f_s2, eta, a3, asym, gamma_l, l);
    //mall.setZero();
    //mall.segment(imin, imax-imin)=m0;
    y_out.segment(imin, imax-imin)= y_out.segment(imin, imax-imin) + m0;
    return y_out;
}
