// -------------------
// ---- Functions used to compute models based on the asymptotics ----
// This contains all the function that describes the Bumped period spacing
// and how they can be used to derived mode amplitudes from inertia ratio
// as they have been developed and tested. This arises from the following publications:
// https://arxiv.org/pdf/1509.06193.pdf (Inertia and ksi relation)
// https://www.aanda.org/articles/aa/pdf/2015/08/aa26449-15.pdf (Eq. 17 for the rotation - splitting relation)
// https://arxiv.org/pdf/1401.3096.pdf (Fig 13 and 14 for the evolution - rotation relation in SG and RGB) 
// https://arxiv.org/pdf/1505.06087.pdf (Eq. 3 for determining log(g) using Teff and numax, used for getting the evolution stage in conjonction with Fig. 13 from above)
// https://iopscience.iop.org/article/10.1088/2041-8205/781/2/L29/pdf (Inertia and Height relation)
// https://arxiv.org/pdf/1707.05989.pdf (Fig.5, distribution of surface rotation for 361 RGB stars)

// ------------------
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "version_solver.h"
#include "data.h"
#include "string_handler.h"
#include "interpol.h"
#include "noise_models.h" // get the harvey_1985 function
#include "solver_mm.h"


Params_synthetic_star mixed_modes_star(const int el, const VectorXd nu_l0, const VectorXd width_l0, VectorXd height_l0, const double delta0l_percent_star)
// Requirements: nu_l0, width_l0, height_l0, must be set into cfg_modes.params
{
	std::random_device rd;
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	//std::uniform_real_distribution<double> distrib(xmin,xmax);
	std::uniform_real_distribution<double> distrib(0 , 1);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen_m(seed); 
	std::normal_distribution<double> distrib_m(0.,cfg_star.sigma_m);

	int en, el, np_min, np_max;
	long double r, tmp, resol, c, xmin ,xmax, delta0l_star;
	VectorXi posOK;
	VectorXd tmpXd, Dnu_p, DPl, ksi_pg, h1_h0_ratio;
		//nu_p_l1, nu_g_l1, 
	VectorXd nu_m_l1, height_l1, height_l1p, width_l1, a1_l1; // Simulate a single harvey profile

	Data_rot2zone rot2data;
	Params_synthetic_star params_out;
	Data_eigensols freqs;

	// Fix the resolution to 4 years (converted into microHz)
	resol=1e6/(4*365.*86400.);

	// ----- l=0 modes -----
	// This section generate l=0 modes following the asymptotic relation of p modes, and make
	// rescaled width and height profiles for the star using the solar width and height profiles

	// Use fmin and fmax to define the number of pure p modes and pure g modes to be considered
	np_min=int(floor(cfg_star.fmin/cfg_star.Dnu_star - cfg_star.epsilon_star));
	np_max=int(ceil(cfg_star.fmax/cfg_star.Dnu_star - cfg_star.epsilon_star));
	np_min=int(floor(np_min - cfg_star.alpha_p_star*std::pow(np_min - cfg_star.nmax_star, 2) /2.));
	np_max=int(ceil(np_max + cfg_star.alpha_p_star*std::pow(np_max - cfg_star.nmax_star, 2) /2.));  // The minus plus is there because (np_max - nmax_star)^2 is always positive
	
	if (np_min < 1)
	{
		np_min=1;
	}
	std::cout << "np_min =" << np_min << std::endl;
	std::cout << "np_max =" << np_max << std::endl;
	
	// ------- l=1 modes ------
	// Use the solver to get mixed modes
	delta0l_star=-el*(el + 1) * cfg_star.delta0l_percent_star / 100.;
	freqs=solve_mm_asymptotic_O2p(cfg_star.Dnu_star, cfg_star.epsilon_star, el, delta0l_star, cfg_star.alpha_p_star, cfg_star.nmax_star, cfg_star.DPl_star, 
								  cfg_star.alpha_g_star, cfg_star.q_star, cfg_star.sigma_p, cfg_star.fmin, cfg_star.fmax, resol, true, false);

	// Filter solutions that endup at frequencies higher/lower than the nu_l0 because we will need to extrapolate height/widths otherwise...
	posOK=where_in_range(freqs.nu_m, nu_l0.minCoeff(), nu_l0.maxCoeff(), false);
	nu_m_l1.resize(posOK.size());
	for (int i=0; i<posOK.size();i++)
	{
		nu_m_l1[i]=freqs.nu_m[posOK[i]];
		if (cfg_star.sigma_m !=0) // If requested, we add a random gaussian qty to the mixed mode solution
		{
			r = distrib_m(gen_m);
			nu_m_l1[i]=nu_m_l1[i]+r;
		}
	}
	
	// Generating widths profiles for l=1 modes using the ksi function
	Dnu_p=freqs.dnup;
	DPl=freqs.dPg; 

	ksi_pg=ksi_fct2(nu_m_l1, freqs.nu_p, freqs.nu_g, Dnu_p, DPl, cfg_star.q_star, "precise"); //"precise" // assume Dnu_p, DPl and q constant
	h1_h0_ratio=h_l_rgb(ksi_pg); // WARNING: Valid assummption only not too evolved RGB stars (below the bump, see Kevin mail 10 August 2019)
	
	height_l1p.resize(nu_m_l1.size());
	for (int i=0; i<nu_m_l1.size();i++)
	{
		tmp=lin_interpol(nu_l0, height_l0, nu_m_l1[i]);
		height_l1p[i]=tmp;
	}

	height_l1p=height_l1p*cfg_star.Vl[0];
	height_l1=h1_h0_ratio.cwiseProduct(height_l1p);
	width_l1=gamma_l_fct2(ksi_pg, nu_m_l1, nu_l0, width_l0, h1_h0_ratio, el);
	
	// Generating splittings with a two-zone averaged rotation rates
	rot2data=rot_2zones_v3(cfg_star.rot_env_input, cfg_star.rot_core_input, " "); //rot_env, rot_c		
	a1_l1=dnu_rot_2zones(ksi_pg, rot2data.rot_env, rot2data.rot_core);


	params_out.nu_l0=nu_l0;
	params_out.nu_p_l1=freqs.nu_p;
	params_out.nu_g_l1=freqs.nu_g;
	params_out.nu_m_l1=nu_m_l1;
	params_out.width_l0=width_l0;
	params_out.width_l1=width_l1;
	params_out.height_l0=height_l0;
	params_out.height_l1=height_l1;
	params_out.a1_l1=a1_l1;
		
	return params_out;
}