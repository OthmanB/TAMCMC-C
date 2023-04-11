# Version history #

### v1.83.X Tweaks ###
	- 1.83.1: Adding the possibility of settting priors values for ferr in the hyper prior section. The syntax is the same as the one used for the common part, except that the first column must be the fref frequency (as before)
    - 1.83.2: In bin2txt, adding as argument the last sample to be exctracted. Warning: This may break codes that call bin2txt
    - 1.83.3: Adding A constant model of Width for AppWidth_v4
	- 1.83.4: Bug Fix in reading model files when there is only one spline column (extra priors) instead of 1 for the initial guess + 1 for the prior name + N for the parameters of the prior. Another bug fix concerns the fact that the RGB v4 model was not using the spline! It has been commented by mistake at some point
	- 1.83.5: Update in tools/getstats : Three new optional parameters were added: first_index, last_index and Sample_period. This now align the capabilities of getstats with the capabilities of bin2txt. 
	- 1.83.6: Update of tools/getmodel: Added data_type as optional parameter. If set to 'data_type=range', it expects the argv[1] to have 3 values only: xmin,xmax,resolution. And it generates the x-axis for the model within getmodel. This differs from 'data_type=file', which keeps the former behavior (expecting a data file as input for argv[1])
	- 1.83.7: Update on model model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4 relying on the ARMM: with outparams = True, the created params.model file returns a new table (at the end of the file) that contains information specific to the mixed modes such as: l, fl_mixed, fl_p, fl_g, ksi_pg, h1_h0, a boolean saying if the mode was included in the model or just calculated and thrown.
		Will need later on to implement that change in all ARMM-related models
	- 1.83.8:  Binary files are now generated into a bin directory. More importantly, Update of the Alm model to add the latest version of the Alm (V2) that includes the possibility of a triangular filter and some bug fixes in the way the integral is performed. This resulted in these changes:
			  1. The Alm directory is now in external/Alm (instead of external/integral) and contains several testing elements in python.
	          2. Update of the models.cpp path to activity.h according to (1)
			  3. Update of the prior in priors_calc.cpp for the aj_fit and ajAlm model in order to (also) forbid theta > 90 - delta/2 when in 'gate' case
			  4. The gauss filter now works (!)
			  5. Changed io_ajfit.cpp in order to allow to have filter_type used as an argument to change the filter type between 'gate', 'gauss' and 'triangle'
		THINGS WHICH WERE NOT CHANGED: ajAlm model (fitting directly the spectrum with the activity model) still have a hardcoded filter_type = 'gate'
		      This will need to be changed eventually

### v1.83 New tool and improvments ###
	- Adding a quick_samples_stats.cpp and quick_samples_stats.h in tools/ that contains functions to compute the mean, median and stddev from a VectorXd
	- Using quick_samples_stats in bin2txt_params.cpp in order to show the median and the stddev while unpacking the binary files. This is for example useful in small models such as the Gaussian model
  - Models with AppWidth now require that numax (and optionaly its error) is provided by the user. Previously if this was not provided, the code was calculating it. But this is highly instable due the usally limited number of modes in RGB stars
  - Appourchaux Gaussian priors have been revised so that these are using 5% numax if err_numax is not provided in the model file. And use err_numax otherwise.

### v1.82 New Model ###
	- Adding a model 'model_ajfit' that allows you to fit data for a2, a4, a6 using Alm in order to determine latitudes of active regions
	- Warning: In order to get compatibility with ajfit models, a critical change was made in the config.cpp::setup : I Made direct use of the (global) data.data_all instead of the (temporary local) data_in. If a subroutine within read_inputs_files() changes data_all, this may have impact on the routines set prior model_ajfit. A close monitoring of the situation is required in the future, to see any unusual behavior in e.g the ranges of the data
	- Bug Fix: 
		 * Corrected the missing nu_cl in the function 'decompose_Alm_fct()' and in 'build_l_mode_ajAlm()'. This affected models 'model_MS_Global_ajAlm_HarveyLike' (direct or indirect)
		 
### v1.81 Improvment ###
		- Adding in the model with activity  'model_MS_Global_ajAlm_HarveyLike' the possibility of not using Alm directly: Using decompose_Alm option we can now converts it to a2,a4,a6 terms. The user can then choose which terms are used and wich ones are not. It allows to
			  control the fit accuracy, in line with the simulations involving aj coefficients
			  	- decompose_Alm = -1 : Use Alm directly (no decomposition)
			  	- decompose_Alm = 0  : use a2, a4, a6
					- decompose_Alm = 1  : use a2, a4
			  	- decompose_Alm = 2  : use a2 only
			By default, decompose_Alm is set to -1 in the io_ms_global.cpp

### v1.8 Fixes and test  ####
	- Fixes the problem of a3 that is always positive in the posterior: It was due to some residual lines of codes in priors_calc.cpp inteded for the obselete models with activity as a power law.
	- Model with activity 'model_MS_Global_a1etaAlma3_HarveyLike' renaming 'model_MS_Global_ajAlm_HarveyLike'  to verification and full testing [TESTING]

### v1.73 New model and tools ###
  - WARNING: I NOTICED THAT THE HFACTOR MODIFICATION MAY BE A PROBLEM IN v3 models. The io_asymptotic is clerarly quite messy and would require cleaning once we 
             have converged toward a nice solution for the fitting
  - Addition: I have added a corrective factor Wfactor, similarly to Hfactor, it changes Wl=(1 + zeta(nu)).W(l=0) to Wl=(1 + Wfactor.zeta(nu)).W(l=0).
  				    It is an optional parameter set to 1 when not provided so that it is retrocompatible with older version. However:
  				    	- It requires to add Wfactor as an argument in the .model file (GU prior is recommended at that stage with Upper bound at 1 to avoid negative widths)
  				    	- It was tested and explicitly implemented only for model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4
  				WARNING: FULL TESTS MADE ONLY FOR THE cubic spline (hermite spline should also work though)
  - Addition: 
  		- tools/quickshow.py: Allows you to quickly visualise a params.model file (either created by the IDL postprocessing file or by the debug 'outparams=true' option inside the models)
  		- tools/Gaussfit_tools/init_fit.py : Adding the capability of using a data file (instead of a sav file) in order to create an initial configuration file for fitting a gaussian. Because usually those data files are already made as the result of a fit, their range is restricted and thus, the user will have to be carefull with the lower harvey profile: The results of the gaussian fit are likely to give weak constrain on it due to the lack of data to support the fit in that region (if the original data file does include low-frequencies)
  		- tools/recale_height.py: A small function that allows you to divide the power of a spectrum by a certain factor (default is factor=1000). This could be usefull to have 
  			better convergence properties of the MCMC as it avoid to have large covariance terms (errors on the height are usually of the order of 20% of the input). The matrix inversion  process when evaluating the best MCMC 'step size' can indeed suffer from important error propagation issue if the covariance terms are too different from approximately 1 to 10.
 	- Changes: The auto prior for delta0l: Use of a narrower prior: +/-2% of Dnu (instead of 5%)
 	
	- Model for RGB stars : model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v4
		 Changes are in external/ARMM, introduction of polyfit.cpp/h in sources/headers, adding external/spline. 
	     * Contrary to v3, it has the following capabilities:
           1. model_type = 1 : The user can choose to use the l=0, shifted with d01 in order to generate l=1 p modes.
                  In this case, we make use of:
             solve_mm_asymptotic_O2from_l0(nu_l0_in, el, delta0l, DPl, alpha_g, q, resol, bool returns_pg_freqs=true, verbose=false, freq_min=fmin, freq_max=fmax)
           2. OR model_type=0: it can use an O2 polynomial:   
               solve_mm_asymptotic_O2p(Dnu_p, epsilon, el, delta0l, alpha_p, nmax, DPl, alpha_g, q, fmin, fmax, resol, true, false); //returns_pg_freqs=true, verbose=false
              Here, Dnu_p, epsilon are calculated using a linear fit of l=0 frequencies.
              alpha_p, nmax are determined using a 2nd order fit of l=0 frequencies.
              The switch between these two option is made using a control parameter at the end of vector of parameters (replacing sigma_limit of v2 and v3 and named model_type)
         	The fit incorporate a list of Nerr parameters situated between ferr=[Nmax+ lmax + Nfl0 + 8] and [Nmax+ lmax + Nfl0 + 8 +Nerr] that are used to correct from the bias of the exact asymptotic model
          These are associated to a list of Nerr list of fixed frequencies between fref=[Nmax+ lmax + Nfl0 + 8 + Nerr] and [Nmax+ lmax + Nfl0 + 8 + 2 Nerr] that provided frequencies of (a) case of bias_type = 0 (cte):
               The reference frequencies fref are used as 'windows' for determining region of 'constant bias': If a mode lies within a window [fref(j),fref(j+1)] of the list, then we apply the correction ferr(j)
               This means that multiple modes may have the same correction factor ferr(j) and that some window may have 0. Be aware of this when post processing: Each window must be cleaned of solutions with 0 modes
           (b) case of bias_type = 1 (Cubic spline, ie strongly correlated as it is twice continuously differentiable):
               The reference frequencies fref are used as 'anchors' for a spline defined by the values ferr.
               The spline fitting is performed here in this function using the spline module from https://github.com/ttk592/spline/ also forked to my repo
           (c) case of bias_type = 2 (Hermite spline, ie less correlated variable as it is once continuously differentiable):
               The reference frequencies fref are used as 'anchors' for a spline defined by the values ferr.
               The spline fitting is performed here in this function using the spline module from https://github.com/ttk592/spline/ also forked to my repo

### v1.72 Bug Fixes and New model ###
	- Bug fixes 
	- Model for RGB stars with constant width (experimental): model_RGB_asympt_a1etaa3_CteWidth_HarveyLike_v3

### v1.71 EXPERIMENTAL Improvment ###
	- Adding a check of the finitness of the new proposal vector. If not successful, it retries to generate a new random vector
	  This will come at an extra computation cost of Nchain*Nparams ~ 500 operations

### v1.70-dev New model + Improvements + Bug Fix  ###
	* Bug Fix:
  	- Some function were putting a a3 coefficient for l=1. a3 does not exist for l=1 and it is corrected [DONE]
  * Improvments:
  	- Use of acoefs.cpp in order to hande the aj decomposition in a cleaner way   [DONE]
  	- Use of a Qlm(l,m) function to avoid mutliple repetition in the code. Qlm(l,m) includes Dnl = 2./3  [DONE]
		- Reorganise existing functions to use the new acoefs functions when necessary  [DONE]
	  - Rearanging and cleaning the functions in build_lorentzian.cpp in order to avoid repeating sections that handle the optimisation [DONE]
  	- Update for Alm: Instead of using directly Alm(theta,delta), we give the user possibility to decompose Alm~F(a2,a4,a6) and to decide which aj to not account for
  	- The logLikelihood is computed only if the logPrior is not -INFINITY. This avoid us to compute the model when it is not necessary ==> performance improvement  [DONE]
  	- Removed the model *acta3 that were obselete [DONE]
  	- Systematic use of eta0 = 3.pi/(G.rho_sun) . (Dnu_sun/Dnu)^2 in models.cpp and io_ms_global.cpp for the centrifugal force instead of eta = eta0.a1^2. Additionnaly eta0 was corrected from a bug [DONE]
  	- Optimisation of the structure of priors_calc.cpp inside MS_global case.
	* New models:
		- Adding acoefs capabilities up to a6: Compose frequencies with a-coefficients and decompose frequencies into a-coefficients [DONE]
		- Model that handle <a1>_l, <a2>_l, <a3>_l,<a4>_l,<a5>_l, <a6>_l in linear function and for MS stars
		- Model that handle a1(nu,l), a2(nu,l), a3(nu,l), a4(nu,l), a5(nu,l), a6(nu,l) in linear function and for MS stars

### v1.65-dev Improvments ###
	* Refactoring the section external code that can compute Glm/Alm function in preparation of the more general solver for any kind of 
          Activity zone. Old function named ylm.cpp and ylm.h were replaced by activity.cpp and activity.h. Note the although functions were
          renamed [xxxx]Alm[xxxx] I currently only set Alm in "gate" mode and did not update the model itself and checked it consistency.
          Ideally, the "gate" or "gauss" mode should be a global parameter
        * Tuning of prior_calc.cpp::priors_Harvey_Gaussian() in order to enable a fit of the Gaussian envelope in a similar fashion that the MCMC part  
          Of my IDL code 'Preset-Analysis-v8.1', ultimately replacing it (for the Gaussian fit at least for the moment)
        * Full support of the Harvey_Gaussian fit.
        * New tool to automatically generate a .data and .model file for the Gaussian envelope fit: init_fit.py
          Most elements and logic is taken from init_fit.pro of Preset-Analysis-v8.1'. Note that it will still need guesses for Amax and numax. Those are in principle provided by my code 'envelope_measure.pro' 
 				* Creation of several function in order to bridge mode_envelope.pro (that gives guesses on Amax and numax) with the init_fit.py function
 				* Adding a tool within init_fit.py that can generate a .model and .data file from the results of at Gaussian fit in order to make a second Gaussian fit (pass 2) that will enhance the noise background guess. This reveals to be sometime usefull due to the fact that priors set automatically in init_fit.py (pass 1) are sometimes poor
	* Fix issue where simulation files that specify trunc_c [Value] cannot be read due to the fact that the program expect trunc_c Fix [Value]

### v1.64-dev Improvment ###
	* New model:
		- Adding a model named 'model_RGB_asympt_a1etaa3_AppWidth_HarveyLike' that: (1) Handle a pure asymtptotic relation and (2) use random quantities as hyperparameters to evaluate the inaccuracy of that model and to still be able to fit the spectrum exactly despite those inaccuracy
	* Safety:
		- Ensure that the user is warned when the command line is typed with indexing starting at 0 (e.g. ./cpptamcmc execute 1 0 1) instead of 1 (./cpptamcmc execute 1 1 2)
### v1.63-dev Improvment ###
	* Improvment:
		- Implementation of clm for l=3 modes. allows to account for a3(n,l=3). Previously approximated by cln=0
### v1.62-dev Bug Fix and Improvments ###
	* Bug Fix:
		- properly linking optimised_xxx_a1a2a3* to a1a2* a1la2* and a1nla2* models (was mistakenly linked to a1etaa3)
	* Improvments: model_MS_Global_a1etaGlma3_HarveyLike() that use Gizon 2002 AN prescription for describing Activity effect on modes
		- Adding build_l_mode_a1etaGlma3(), optimum_lorentzian_calc_a1etaGlma3() and model_MS_Global_a1etaGlma3_HarveyLike()  [100%]
		- Adapting model_def and model list  [100%]
		- Adapting prior_calc.cpp for the new model [100%] 
		- testing [100%]: Note that I let epsilon to take negative values. I suspect this is important to allow a2_AR<0
		
### v1.61-dev Improvments ###
	* Improvments:
		- Adding an optional boolean switch in the list of models to allow outputs of the mode parameters in the form of an ascii table (used in getmodel)
		- Adding optional boolean switch in model_def.cpp for call_model() and call_model_explicit : This to ensure usability by getmodel 
		- Adding a writing function for ascii table outputs
		- Update of getmodel: It will output an ascii table of all the best fit parameters

### v1.6-dev Improvments ###
  * New models:
  	* Adding core functions in build_lorentzian.cpp to handle a2 coefficient (see Gizon,2002, https://www2.mps.mpg.de/homes/gizon/PAPERS/an_323_251.pdf) 
    Need to add a2 as a polynomial of O1 to handle a2=cte or a2=linear(nu) or a2=quadratic(nu) for a1l fit and a1 fits (replaces eta parameter)
  	* Add model in models.cpp that handle a2 for MS stars in the case of a1 is constant 
  	* Update io_ms_global.cpp to initialise properly model parameters for the new models
  	* Update priors_calc.cpp to handle a2 priors
  	* update models_list.ctrl 
  	* update model_def.cpp in the switch cases
  	* update the default_errors.cfg

  * Bug fixes:
  	- Update of the CMakeList.txt with some (missing) -fopenmp option (guarantees omp is used whenever possible)
    - Removed the 4 threads hard-coded for the ARMM solver and explicitly impose num_threads = Nchains in MALA.cpp to paralelise the Tempered chains. It worked better like this
    - a3/a1 and a2/a1 limits where systematically wrongly set: priors_calc.cpp was set such that a1 was extracted from the slot of only classic models. A new extra_priors[5] switch allows to mitigate the issue
  * Prospect for subversions:
  	- Extend testing and implementation of the a2 to a1l, a1nl and a1n MS models. Currently models.cpp do have those, but priors_calc, model_def and io_ms_global.cpp are incomplete and force the stop of the program if one tries to use them

### v1.57-dev New model ###
  * Adding a model with free widths for asymptotic mixed modes fitting (model model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3)
  * Adding the parallelised version of ARMM (with 4 threads)

### v1.56-dev Minor improvments ###
  * Change of the sigma for the gaussian smoothness condition for fl1p modes for asymptotic models
  * Adjustments of default_errors.cfg for asymtptotic models
  * Enlargement of the range for the fl1p modes in asymptotic models
  * Impose abs() when generating initial guesses for the Width using Appourchaux relation in asymptotic models
  * Removal of the smoothness condition for the l=2 in asymptotic models
  * Removal of the d02 prior in asymptotic models

### v1.55-dev New models and improvments ###
	* Bug Fix:Fix bug when initializing delta0l with fix_auto [DONE]
	* Adding : sigma_Hl1 (replaces sigma_p) [DONE]
	* Adding a model that has the fl1p modes as free parameters... priors and initial guesses are automatically set [DONE]
	* Improvment: taging modes for asymtptotic model as RGB using error_default.cfg standard: Frequency_RGB_l [DONE]
	* Upgrading Dnu determination in priors_calc.cpp using l=0 mode fitting


### v1.52-dev Improvoment ###
	* Bug fix: There was multiple NaN that were still not corrected in v1.51. It was comming from the fact the Hl0 interpolation at fl1 frequencies 
			   Could lead to negative Hl1 heights sometimes. The fix was to enforce the Hl1 to be positive when generating the model, effectively avoiding NaN in the likelihood
	* Optimisation: Use of the goto statement when f=-INFINITY priors_calc.cpp so that we skip any calculation that would follow if f=-INFINITY
	* Improvment: 
		- Enforcing positivity of Width parameters when those have Gaussian priors
		- Adding a small suite of tool in IDL (data_extract_IDL) into the tools directory. This allows to (a) extract samples (all or a subset defined by the user) and save them into an IDL sav file. For this, bin2txt must be compiled and present into the tools directory. (b) Generate models using the getmodel binary program to plot the best fit on the top of data and save this in an eps file.
		- getmodel::check_retrocompatibility() function is now having status=2 by default, which means that by default we assume that models follow a compatible format with any code abobe 1.3. This means that getmodel will no longer need to be updated to handle new models. 
		
### v1.51-dev ###
	* Adding an upper limit to the random quantities generated by sigma_m parameter
	* Several bug fixes
	* Remaining problem: 
		(1) All of the frequencies seems to be scarcely sampled... default errors too big? Will need to try to use a separate Frequency_l (Frequency_RGB_l) and set it properly whenever needed (a threshold on Dnu)
		(2) We some extremely strange NaN on the likelihood when testing the full modeling of KIC 008751420. This problem does not appear when using a narrow range of p modes. What is even more strange is that this happens while the proposed position seems to be the same as the initial position (according to the debug message shown when NaN arises). It could mean that there is a serious issue in the random generation or in the core of the code that really need to be investigated.
		
### v1.5-dev [100% Implemented... Need testing] ###
	* Improvments: adding the possibility to fit the spectrum with the asymptotic relation of the mixed modes (l=1 and l=2)
	  in the case of RGB models (handled by a new file, io_rgb_global.cpp)
	  	- update models.cpp to include the model [100%]
	  	- convert the mixed mode solver functions from python to c++ by starting by the core solver function (that can be tested within python by a wraper) [100%]
	  	- update of config.cpp [100%]
	  	- update of model_def.cpp [100%]
	  	- create and write io_rgb_global.cpp [100%]
	* Bug fix: removing useless '...' into the linking process of the CMakelist.txt, that were leading to show warnings during the compilation [100%]

### v1.4.31-dev Bug fix: 30/09/2020 ###
    * getmodel was not properly updated to handle 1.3.3 update on handling the models names using a file (.list file). 
      The consequence is that no model was recognized and the getmodel tool could not operate. 
      The fix consist on making an overload fonction Config::convert_model_fct_name_to_switch(const std::string model_name, const Data_Basic models_ctrl)
      that take for argument the content of a list file. The content of a list file is read by the usual Config::read_listfiles() and the models_ctrl.list file that contain the lists must be in the execution directory of the getmodel program or of the program that run itself getmodel
      Note that this file can be found in the Config/default/ directory

### v1.4.3-dev improvements: 20/08/2020 ###
		* For local fitting: Adding the option for an automatic set of the Height prior. In auto mode, A Jeffrey's prior is used, with an upper bound given by the input value multiplied by a factor X given by the user (typically 3 is fine). The lower bound is fixed as Y times lower than the actual input  
		The reason is that for RGBs, a single common value of the upper bound for the Jeffrey prior lead to a too large dynamic range for low SNR modes (e.g. low freq or l=3).
		In the .model file, the new prior can be set by setting e.g. (here X=3 and Y=10):
			               Height            Fix_Auto          10.000000      3.000000
		Instead of typically:
		                   Height            Jeffreys          1.000000          1                 10000.00
		The change is also valid for amplitude fitting (instead of height). In that case, the Fix_Auto option will multiply the max height (input*user_value) by \Delta\nu/3 to get the upper bound. The line must be like this:
						   Amplitude            Fix_Auto          10.000000    3.0000
### v1.4.23-dev Bug Fix: 02/04/2020 ###
        * The names for the new inclination parameters in io_ms_global.cpp contained a space. This posed a problem when parsing the string with bin2txt. 
          The change consist in removing the space. e.g.: 'Inc: H1,0' becomes 'Inc:H1,0'
  
### v1.4.22-dev Bug Fix: 17/03/2020 ###
        * In the rare eventuallity that the user was giving an inclination input leading to the Sum of H(n,l,m) over m to be 1,
          A round-off error was making the log_uniform(0, 1, H(n,l,m)) prior (see line 64 of priors_calc.cpp) to be infinity.
          The fix was to add 1e-10 to the upper boundary. 
### v1.4.21-dev Bug Fix: 16/03/2020 ###
	* Fix name of variables rthat were leading to crash due to forbidden caracters ("=") in the variable names

### v1.4.2-dev improvements: 17/02/2020 ###
	* Adding the possibility to fit Heights of (l,m) instead of the stellar inclination
	  For that purpose, the following is going to be implemented:
		- New models for the io_MS_Global class of models:  model_MS_Global_..._Classic_v2 and model_MS_Global_..._Classic_v3 [Status: 100% developped, 99% tested]
		- A New model for the io_local class of models: model_MS_local_basic_Hnlm [Status: 100% developped, 99% tested]
	
### v1.4.1-dev improvements and fixes ###
    * Bug fix on local fitting models (io_local.cpp):
    	- When providing only a single degree for the local fitting case, the code crashes as it assumes that lmin=0. 
          This limitation has been fixed
        - Priors are mostly improperly set if Dnu is not provided as an argument on the top of the .model file (with the exclamation mark and followed by the value)
          Corrected by:
          		- For Widths: Warning the user that the default is in that case a maximum width of 20 microHz (instead of the dynamically defined max width, based on Dnu/3)
          		- For Frequencies: The GUG priors was using 0.01*Dnu for the Gaussian edges. Now it uses 0.01*(window_max(n,l) - window_min(n,l)). This is valid only for local models. models relying on io_ms_global.cpp still use 0.01*Dnu.
        - Remaining bug on the noise background: The white noise detection was done improperly if the user provded -2 values in the last line of the section '# Noise  parameters: A0/B0/p0, A1/B1/p1, A2/B2/p2, N0'. These were not removed in the noise parameter vector.
          This is fixed.
    * General bug fix: When reading .model files, the relax parameter for Heights an Widths were reversed. This affect all previous version since the origin of the MCMC code.
    		This is fixed.

### v1.4.0-dev Local fit implementation ###
	* Adding a whole new serie of function to handle a local fit. Warning: Changes require slight changes in the .model file. The various fitting ranges have to be defined 
          using the '*' marker [DONE] [TESTED]
 
### v1.3.3-dev New model and improvements ###
	  * Add two models with fit of the Width from Appourchaux 2012 instead of individual widths. 
	    This is made in a transparent way, using the input widths for defining the guesses. Thus, no need to change drastically the structure of .model files
	  	Priors on the parameter of the model are hard-coded at the moment. See io_ms_global.cpp for their value (typically 10-20% of the expected values) [DONE] [NEED TEST ON DALMA]
	  * When starting a new process, a 'version.txt' file is writtent in the object directory. This file gives the code version used for the processing [DONE] [TESTED]
	  * Adding a hard-coded limit on a3/a1 < 0.2. This to avoid very unphysical solutions (that lead to a star rotating in the opposite direction in the pole [DONE] [TESTED]
	  * Redesigning the way likelihoods, models, priors and prime priors are handled: Instead of having them hard coded, they are now defined along with their switch case into *.list 
	    files within Config/default. This significantly simplifies the implementation of new models.
	    
### v1.3.2-dev Improvements ###
      * Adding the possibility to fit amplitudes instead of Height by specifying  fit_squareAmplitude_instead_Height   [bool value]   into the .model file [DONE] [TESTED]
      * Add the possibility to change the priors for Height/Amplitudes/Width in the .model file [DONE] [NEED THOROUGH COMPARATIVE TESTING WITH EARLIER STABLE VERSION] 
        Note that we are permissive regarding the denomination Height or Amplitude. However, an explicit fit of squaredAmplitude could be requested. See example of .model file
      * Add the possibility to change frequency prior between GUG or Uniform. [DONE] [NEED THOROUGH COMPARATIVE TESTING WITH EARLIER STABLE VERSION] 
        For GUG the syntax is
      			Frequency     GUG     -1    -1    -1     [sigma1]    [sigma2]       (-1 are replaced by values of the table of frequencies)
      	For Uniform the syntax is
      			Frequency     Uniform  -1  -1   -1      (-1 are replaced by values of the table of frequencies)
      * Reorganising io_ms_global using new functions: initialise_params(), fill_param() and add_params() 
        This reduced the size of this program significantly by removing redoundant commands
        This allows much easier generation of a model. Note that the function could be used to create any kind of vector of parameters
      * Minor esthetic improvements in the text outputs
      * Few minor typos and fixes into the .md files [DONE]
      * Starting a basic documentation in tex [DEV] [POSTPONED: Replace by Wifi on Github]
        
### v1.3.1 Minor improvements/Bug fix ###
	* Correcting a typo in getmodel.cpp that prevented the compilation of getmodel tool
	* Further cleaning of useless lines, in particular in outputs.cpp
	* Removal of a debug text shown when reading the header of the parameters
	* Adding model_MS_Global_a1etaa3_Classic which allows to directly fit a1 and inclination (rather than sqrt(a1).cosi and sqrt(a1).sini)
	* Adding adaptive scheme for MS_Global_a1eta3: If the user provide splitting_a1 and inclination instead of the expected parameters (sqrt(a1)cosi, sqrt(a1)sini), then 
	  the program convert and adjust the inputs so that they are compatible with the model
	* Removal of some lines that were added in v1.3 and that stop the program if using 4-lines format for the s2 inputs (instead of 3 in simulations)
	* Adding a failsafe if the user request to read a y-data column that does not exists
      * FUTURE:
  		* Write comprehensive md files /doc file with compilation configuration and execution instructions 
  		* Verifying the Gaussian fit case
  		* Handling of complex space fitting
  		* Continue clean handling of errors perhaps using 'Either'? see https://hackernoon.com/error-handling-in-c-or-why-you-should-use-eithers-in-favor-of-exceptions-and-error-codes-f0640912eb45
  				* Create a dedicated message.cpp / message.h that handles messages and errors (move Diagnostics::file_error in there)		 	
  				* Better structure for output.cpp ... use header to save constants


###  v.1.3 Major Improvements ###
	* Detection of openmp and option to deactivate it in order to compile with clang that is not compatible with openmp
	* Detection of GSL and option to deactivate it if this library is not available on the considered machine
	* Detection of GNUPLOT. It is mandatory
	* Option to compile on DALMA (unified cmakefile)
	* Option to show the version of the program/compiler added.
	* Variables start_index_ids and last_index_ids are now obselete. To select process to be run, these are now passed as argument of the executable
	* Better handling of io_MS_Global 
		* In the config_default.cfg, only prior_fct_name is now required to specify which model should be used. 
		When prior_fct_name=io_ms_Global (capital letters matter!) the code calls io_ms_global.cpp, that read the .model file. The .model file
		must now contain a line that specify the model name (common variable model_fullname).  
		Basically, prior_fct_name can act as a switch between models in the future
		* Use of an error function for repetitive error messages 
		* The variable average_a1nl is now useless in .model files. This due to the proper handling of the model names in io_ms_global.cpp
		* The c constant that controls the Lorentzian truncation is now replaced by a variable called trunc_c that could be added in the .model file
		  If set to -1, trunc_c=10000., a very large value that is equivalent to no truncation at all.
	* Cleaner Diagnostic.cpp / config.cpp / config_preset.cpp / MALA.cpp
		* Use of an error function for file errors at opening time. 
		* Automatic handling of missing gsl (a single diagnostic.cpp) 
		* Function strsplit, strtrim, dbl_to_str, VectXd_to_Vec, str_to_dblarr, str_to_Xdarr, int_to_str are inside string_handler.cpp and not duplicated 
	* bin2txt and getmodel:
		. Complete integration within the main project, including compilation.
		. Version option support for bin2txt and getmodel 
	* getstats integration within the main project 
		. Adding function that check consistency with MS_Global models that predate 1.3.0
		. Adding function that adapt input vector if plength=10 and the we deal with an MS_Global model

### v1.2.3 Improvements ###
	* Added functionalities:
		The keyword 'average_a1nl' was added. It allows you to decide wether to fit a1(n), a1(n,l) or a1(l).
		NOTE THAT YOU DO NOT HAVE A CONTROL ON THE FIXED/FREE VARIABLES THERE. e.g. if you want to fix a1(l=1) = a11 and let free a1(l=2), you must do it by code edition
		
		The syntax should be as follow: average_a1nl     bool    [0/1]   [0/1]
        
        Ex 1: average_a1nl     bool    1    1" << std::endl;
              This uses model_MS_Global_a1etaa3_* which assumes a1(n,l) = a1
        
        Ex 2: average_a1nl     bool    1    0" << std::endl;
              This uses model_MS_Global_a1n_etaa3_* which assumes a1(n,l) = a1(n)
        
        Ex 3: average_a1nl     bool    0    1" << std::endl;
              This uses model_MS_Global_a1l_etaa3_* which assumes a1(n,l) = a1(l)
        
        Ex 4: average_a1nl     bool    0    0
              This uses model_MS_Global_a1nl_etaa3_* which assumes a1(n,l) are all free

### v1.2.2 Minor Changes ### 
	* Corrected issues:
		* Problem: Initial values for labels and units (empty vector of strings) was not allowing more than two columns input data. This is a problem when the user wish to
		            use a 'best_models_fit.ascii' file generated by the IDL_PostMCMC program
        * Solution: The program was modified such that I initialise the labels and units vector with 5 columns instead of 2.

### v1.2.1: Major Changes ### 
    * New functionalities:
   		* The program now uses openmp to parallelise the parallel chains computation. The number of used thread is therefore 
   		   controlled by the OMP_NUM_THREADS=X with X the number of used cpus. Gain are typically optimal when OMP_NUM_THREADS=Nchains/2
   		   Gains are null if OMP_NUM_THREADS>Nchains.
   		* Added scripts for running the program on NYUAD-DALMA (slurm scripting)
   		* Multiples instances can now run using the same executable (required by parallelisation)
   		* The program is now executed as a command line in order to process a specific range of object of the 'ids' list in config_presets.cfg
   		   The command line is of the following syntax: ./tamcmc execute [0/1] [first line number to process] [last line number to process]
   		   Warning : The line id index begins at 0 while the the argument of the command line are begining to 1.
   		   Warning2: The variables 'start_index_ids' and 'last_index_ids' in config_presets.cfg are now obselete
   		             These will be removed on a future version 
   	* Corrected issues:
   		* The clock was not working properly due to parallelisation. This has been rewritten by Benoit Marchand
   		* There was a bug on the models.cpp functions that was fixing the l=0 asymetry to eta (thus to almost 0)
   		* The variable do_backup_cfg_files was creating a backup of the data file while this is controlled by do_backup_input_file
   		
   	* Identified issues:
   		 * Problem: In some rare occasions, the following error appears:
			terminate called after throwing an instance of 'std::ios_base::failure[abi:cxx11]'
			what(): cannot open pipe gunplot -persist: isotream error
		    Reason: The problem is due to memory limits. If the system runs out of memory, gnuplot fails
		    Solution: Reduce Nbuffer so that memory usage remains within your system specifications. A future
		    	      will be more optimised in terms of memory usage so that this kind of situation may be more
		    	      unlikely to face
		    	      
		* Problem: In some rare occasions, 'Nan' are return in the vector of parameters. These models are rejected automatically and do not prevent the proper execution of the code

		* Problem: In some rare occasions, build_lorentzian encouter 'Nan' for some parameters, which make the code to unexpectely stop.
		   Solution: TEMPORARY SOLUTION IS to restart the analysis from where it stopped (restore = 3)

		* Problem: When in restore=3 and if Nbuffer is changed relative to the old run, saving of all variables might not always
		            happens when it should.
		   Reason: It is likely due to the fact the logical test to decide when to save is not accounting well from i0
		   Solution: account for i0 when deciding when to save

		* Problem: Verbose missing when writing the evidence
		   Reason: No verbose yet
		   Solution: Implement a verbose
		   
	    * Problem: No evidence written if all diagnostics are OFF
		   Reason: Investigate the issue

### v1.1.3.4: Minor improvements ### 	
	* Added some verbose and security check in build_lorentzian.cpp. The goal is to see if the bug we encounter once every
	  ~10 Million samples (problem 1 in known existing issues) is coming from there.
	  
	* Corrected issues:
		1. Problem: when restore=3, the acceptance rate plots is not correct.
		   Reason: The exchange file that contains the rejection rate is not defining correctly the x-axis
		   Solution: modify outputs::reject_rate() so that it defines an x-axis which account for the samples that were already processed during the last execution (variable outputs::buf_restore.Nsamples_sofar)

	* Known existing issues:
		* Problem: In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   Solution: Contrary to what I was thinking initially, this error seems not related to a plot function. But rather to the truncation approximation of the Lorentzian.
		   REASON UNKNOWN. CIRCUNSTANCES: When truncating the Lorentzians?

		* Problem: In some rare occasions, the following error appears:
			terminate called after throwing an instance of 'std::ios_base::failure[abi:cxx11]'
			what(): cannot open pipe gunplot -persist: isotream error
		   Solution: I posted the problem in stack overflow
		
		* Problem:  The variable do_backup_cfg_files seems also to create a backup of the input file (normally controlled by do_backup_input_fil)
		   Reason: Logic condition is incorrect
                   Solution: Put the correct logic condition

		* Problem: Verbose missing when writing the evidence
		   Reason: No verbose yet
		   Solution: Implement a verbose

		* Problem: Missing parallelisation
		   Solution: Wait that Benoit Lemarchand review the code and implement the parallelisation in it. Deadline given by Benoit is 15/01/2017

###  v1.1.3.3: Minor improvements ### 
	* Corrected issues:
		* Problem: Error messages appear in machines which have gunplay installed, but not gnuplot-qt. 
		   Reason: in diagnostics.cpp, we switch from a ps terminal into a qt terminal to reinitialise the terminal
 		   Solution: remove the line set term 'qt' as it is not necessary to reset the terminal type.
		   tested: Requires to be tested on othxeon12
		* Problem: When handling a list of process it is hard to keep track of which have failed and which have finished and which are ongoing
		   Solution: Generate a simple output log file on the following text format:
			[Number of the process]   [Name of the process]     [processing_status]
		        Use the following codes for processing_status: No processing, Pending, Ongoing , Finished   	 

	* Known existing issues:
		* Problem: In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   Solution: Contrary to what I was thinking initially, this error seems not related to a plot function. But rather to the truncation approximation of the Lorentzian.
		   REASON UNKNOWN. CIRCUNSTANCES: When truncating the Lorentzians?

		* Problem: In some rare occasions, the following error appears:
			terminate called after throwing an instance of 'std::ios_base::failure[abi:cxx11]'
			what(): cannot open pipe gunplot -persist: isotream error
		   Solution: I posted the problem in stack overflow
		
		* Problem: when restore=3, the acceptance rate plots is not correct.
		   Reason: The exchange file that contains the rejection rate is not defining correctly the x-axis
		   Solution: modify outputs::reject_rate() so that it defines an x-axis which account for the samples that were already processed during the last execution (variable outputs::buf_restore.Nsamples_sofar)

		* Problem:  The variable do_backup_cfg_files seems also to create a backup of the input file (normally controlled by do_backup_input_fil)
		   Reason: Logic condition is incorrect
                   Solution: Put the correct logic condition

		* Problem: Verbose missing when writing the evidence
		   Reason: No verbose yet
		   Solution: Implement a verbose

		5. Problem: Missing parallelisation
		   Solution: Wait that Benoit Lemarchand review the code and implement the parallelisation in it. Deadline given by Benoit is 15/01/2017


### v1.1.3.2: Minor improvements ###
	* Corrected issues:
		- Fixed an issue with the plots of the models (diagnostics.cpp): The smooth coeficient was improperly calculated such that plots of models were not readable in some cases.
		- Fixed an issue with the plots of the likelihood (diagnostics.cpp): When the likelihood was large, a truncation in the exchange file was making the plots awkward.

	* Added functionalities:
		* Handling the outputs is now made easier:
			- In config_presets.cfg, the new variable 'cfg_out_dir=' allows you to set the main directory for all outputs files (diags, diags/pdfs, outputs, restore). 
			  No need anymore to use config_default.cfg for this: Outputs.output_dir, Outputs.restore_dir, Diagnostics.output_dir are bypassed and superseeded by cfg_out_dir
			- Within cfg_out_dir, output directories are automatically created. The following tree is checked and created/completed is required:
				.../[table_ids(i,0)]			
						.../outputs
						.../diags
						     .../pdfs
						.../restore
			   With table_ids specifies the name of the object that is analysed, according to the table at the bottom of config_preset.cfg
		* [BETA TEST] Saving averaged covariance matrix over Nbuffer or rest (depends how many samples remain)... these might provide more reliable covariance matrix out of the Learning process.
		   REQUIRES THOROUGH TESTING BEFORE VALIDATION

		* Added the possibility to backup the input configurations (*.model, *.cfg) and/or the input data (*.data) into the output directory. The backup is made by zipping the files in 'inputs_backup/'
		   This is controled by two new variables in config_default.cfg, do_backup_cfg_files and do_backup_input_file

		* Added the evidence calculation. Outputs are handled by diagnostics.cpp and are save in the diags directory (text file)

		* Added the possibility of showing/not showing the original data (ie, raw data without smooth) in the plots of models. The control variable is 'show_original_data' in config_default.cfg
		   Setting show_original_data is recommended for very noisy data as it improves the rendering of the plots.

	* Known existing issues:
		* In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   This error seems to occurs only in 'nohup' mode and might be related to plot functions. 
		   REASON UNKNOWN. CIRCUNSTANCES: Probably only in 'nohup' mode

		* Missing parallelisation

### v1.1.2: Minor improvements ### 
	* Corrected issues:
		* Problem: In the Misc/IDL_PostMCMC/, the index of the pdf for Height/Width/Amplitude begins at 2.
		   Cause: wrong indexation in the MS_GLOBAL_HEIGHT_L, MS_GLOBAL_AMPLITUDES and MS_GLOBAL_WIDTH
		   Solution: (1) put index=index+1 at the end of the loop and (2) replace index+1 by index when formating the number prior generating the filename
		* Problem: Diagnostic plots are not properly shown (edges missing)
		   Cause: Wrong setup for the margnins in diagnostics.cpp
		   Solution: Change the bottom and left margin from 0.02 (or 0.03) to 0.06

	* Added functionalities:
		* Possibility to switch off/on the smoothness condition by adding 'freq_smoothness  bool  [1/0]   [smoothness coefficient]' in the .model file.
		   Example: 
			- Turn ON the smoothness, with a smoothness coeficient of 1.5 microHz/n^2:  freq_smoothness  bool  1   1.5
			- Turn OFF the smoothness:  freq_smoothness  bool  0   1.5
		* The prior Jeffrey_abs() was added to the list of possible prior. Example of use: For the mode asymetry.

	* Known existing issues:
		* Missing functionality calculating the evidence
		* In some rare occasions, the following error appears:
			tamcmc.out: eigen/Eigen/src/Core/MapBase.h:148: Eigen::MapBase<Derived, 0>::MapBase(Eigen::MapBase<Derived, 0>::PointerType, Eigen::MapBase<Derived, 0>::Index, Eigen::MapBase<Derived, 0>::Index) [with Derived = Eigen::Block<const Eigen::Matrix<double, -1, 1>, -1, 1, false>; Eigen::MapBase<Derived, 0>::PointerType = const double*; Eigen::MapBase<Derived, 0>::Index = long int]: Assertion `(dataPtr == 0) || ( nbRows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == nbRows) && nbCols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == nbCols))' failed.
		   This error seems to occurs only in 'nohup' mode and might be related to plot functions. REASON UNKNOWN. CIRCUNSTANCES: Probably only in 'nohup' mode
		* Saving averaged covariance matrix over Nbuffer... these might provide more reliable covariance matrix out of the Learning process
		* Missing parallelisation


### v1.1.1: Fully functionnal release for Burn-in/Learning/Acquire mode ### 
	* Corrected issues: 
		  * Problem: In the Acquire phase, the sampling proceed for few ~100 samples but then the chain get stuck at one position in the parameter space
			 Cause: Seems to come from the Chain Mixing... If We impose that there is no mixing, then the problem disapear. Further examination of the my program
				shows that only the value of the parameters of the coldest chain were returned.
			 Solution: Impose that the parameters (via params) are initialized to the restore value for all parallel chain (changes in model_def.cpp)
		  * Problem: The initial acceptance rate is rather low for high temperature chains
		  	 Cause: The initial covariance matrix was not properly scaled for chains of T>1. 
		  	 Solution: Corrected by changing in MALA::init_proposal from,
		  	 	initial sigma[m]=std::pow(2.38,2)*(0.1*m+1.)/Nvars into, 
		  	 	sigma[m]=std::pow(2.38,2)*std::pow(Tcoefs[m],0.25)/Nvars
		  	 	This scaling ensures that sigma is always properly initialised accordingly to errors given in errors_default.cfg
		  * Problem: When loading from a restore file, huge quantities of "Warning: Hit epsilon1 in p1"
		  	 Cause: Corrolary of problem 1 and 2. 
		  	 Solution: Solving problem 1 and 2 solves this issue.
		  * Problem: The Learning is made by swapping chains. This should be avoided on the preset mode
		  	 Cause: The dN_mixing variable is updated to dN_mixing = Nsamples + 1
		  	 Solution: Hardcoding of dN_mixing = Nsamples + 1 in Learning mode
		  * Problem: In binary mode, the header does not specify the total number of samples saved so far.
			  Cause: Need to add a line in the header with that information
			  Solution: Added a line in the header file, with keyword 'Nsamples_done'

### v1.1: Second release with added functionalities ### 

### v1.0: Initial release with basic working configuration ### 

	* Known existing issues:
		* Problem with the diagnostics plots: These are not properly shown on some postscript pages. A redefinition of the drawing area is required.
		* Missing functionality calculating the evidence
		* No way to switch ON/OFF the smoothness condition on frequencies
		* Missing parallelisation

	* Added functionalities:
		* In main program: Added an indicator of the swaping rate for the parallel chains
		* In Misc: small routines that can read the binary files have been added (BETA VERSION)
