# README #


### TAMCMC C++ repository ###

Programs to perform asteroseismic analysis using a MCMC optinmisation algorithm. The main codes are written in C++. Some additional code is required to perform the postprocessing in e.g. python (not provided here). 
The algorithm implementation is based on the IDL implementation of my code described in Benomar et al. (2009), see http://adsabs.harvard.edu/abs/2009A%26A...506...15B
This implementation is itself based on the Metropolis-Hasting-Langevin scheme described in Atchade Y. (2006), see https://link.springer.com/article/10.1007/s11009-006-8550-0

* Version 1.4.1, read the file CHANGELOG.md for a full list of the code history.
	The current implementation focuses on the asteroseismic analysis (individual pulsations) of main sequence stars. For other stars, one would need to write
	dedicated functions, but the code is modular so that only minimal effort is required to do so. The version 1.4.0 distinguish itself from previous version by the fact you have the possibility to carry out local fit through the model_name = 'model_MS_local_basic' 
	
* Future implementations should incorporate C++ version of my codes (currently in Python or IDL) to prepare the configuration files (*.model and *.data)

### How do I get set up? ###

* All setup is made through configuration files of extension .cfg. So far, three files are required to configure the program: config_presets.cfg, config_default.cfg, errors_default.cfg. Those can be found in the subdirectory Config/.
In addition, an analysis setup file is required (.model) and a data file (.data). All those files are ASCII files.

* Dependencies: Eigen, Boost, cmake,  GSL (optional), OpenMP (optional)

* No exhaustive documentation is yet available, but these examples might be self explanatory. A complete documentation is scheduled to be released soon,

* Deployment instructions: 
Use cmake to install the program. For the main directory of the program:
```
	mkdir build
	cd build
	cmake ..
	make
```
Some options exists to deactivate functions (eg. GSL and/or OpenMP):
   ```-DWITH_OPENMP=OFF``` : turn off implementation of parrallel computation with OpenMP
   ```-DWITH_GSL=OFF``` : turn off the gnuplot capabilities
   ```-DCMAKE_C_COMPILER=gcc``` : Use a specific installation of your C compiler 
   ```-DCMAKE_CXX_COMPILER=g++``` : Use a specific installation of your C++ compiler 

### Quick start and first test ###
   0. You might review the details on the code in the Github wiki: https://github.com/OthmanB/TAMCMC-C/wiki (still under construction)
   Here I will just quickly state the important steps of an analysis, with specifics about running a quick on one of the provided models.

   1. Review the default configuration file ```[prog_root]/Config/default/config_default.cfg```. For a basic use, there is probably nothing to change there.

   2. Review and adapt the paths to the inputs files (*.model and *.data files). You will need to also provide the list of objects (ie, stars) that you want to process.
      See the provided config_presets.cfg as an template.

   3. Review the phase of the analysis (B, L, A). The template show a typical setup, but the number of samples is your choice and depends on the model complexity.
      The current template (for v1.4.0) is adequate for a test and for a local fit. L ~ 700 000 might be more suitable for a real analysis, with A ~ 50 000 on small parameter space. A ~ 500 000 - 2 000 000 may be required when fitting ~ 100 parameters. Consult the Github wiki to understand the parameters. But for a test, no need to touch anything except the paths to the inputs and output files. 

   3. Once compiled, the main file of the code (```cpptamcmc```) should be moved from the build directory, to the root directory of the program [prog_root].
      You need to execute it. The syntax for a call is reminded at execution of ./cpptamcmc
      for a test, just type: ```./cpptamcmc execute 1 1 1 1``` . This will execute the first model of the list on config_presets.cfg and start the analysis of the slice 1 in the model file (first range of frequency that was requested to analyse by the user).
      If you just want to proceed the analysis on the first slice (ignoring other frequency ranges), then ```./cpptamcmc execute 1 1 1 1 2``` (there is a slight bug here as it would be better to have 1 1 1 1 1 for consistency with past choice. It will be corrected later). 

   4. Once the whole analysis is finished, you will get a lot of files in ```[your_output_directory]/[the_name_of_the_analysed_object]```. The subdirectory ```diags``` contains plots and some diagnostics files. ```inputs_backup``` contains some zip of the configurations that was used. ```restore``` contains a snapshot of the MCMC chain state that can be used to restart an analysis. More importantly ```outputs``` contains the results of the analysis. 
      By default (also highly recommended), it is in a binary format. You can extract outputs from those files using the tools that were compiled at the same time than cpptamcmc. 
      The most important of those tool is certainly 'bin2txt' because it converts the binary list of MCMC samples for each parameters into plain ASCII files. Run ```./bin2txt``` to have a set of instruction on how to use it.
   
That's it, with that, you should have some fun playing around!
  
### Contribution guidelines ###

No external contribution is expected. This project is constantly improved, so please contact me if you need to see some worthy functionnality implemented. 

### Who do I talk to? ###

* Owner: Othman Benomar (Project Associate Professor at NAOJ / NYUAD Visiting Scientist)
* Contact: ob19@nyu.edu  /  othman.benomar@nao.ac.jp
