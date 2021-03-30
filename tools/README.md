### Purpose of getmodel ###
	Read a simple input file that contains length and a list of vector of parameters. From this, it builds a model, provided that the user gave:
			- the model name
			- the data filename, in the same format as when given to TAMCMC
	The model is written in an output ASCII file. That output can be user-specified, or a default one


### Purpose of bin2txt ###
	Convert binary outputs of the program into ASCII outputs (one file per parameter/variables), including constants

### Purpose of read_stats ###
    Convert binary files with the Likelihood/Prior/Posterior into ASCII outputs.
    
### Purpose of data_extra_IDL ### 
    This is a suite of IDL functions that allows you to convert the binary outputs from the MCMC fitting results into a serie of sav files along with a pdf on the form of an eps image for each of the parameters that were considered for the MCMC fitting. It also show the best fit on the top of the data. This program calls bin2txt and getmodel that must be in one directory below the data_extra_IDL directory (directory that is created when unzipping the file). 
    This code is usefull for a first analysis of the outputs.
    
### Purpose of a2_nl.py ###
    To be used only for models that handle a2 coefficient using a polynomial function. This program translates the a2 polynomial coeficients extracted using data_extract_IDL.zip into a2(nu_nl) with a proper handling of the error propagation. It makes a jpg plot of the this.

### Purpose of nusini_nucosi_disociator.py ###
    To be used with models that fit sqrt(a1).sini and sqrt(a1).cosi: It will make those two correlated quantities separated : inc , a1
