@ascii2sav.pro
@build_histogram_v2
@histograms_bin2txt_params.pro
@estimate_1_sigma_error
@MS_Global_fitplot
@Gaussian_fitplot
@Local_fitplot
@Echelle_Diagram
@fimp
@nimp
@fsc_color
pro showbestfit

	c=file_search('*.pro', /FULLY_QUALIFY_PATH ) ; It is assumed that the pro files are there
	b=byte(c[0])
	p=max(where(b eq 47))
	cpath=strtrim(b[0:p-1],2)
    ;dir_OS=cpath + '/../../' ; USE TO GET RESULTS DIRECTLY INTO THE PRODUCTS DIRECTORY INSIDE THE TAMCMC PROGRAM
	dir_OS='/Users/obenomar/Work/trash/novalue/aj_ln/'
    
    dir_outputs=dir_OS + '/outputs/'
    dir_inputs=dir_OS + '/inputs/'
    dir_out=dir_OS + '/products/'
	modelname='model_MS_Global_aj_HarveyLike'
	phase_list='A*'
	index0_list=0;
	;last_index=-1;
    ;modelname='model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3' ; MODEL FOR DEEPINVEST
    ;modelname='model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3'

;    dir_OS='/Users/obenomar/tmp/Simulator/'
;    dir_outputs=dir_OS + 'data/output_gaussfit/'
;    dir_inputs=dir_OS + 'data/input_gaussfit/'
;    dir_out=dir_OS + 'data/products_gaussfit/'
;    phase_list='L*'
;    index0_list=10000
;    modelname='model_Harvey_Gaussian'
;
;    dir_outputs=dir_OS + 'data/output_gaussfit_S1/'
;    dir_inputs=dir_OS + 'data/input_gaussfit_S1/'
;    dir_out=dir_OS + 'data/products_gaussfit_S1/'

;	dir_outputs=dir_OS + 'data/output_gaussfit_30-Jun-2021/'
;    dir_inputs=dir_OS + 'data/input_gaussfit_30-Jun-2021/'
;    dir_out=dir_OS + 'data/products_gaussfit_30-Jun-2021/'
;	dir_outputs=dir_OS + 'data/output_gaussfit_30-Jun-2021_fulldata/'
;    dir_inputs=dir_OS + 'data/input_gaussfit_30-Jun-2021_fulldata/'
;    dir_out=dir_OS + 'data/products_gaussfit_30-Jun-2021_fulldata/'

;	dir_outputs=dir_OS + 'data/output_gaussfit_01-Jul-2021/'
;    dir_inputs=dir_OS + 'data/input_gaussfit_01-Jul-2021/'
;    dir_out=dir_OS + 'data/products_gaussfit_01-Jul-2021/'

    ;modelname='model_RGB_asympt_a1etaa3_AppWidth_HarveyLike'
    ;dir_OS='/Volumes/home/' ; On the mac 
    ;dir_OS='~/obnas-shared/dataonly/' ; On the virtual machine    
    ; ##### 1 #####
    ;modelname='model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3'
    ;dir_outputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth/test/outputs/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;dir_inputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth/test/inputs/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;dir_out=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth/test/products/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;index0_list=[170000, 0, 0 , 0 , 0, 5000, 0, 0, 10000, 15000, 300000, 15000, 0, 250000]
    ;phase_list=['B0*', 'B0*', 'B0*', 'B0*', 'B0*', 'B1*', 'B0*', 'B0*', 'B1*', 'B1*', 'B0*', 'B1*', 'B0*', 'B0*']
    ;dir_out = cpath + '/products/'
    
    ; ##### 2 #####
    ;modelname='model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3'
    ;dir_outputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth-2/test/outputs/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;dir_inputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth-2/test/inputs/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;dir_out=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth-2/test/products/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;index0_list=[10000, 10000, 40000, 75000, 15000, 15000, 45000]
    ;phase_list=replicate('B0*', 7)
    
    ; ##### 3 #####
    ;modelname='model_RGB_asympt_a1etaa3_freeWidth_HarveyLike_v3'
    ;dir_outputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.57-Siddarth-freeWidth/test/outputs/Siddarth_Mail-Dec20-INCOMPLETED-freeWidth/'
    ;dir_inputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.57-Siddarth-freeWidth/test/inputs/Siddarth_Mail-Dec20-INCOMPLETED-freeWidth/'
    ;dir_out=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.57-Siddarth-freeWidth/test/products/Siddarth_Mail-Dec20-INCOMPLETED-freeWidth/'
    ;phase_list=['B0*','B0*','B0*','B0*','B0*','B0*','B0*','B0*','B0*','B0*','B0*','B0*','B0*','B0*']   ; check status.log in the outputs directory to a list of the status
    ;index0_list=[120000, 0, 0, 30000, 30000, 340000, 10000, 0, 350000, 375000, 300000, 350000, 15000,375000] 

	;; ##### 4 #####
    ;modelname='model_RGB_asympt_a1etaa3_AppWidth_HarveyLike_v3'
    ;dir_outputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth-001867706_deepinvest/test/outputs/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;dir_inputs=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth-001867706_deepinvest/test/inputs/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;dir_out=dir_OS + '2020/Mixed-modes/Dalma_processes/TAMCMC-C-1.56-Siddarth-001867706_deepinvest/test/products/Siddarth_Mail-Dec20-INCOMPLETED/'
    ;phase_list=['B0*','B1*']   ; check status.log in the outputs directory to a list of the status
    ;index0_list=[250000, 55000] 

    ;phase='A*'
    ;phase='L*'
 
    dir_filter='*' ; Used to choose which directory should be processed
    ;dir_filter='kplr010963065*'

    Nb_classes=100.;
    ;index0=70000. ; index of the first entry that is kept
    keep_period=4. ; Keep 1 out of keep_period
 	
    dirs=file_search(dir_outputs + dir_filter) ; lists the directories
    Ndirs=n_elements(dirs)
    Ok=intarr(Ndirs)
    
    i0=0
	imax=Ndirs-1
    for i=long(0), Ndirs-1 do begin
    	print, '[' + strtrim(i,2) + '] ' + dirs[i]
    endfor
    print, "type '.cont' to proceed"
    stop
    if n_elements(phase_list) eq 1 then phase_list=replicate(phase_list[0], Ndirs)
  	if n_elements(index0_list) eq 1 then index0_list=replicate(index0_list[0], Ndirs)  
    for i=long(i0), long(imax) do begin
    	phase=phase_list[i]
    	index0=index0_list[i]
    	print, '    ------- Processing ' + dirs[i] + ' --------'    
    	b=byte(dirs[i])
    	pos=max(where(b eq 47 OR b eq 92)) ; detect slashes
        if pos[0] eq -1 then starID=dirs[i] else starID=strtrim(b[pos+1:*], 2)
    	done=PostMCMC_showbestfit(dir_outputs, dir_inputs, dir_out, modelname, starID, phase, Nb_classes, index0, keep_period)
    	print, '    -------------------------------------------'
		OK[i]=done
	endfor
	posNotOK=where(OK eq 0)
	if n_elements(posNotOK) gt 1 then begin
		print, "List of process that could not be processed:"
		for i=0, n_elements(posNotOK)-1 do print, dirs[posNotOK[i]]
	endif else begin
		if posNotOK[0] eq -1 then print, 'All done'
		if posNotOK[0] ne -1 then begin
			print, "List of process that could not be processed:"
			print, dirs[posNotOK[0]]
		endif
	endelse
end

; Used to interpret the results from all MS_Global models
; starID: the id as define in the config_presets.cfg
; phase: phase (B, L, A) as defined in the config_presets.cfg
; Nb_classes: Number of classes for the histograms (default = 200)
; index0: first index for the kept samples (default = 0)
; keep_period: periodicity for the the kept samples (default = 1  ==> all samples are kept)
; delete_txt_outputs: if 1, deletes the ascii files for the samples output files. Only sav files will be kepts (default = 0)
function PostMCMC_showbestfit, root_outputs, root_inputs, dir_out, modelname, starID, phase, Nb_classes, index0, keep_period
	done=1

	subdir=''	
	;dir_bin2txt='../' ; directory where the function that converts binaries into ascii is.
	;dir_getmodel='../' ; directory where the function that generate the models is
	dir_bin2txt='cpp_prg/' ; directory where the function that converts binaries into ascii is.
	dir_getmodel='cpp_prg/' ; directory where the function that generate the models is

	; --- Defining the directory/files using the strict rule for managing inputs/outputs ----
	dir_IDL_out=dir_out + starID +'/'
	binresultdir=root_outputs + starID + '/outputs/'
	diagsdir=root_outputs  + starID + '/diags/'
	root_filename=binresultdir + starID + '_' + phase + '_params' ; here phase may contain and asterix
	test=file_search(root_filename + '_chain-0.bin')
	if test eq '' then begin
		print, 'Could not find files compatible with requested phase'
		print, 'Change the phase name. The program will stop now'
		done=0
		;goto, bypass
		;stop
		goto, bypass
	endif else begin
		root_filename=detect_root_name(test) ; ensure that we use a valid root_filename
	endelse
	;stop
	data_file=file_search(root_inputs + starID + '.data' )
	if data_file eq '' then begin
			print, 'Warning: Input data file not found.' 
			print, 'Check that the given directory and file exist'
			print, 'The program will stop now'
			stop
	endif
	
	check_dir=file_search(binresultdir) ; check if the dir has no capital letter
	if check_dir[0] eq '' then begin
			print, 'Warning: Output data files not found in the specified directory.' 
			print, 'Check that you provided the correct output directory'
			print, 'The program will stop now'
			stop
	endif
	
	f=file_search(dir_IDL_out + '')
	if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + ''

	; ----- Convert binary files in text files, easier to read by IDL AND plot their pdf using nice1D_hist ----
	skip_plot=0 ; To accelerate debuging, put this to 1
	print, 'Convert binary files in text files and then into sav files, easier to read by IDL and plot their pdf...'
	f=file_search(dir_IDL_out + 'Files/')
	if f[0] eq '' then spawn, 'mkdir ' + dir_IDL_out + 'Files/'	
	hists_info=histograms_bin2txt_params(dir_IDL_out+'Files/',  Nb_classes, root_filename, dir_bin2txt,  modelname,index0=index0, keep_period=keep_period, skip_plot=skip_plot)
	parameters_length=read_plength(dir_IDL_out +'Files/plength.txt')

	;; ---------------------------- MANUAL FILTERING ---------------------------------
	;; ------ WARNING TOUCH THIS ONLY IF YOU WANT TO ISOLATE SEPARATE SOLUTIONS ------
	;; --- IF YOU DON't UNDERSTAND WHAT THIS IS, COMMENT IT TO AVOID BAD SURPRISES ---
	;index=15 ; DP for Warrick binary star, spectrum 00000004
	;; We isolate 4 different solutions: DP~[150,200], [220, 280], [300, 350], [500, 550]
	;ranges=dblarr(4, 2)
	;ranges[0,*]=[150,200]
	;ranges[1,*]=[220,280] ; This range is somewhat empty
	;ranges[2,*]=[300,350]
	;ranges[3,*]=[500,550]
	;all_synthese=filter_multimodalities(dir_IDL_out+'Files/', index, ranges, dir_IDL_out+'Files/', Nb_classes)
	;sig_pos=[19,20] ; positions where there is variance hyperparameters
	;all_synthese[*,*,sig_pos[0]]=0 ; forces variance to 0 to avoid wrong plots
	;all_synthese[*,*,sig_pos[1]]=0
	;ind_plot=0 ; Define which solution this function will plot
	;hists_info={stat_synthese_unsorted:reform(all_synthese[ind_plot, *,*])}
	;;stop
	;; -------------------------------------------------------------------------------
	;; -------------------------------------------------------------------------------
	
	; ---- Save the model and the median parameters----
	print, 'Modeled spectrum and median values...'	
	val_med=reform(hists_info.stat_synthese_unsorted[3, *])
	;val_med_m1s=val_med ; bypassing my attempt to show uncertainty because the solution I had implemented was not working well
	;val_med_p1s=val_med ; bypassing my attempt to show uncertainty because the solution I had implemented was not working well
	;val_med_m1s=reform(hists_info.stat_synthese_unsorted[2, *]) ; median - 1sigma ; THIS DOES NOT WORK DUE TO CORRELATIONS IN THE HYPERSPACE
	;val_med_p1s=reform(hists_info.stat_synthese_unsorted[4, *]) ; median + 1sigma ; THIS DOES NOT WORK DUE TO CORRELATIONS IN THE HYPERSPACE
	; Write a configuration file suitable for the compute_model.cpp function
	params_cfg= dir_IDL_out + 'best_models_params.txt'
	file_out=dir_IDL_out + 'best_models_fit.ascii'
	file_psfit=dir_IDL_out + 'best_models_fit.eps'
	openw, 3, params_cfg
		str='# This file contains in the first line, the parameters structure (plength). All following lines, correspond to a single vector of parameters'
		printf, 3, str
		str=''
		for i=long(0), n_elements(parameters_length)-1 do str=str+ '   ' + strtrim(parameters_length[i],2)
		printf, 3, str
		str=''
		for i=long(0), n_elements(val_med)-1 do str=str+ '   ' + strtrim(val_med[i],2)
		printf, 3, str
		;str=''
		;for i=long(0), n_elements(val_med_m1s)-1 do str=str+ '   ' + strtrim(val_med_m1s[i],2)
		;printf, 3, str
		;str=''
		;for i=long(0), n_elements(val_med_p1s)-1 do str=str+ '   ' + strtrim(val_med_p1s[i],2)
		;printf, 3, str
	close, 3
	; Use the in-built function of TAMCMC to get the median model
	print,' Executing: ./getmodel ' + data_file + ' ' +  params_cfg + ' ' + modelname + ' ' + file_out
	spawn, dir_getmodel +'./getmodel ' + data_file + ' ' +  params_cfg + ' ' + modelname + ' ' + file_out
	f=file_search('params*.model') ; File generated by getmodel into the execution directory since v1.61 ==> Allow to have a direct output formated table of the modes
	if f[0] ne '' then begin
		print, '   Note: getmodel version > 1.61: Formated outputs from the MCMC available in the output directory'
		for i=0, n_elements(f)-1 do spawn, 'mv ' + f[i] + ' ' + dir_IDL_out + '/best_fit_params_' + strtrim(i,2) + '.model'  
	endif	
	pos_fl0=total(parameters_length[0:1])
	Nfl0=parameters_length[2]
	; read the file that was just created
	Ncols=detect_Ncolumns(file_out, skip=0)
	model_bestfit=read_Ncolumns(file_out, Ncols, 5d5, skip=0, ref_N=0)
	fmin=min(model_bestfit[0,*])
	fmax=max(model_bestfit[0,*])	
	if modelname ne 'model_Harvey_Gaussian'  then begin
		if modelname ne 'model_MS_local_basic' then begin
			fit=linfit(findgen(Nfl0), val_med[pos_fl0:pos_fl0+Nfl0-1])
		endif else begin
			fit=dblarr(2)
			; We need to get an actual range that correspond to the actual frequencies that are fitted
			range=guess_localfit_freqrange(params_cfg, parameters_length)
			fmin_loc=range[0]
			fmax_loc=range[1]
		endelse
		if fit[0] ne 0 then begin ; Dealing with a global fit
			;stop
			MS_Global_fitplot, model_bestfit[0, *], model_bestfit[1, *], model_bestfit[0, *], model_bestfit[2, *], fit[1], file_psfit, fmin, fmax
		endif else begin ; Dealing with a local fit
			Local_fitplot, model_bestfit[0, *], model_bestfit[1, *], model_bestfit[0, *], model_bestfit[2, *], file_psfit, fmin, fmax, fmin_loc, fmax_loc
		endelse
	endif else begin
			Gaussian_fitplot, model_bestfit[0, *], model_bestfit[1, *], model_bestfit[0, *], model_bestfit[2, *], file_psfit, fmin, fmax
	endelse
	model_noise=do_residuals(dir_IDL_out, dir_getmodel, data_file, modelname, val_med, parameters_length)
	if modelname eq 'model_Harvey_Gaussian' then begin
		f=file_search(dir_IDL_out + '/NewGaussfit_Guess')
		if f[0] eq '' then spawn, 'mkdir ' + dir_out+'/NewGaussfit_Guess'
		file_out=dir_out+'/NewGaussfit_Guess/'+starID+'.guess.txt'
		guess_numax_from_noisefit, model_noise[0,*], model_noise[1,*], model_noise[2,*], val_med, file_out, [10, max(model_noise[0,*])]
	endif
	if modelname ne 'model_Harvey_Gaussian' then begin
		; ----- Echelle diagram --s-	
		file_ps_echelle=dir_IDL_out + 'Echelle_Diagram.eps'  
		nimp,name=file_ps_echelle,/paper,/eps
			; The input here must be : freq, spectrum , noise model. NOTE: THE NOISE MODEL IS ALWAYS REMOVED INSIDE show_ech_diag_CPP
			;show_ech_diag_CPP, dir_IDL_out + '/best_fit_params_0.model', model_noise, ps=1, shifts=shifts, trunc_spec=trunc_spec
			show_ech_diag_CPP, dir_IDL_out + '/params_0.model', model_noise, ps=1, shifts=shifts, trunc_spec=trunc_spec
		fimp
		file_ps_echelle=dir_IDL_out + 'Echelle_Diagram_residuals.eps' 
		nimp,name=file_ps_echelle,/paper,/eps
			; The input here must be : freq, spectrum , noise model. NOTE: THE NOISE MODEL IS ALWAYS REMOVED INSIDE show_ech_diag_CPP
			;show_ech_diag_CPP, dir_IDL_out + '/best_fit_params_0.model', model_bestfit[0:2,*], ps=1, shifts=shifts, trunc_spec=trunc_spec
			show_ech_diag_CPP, dir_IDL_out + '/params_0.model', model_bestfit[0:2,*], ps=1, shifts=shifts, trunc_spec=trunc_spec
		fimp
		save, model_bestfit, model_noise, filename=dir_IDL_out + '/spectrum_models.sav'
		;stop
	endif

	bypass:
	return, done
end

; A function that take the best fit (with or without gaussian)
; and evaluate the position of numax using the residuals
; This allows you to have a new parameter vector for a subsequent
; fit that would have good guesses for the noise and therefore
; a more accurate guess for numax than what can be done with init_fit.py or init_fit.pro
; This is therefore very usefull when a fit failed to detect the gaussian while it is evident
; for a human that there is one.
; This will however not work if there is 'polution' such as peaks due to binaries
; because the guesses are made based on the max of the residual
pro guess_numax_from_noisefit, freq, spec, model, params_fit, file_out, xrange
	noise_params=params_fit
	noise_params[7]=0

	pos=where(freq ge xrange[0] AND freq le xrange[1])
	y=spec[pos]
	x=freq[pos]
	scoef=0.75/(freq[2]-freq[1])  ; Smooth over 3 microHz
	res=smooth(spec, scoef, /edge_truncate)/model
	res=res[pos]
	posmax=where(res eq max(res))
	numax_guess=x[posmax]
	Amax_guess=res[posmax]

	xmin=0
	xmax=max(x)
	ymin=0
	ymax=max(res)*1.2
	nimp,name=file_out+'.eps',/paper,/eps
	plot, x, res, /nodata, background=fsc_color('White'), color=fsc_color('Black'), thick=2, $
			xr=[xmin,xmax], yr=[ymin,ymax], /xst, /yst
	oplot,x, res, color=fsc_color('Dark Gray')
	oplot,x, replicate(1, n_elements(x)), linestyle=2, thick=3, color=fsc_color('Black')
	oplot,[numax_guess,numax_guess], [xmin, xmax], color=fsc_color('red'), thick=2, linestyle=2
	fimp
	guess_params=noise_params
	guess_params[7]=Amax_guess
	guess_params[8]=numax_guess
	openw, 3, file_out
		str='# This file contains a single line that gives results from guess_numax_from_noisefit(), function that provides new guesses for a Gaussian fit using a previous fit'
		printf, 3, str
		str=''
		for i=long(0), n_elements(guess_params)-1 do str=str+ '   ' + strtrim(guess_params[i],2)
		printf, 3, str
	close,3
	;stop
end

function do_residuals, dir_IDL_out, dir_getmodel, data_file, modelname, val_med, parameters_length
	params_cfg_noise= dir_IDL_out + 'noise_model_params.txt'
	file_noise_fit=dir_IDL_out + 'noise_models_fit.ascii'
	openw, 3, params_cfg_noise
		str='# This file contains in the first line, the parameters structure (plength). All following lines, correspond to a single vector of parameters'
		printf, 3, str
		str=''
		for i=long(0), n_elements(parameters_length)-1 do str=str+ '   ' + strtrim(parameters_length[i],2)
		printf, 3, str
		str=''
		val_med_noise=val_med
		val_med_noise[0:parameters_length[0]-1]=0 ; put heights to 0
		s=n_elements(parameters_length)
		if modelname ne 'model_Harvey_Gaussian' then begin
			if s gt 10 then val_med_noise[total(parameters_length[0:s-1])-1:*]=0 ; Whatever is after the inclination must be put to 0... this to deal with models that may have modes after inclination
		endif else begin
			val_med_noise[7]=0 ;Amplitude of the Gaussian set to 0
		endelse
		for i=long(0), n_elements(val_med_noise)-1 do str=str+ '   ' + strtrim(val_med_noise[i],2)
		printf, 3, str
	close, 3
	; Use the in-built function of TAMCMC to get the median model with noise only
	spawn, dir_getmodel +'./getmodel ' + data_file + ' ' +  params_cfg_noise + ' ' + modelname + ' ' + file_noise_fit
	f=file_search('params*.model') ; File generated by getmodel into the execution directory since v1.61 ==> Allow to have a direct output formated table of the modes
	if f[0] ne '' then begin
		print, '   Note: getmodel version > 1.61: Formated outputs from the MCMC available in the output directory'
		for i=0, n_elements(f)-1 do spawn, 'mv ' + f[i] + ' ' + dir_IDL_out + '/noise_fit_params_' + strtrim(i,2) + '.model'  
	endif
	; read the file that was just created
	;Ncols=5 ; col[0]=freq, col[1]=spec_reg, col[2]=median_model_noise
	Ncols=detect_Ncolumns(file_noise_fit, skip=0)
	model_noise=read_Ncolumns(file_noise_fit, Ncols, 5d5, skip=0, ref_N=0)

	return, model_noise
end

; Uses the information on the file fed into the getmodel (txt best fit file) 
; In order to evaluate an empirical range for the fitted frequencies
function guess_localfit_freqrange, params_cfg, plength
	low=''
	high=''
	med=''
	a=''
	openr, 3, params_cfg
		readf, 3, a ; read fist line ==> ignore as will be a comment
		readf, 3, a ; read second line ==> ignore as plength is known
		readf, 3, med ; read third line ==> median vector ignore
		readf, 3, low ; read fourth line ==> -1sigma vector ==> lower freq bound
		readf, 3, high ; read fifth line ==> +1sigma vector ==> uper freq bound
	close, 3
	lvec=split_txt(low)
	hvec=split_txt(high)

	; Extract all of the frequencies
	p0=total(plength[0:1])
	p1=total(plength[0:5])-1
	fmin=min(lvec[p0:p1])
	fmax=max(hvec[p0:p1])

	p0width=total(plength[0:6])
	wmax=max(hvec[p0width:p0width+plength[7]-1]) ; look at the maximum width to better evaluate fmin/fmax

	fmin=fmin-wmax*4
	fmax=fmax+wmax*4
	;stop
	return, [fmin,fmax]
end


; A short function that take a text line and split it assuming spaces as separator
function split_txt, line
	uu=strsplit(line)
    N_uu=N_elements(uu)-1
    vec=dblarr(N_uu+1)
    for j=0,N_uu-1 do begin
        vec(j)=float(strmid(line,uu(j),uu(j+1)-uu(j)-1))
    endfor
	vec(N_uu)= float(strmid(line,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))

	return, vec
end

; A short function that read the plength.txt files generated by bin2txt
function read_plength, file
	a=''
	i=0.
	openr, 3, file
		while (eof(3) ne LOGICAL_TRUE(1)) do begin
			 readf, 3, a
			 if i eq 0 then plength=double(a) else plength=[plength, double(a)]
			 i=i+1
		endwhile
	close,3
	return, plength
end

; Determine how many columns exists in a file
function detect_Ncolumns, file, skip=skip

if n_elements(skip) eq 0 then skip=1 ; defaut we skip one line only

openr, 3, file

	
	a=''
	i=0d
    while i le skip +1 do begin
      	if i lt skip then readf,3,format='(q)'
        if i ge skip then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)
        endif
        i=i+1
    endwhile
close, 3

	;stop
return, N_uu
end

; N: number of columns
; K: maximum number of lines
function read_Ncolumns, file,N, K, skip=skip, ref_N=ref_N

if n_elements(skip) eq 0 then skip=1 ; defaut we skip one line only
if n_elements(ref_N) eq 0 then ref_N=1 ; defaut we identify the zero non-used tab elements with column 1
if n_elements(spectrum) eq 0 then spectrum=1
openr, 3, file

	param=dblarr(K,N)
	a=''
	i=0d
      while EOF(3) ne 1 do begin
      	if i lt skip then readf,3,format='(q)'
        if i ge skip then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)-1
          	for j=0,N_uu-1 do begin
          		param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          	endfor
			param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		endif
		i=i+1
      endwhile

close,3
param0=param
test=where(param[*,ref_N] ne 0)

param=param[test,*]
param=transpose(param)

return,param
end


function file_syntax, in, core, add_extension=add_extension
	if n_elements(add_extension) eq 0 then add_extension=1

	if in lt 10 then file_out=core+'00'+strtrim(round(in),1)
	if in ge 10 AND in lt 100 then file_out=core+'0'+strtrim(round(in),1)
	if in ge 100 then file_out=core +strtrim(round(in),1)
	if add_extension eq 1 then begin
		file_out=file_out + '.sav'
	endif
	
return, file_out
end

function format_filename, in, idl_format
	;NOTE ON IDL_FORMAT: This is for the very old anaylsis (before 2014). It is discontinued
	; but for sake of compatibility, I had to but it here as a dummy variable
	return, file_syntax(in, '', add_extension=0)
end

function addzeros, in

	if in lt 10 then out='00'+strtrim(round(in),2)
	if in ge 10 AND in lt 100 then out='0'+strtrim(round(in),2)
	if in ge 100 then out=strtrim(round(in),2)
	
	return, out
end

; cut the file name just before _chain-*.bin
function detect_root_name, file_chain_bin

	b=byte(file_chain_bin)
	pos=max(where(b eq 95)) ;detect last '_'
	name=strtrim(b[0:pos-1],2)
	return, name
end

function interpret_varnames_Cpp, varname
	variable_name=strarr(4)
	variable_name[0]=varname
	return, variable_name
end
