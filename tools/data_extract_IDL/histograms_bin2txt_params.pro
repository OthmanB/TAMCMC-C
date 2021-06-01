@nice_hist1D
@read_bin2txt_files
; A small program that plot the pdfs, after having converted binary files into text files (if requested)
; It returns the histograms and a summary of the statistical information (min_samples, med-2s, med-1s, med, med+1s, med+2s, max_samples, max_pdf)
function histograms_bin2txt_params, dir_out, Nb_classes, root_filename, dir_bin2txt, modelname, parameters_length, index0=index0, keep_period=keep_period, delete_txt_outputs=delete_txt_outputs, skip_plot=skip_plot

	tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

	if n_elements(index0) eq 0 then index0=0; 
	if n_elements(keep_period) eq 0 then begin ; Keep all samples by defaut
		keep_period = 1
	endif
	if n_elements(make_hist) eq 0 then make_hist=1
	if n_elements(delete_txt_outputs) eq 0 then delete_txt_outputs=1 ; by default, erase outputs in ASCII format... keep only the .SAV

	print, 'Converting the binary file into a serie of text file by using the built-in C++ function of TAMCMC...'
	;spawn, '../bin2txt/./bin2txt_params.out ' + root_filename + ' 0 ' + dir_out + ' ' + strtrim(round(index0), 2) + $
	;		' ' + strtrim(round(keep_period),2)  + ' 0' ; The last '1' means that we DO NOT replicate the constant in the ascii file
	spawn, dir_bin2txt + './bin2txt ' + root_filename + ' 0 ' + dir_out + ' ' + strtrim(round(index0), 2) + $
			' ' + strtrim(round(keep_period),2)  + ' 0' ; The last '1' means that we DO NOT replicate the constant in the ascii file

	files=file_search(dir_out + "*.ASCII") ; list all created files

	; ---- Convert and transfer the reduced samples into the IDL working directory (user-specified) ----
	Nsamples=ascii2sav(dir_out, delete_txt_outputs)
	files=file_search(dir_out + "*.sav") ; list all .sav files

	vals=evaluate_a1_inc(modelname, dir_out, 0) ; This function calculate and save a1 and inc from sqrt(a1).cosi, sqrt(a1).sini if required
	normalize=1

	; ----- All Pdfs and distributions information -----
	print, 'Generating the pdf from the samples...'
	i=0.
	N=0
	print, 'Determning the number of samples...'
	while (N lt 2) AND i lt n_elements(files) do begin ; Finding the number of parameters using first free parameter
		restore, files[i]
		N=n_elements(param)
		i=i+1.
	endwhile	
	if i ge n_elements(files) then begin
		print, 'No free parameters found!' 
		print, 'Cannot proceed further. The program will stop now'
		stop
	endif
	
	struc=plot_pdf(files, N, Nb_classes, dir_out)

	;hist_empty=build_histogram(randomu(seed, N),Nb_classes)
	;hist_all=dblarr(n_elements(files), 2, n_elements(hist_empty[0,*])) ; all the pdfs for all variables
	;Stat_Synthese_unsorted=dblarr(8, n_elements(files)) 
	;for i=long(0), n_elements(files)-1 do begin
	;	print, 'processing file: ' + files[i] + '...'
	;	restore, files[i]
	;	print, '   - build_histogram ...'
	;	if n_elements(param) eq N then begin
	;		hist=build_histogram(param,Nb_classes)
	;		hist_all[i, *,*]=hist
	;		distrib_stats=estimate_1_sigma_error( hist[0,*],hist[1,*],68.3,2,tab_critere)
	;		Stat_synthese_unsorted[0:6,i]=distrib_stats[3:*] ; les deciles, quartiles,mediane,...
	;		Stat_synthese_unsorted[7,i]=distrib_stats[0] ; le max
	;	endif 
	;	if n_elements(param) eq 1 then begin ; Case of a constant... we generate a pseudo histogram around the constant value
	;		hist=build_histogram(replicate(param, N),Nb_classes)
	;		hist_all[i, *,*]=hist
	;		;distrib_stats=estimate_1_sigma_error( hist[0,*],hist[1,*],68.3,2,tab_critere)
	;		Stat_synthese_unsorted[*,i]=param[0] ;distrib_stats[3:*] ; les deciles, quartiles,mediane,...
	;		;Stat_synthese_unsorted[7,i]=param[0];distrib_stats[0] ; le max
	;		;stop
	;	endif
	;	if n_elements(param) ne 1 AND n_elements(param) ne N then begin
	;		print, 'Incorrect number of samples in file: ', files[i]
	;		print, 'Cannot proceed as the histogram is of fixed size'
	;		print, 'Debug required'
	;		print, 'The program will stop now'
	;		stop
	;	endif
	;
	;	b=byte(files[i])
	;	posslash=max(where((b eq 47) OR (b eq 92))) ; detect last '/' OR '\'
	;	posdot=max(where(b eq 46)) ; detect last dot
	;	file_core=strtrim(b[posslash+1:posdot-1],2)
	;	if skip_plot eq 0 then begin
	;		nimp,name=dir_out + file_core +'.eps',/paper,/eps ; Savita's code
	;			nice_hist1D, hist, ps=1, file_out=file_out, xr=xr, $
	;				title=title, xtitle=variable_name[0], col=col, show_stats=1, $
	;				legendsize=legendsize, legend_precision=legend_precision
	;		fimp
	;	endif
	;	;if n_elements(param) eq 1 then stop ; THIS IS A DEBUG STOP ONLY. REMOVE IF YOU DO NOT KNOW THE FUNCTION OF THIS STOP
	;endfor

;struc={stat_synthese_unsorted:dblarr(8, n_elements(files)), hists:dblarr(n_elements(files), 2, N/Nb_classes), Nb_classes:0., Nsamples:0.}
;struc.stat_synthese_unsorted=stat_synthese_unsorted
;struc.hists=hist_all
;struc.Nb_classes=Nb_classes
;struc.Nsamples=Nsamples

return, struc
end

; This function can be used after histograms_bin2txt_params() in order to reprocess all 
; of the data by isolating multiples solutions that may exist for one parameter.
; dir_files_samples: The directory that contains the sav files created by histograms_bin2txt_params()
; ind: The index pointing to the variable that requires filtering. e.g if you want to isolate a solution in DP~220-250 from another at DP~500-550
;        You will need to know in advance in which file DP is stored, hence knowing its index
; ranges: A 2D table that contain a list of low/high boundaries that defines the different zone that are whished to be isolated
; dir_out: A root directory that will contain the ouputs. Note that one subdirectory per range will be created, with the following syntax:
;          [index]_[range_min]-[range_max]
;          The directory will also have a single file with the probability associated to the values within that range. For that it will perform the
;		   ratio of all counted values within the range by the total number of samples for all ranges
; Returns: a list of structures 
function filter_multimodalities, dir_files_samples, ind, ranges, dir_out, Nb_classes
	; Restoring the variable that will be used for filtering
	restore, dir_files_samples  + format_filename(ind, 0) +'.sav'
	var=param
	Nsamples=n_elements(var) ; Total number of samples. Required to evaluate the probabilities for the models

	sum_file='summary.sav'

	files=file_search(dir_files_samples + "*.sav") ; list all .sav files that contain all of the samples
	for r=0, n_elements(ranges[*,0])-1 do begin ; Get the indexes of the samples that need to be kept
		pos_r=where(var ge ranges[r, 0] AND var le ranges[r,1])
		subdir=strtrim(round(ind),2) + '_' + strtrim(string(ranges[r,0], format='(f6.2)'),2) + '-' + strtrim(string(ranges[r,1], format='(f10.2)'),2)
		f=file_search(dir_files_samples+subdir)
		if f[0] eq '' then spawn, 'mkdir ' + dir_files_samples+subdir
		print, 'pos_r = ', n_elements(pos_r)
		print, 'Nb_classes =', Nb_classes
		if n_elements(pos_r) ge 2*Nb_classes then begin ; Proceed to saving and plotting only if we have enough data point
			Proba=1.*n_elements(pos_r)/Nsamples
			;stop
			; Generate new sav files that are limited to a single solution
			for i=0, n_elements(files)-1 do begin
				restore, dir_files_samples  + format_filename(i, 0) +'.sav'
				param=param[pos_r]
				; Identifying the name of the file for the variable
				b=byte(files[i])
				posslash=max(where((b eq 47) OR (b eq 92))) ; detect last '/' OR '\'
				posdot=max(where(b eq 46)) ; detect last dot
				file_core=strtrim(b[posslash+1:posdot-1],2)
				save, variable_name, param, index, filename=dir_out+subdir+'/'+file_core + '.sav'
			endfor
			; Perform all of the plots for a given solution
			files_r=file_search(dir_out+subdir + '/*.sav')
			files_r=files_r[where(files_r ne dir_out + subdir + '/' + sum_file)]
			struc=plot_pdf(files_r, n_elements(pos_r), Nb_classes, dir_out+subdir + '/')
			save, struc, Proba, filename=dir_out+subdir+'/' + sum_file
			if r eq 0 then begin
				all_stat_synthese=dblarr(n_elements(ranges[*,0]), n_elements(struc.stat_synthese_unsorted[*,0]), n_elements(struc.stat_synthese_unsorted[0,*]) )
			endif
			all_stat_synthese[r,*,*]=struc.stat_synthese_unsorted
		endif else begin
			if pos_r[0] ne -1 then begin
				Proba=1.*n_elements(pos_r)/Nsamples
			endif else begin
				Proba=0
			endelse
			openw, 3, dir_out+subdir +'/error.txt'
				printf,3, '# Warning: No enough valid value found in the specified range'
				printf,3, '# range =' , strtrim(ranges[r,0], 2) + '   -   ' + strtrim(ranges[r,1],2)
				printf,3, '# Probability for that range is ' + strtrim(proba,2)
				printf,3, '# If you really wants plots and sav files, please consider reducing Nb_classes to a value below the number of data points'
				printf,3, '# Nsamples after filtering =', strtrim(n_elements(pos_r),2)
				printf,3, '# Current value for Nb_classes =', strtrim(Nb_classes,2)
			close,3
		endelse
	endfor
	return, all_stat_synthese
end

; Function that performs the plotting of the pdfs
; files: The list of sav files that will be read and that contain all of the samples
; N: The total maximum number of samples in each files
; Nb_classes: The number of classes used when plotting the pdf
; tab_critere: The talbe that controls the confidence intervals that will be calculated and shown in the figure
; dir_out : The output directory for the plots in eps format.
function plot_pdf, files, N, Nb_classes, dir_out
	tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma
	hist_empty=build_histogram(randomu(seed, N),Nb_classes)
	hist_all=dblarr(n_elements(files), 2, n_elements(hist_empty[0,*])) ; all the pdfs for all variables
	Stat_Synthese_unsorted=dblarr(8, n_elements(files)) 
	for i=long(0), n_elements(files)-1 do begin
		print, 'processing file: ' + files[i] + '...'
		restore, files[i]
		print, '   - build_histogram ...'
		cte_is_true=0 ; This is a required check to verify if the parameter is a constant in some cases
		if n_elements(param) gt 1 then if stddev(param) eq 0 then cte_is_true=1
		if n_elements(param) eq N AND cte_is_true eq 0 then begin
			hist=build_histogram(param,Nb_classes)
			hist_all[i, *,*]=hist
			distrib_stats=estimate_1_sigma_error( hist[0,*],hist[1,*],68.3,2,tab_critere)
			Stat_synthese_unsorted[0:6,i]=distrib_stats[3:*] ; les deciles, quartiles,mediane,...
			Stat_synthese_unsorted[7,i]=distrib_stats[0] ; le max
		endif 
		if n_elements(param) eq 1 or cte_is_true eq 1 then begin ; Case of a constant... we generate a pseudo histogram around the constant value
			hist=build_histogram(replicate(param[0], N),Nb_classes)
			hist_all[i, *,*]=hist
			;distrib_stats=estimate_1_sigma_error( hist[0,*],hist[1,*],68.3,2,tab_critere)
			Stat_synthese_unsorted[*,i]=param[0] ;distrib_stats[3:*] ; les deciles, quartiles,mediane,...
			;Stat_synthese_unsorted[7,i]=param[0];distrib_stats[0] ; le max
			;stop
		endif
		if n_elements(param) ne 1 AND n_elements(param) ne N then begin
			print, 'Incorrect number of samples in file: ', files[i]
			print, 'Cannot proceed as the histogram is of fixed size'
			print, 'Debug required'
			print, 'The program will stop now'
			stop
		endif

		b=byte(files[i])
		posslash=max(where((b eq 47) OR (b eq 92))) ; detect last '/' OR '\'
		posdot=max(where(b eq 46)) ; detect last dot
		file_core=strtrim(b[posslash+1:posdot-1],2)

		nimp,name=dir_out + file_core +'.eps',/paper,/eps ; Savita's code
			nice_hist1D, hist, ps=1, file_out=file_out, xr=xr, $
				title=title, xtitle=variable_name[0], col=col, show_stats=1, $
				legendsize=legendsize, legend_precision=legend_precision
		fimp
		;if n_elements(param) eq 1 then stop ; THIS IS A DEBUG STOP ONLY. REMOVE IF YOU DO NOT KNOW THE FUNCTION OF THIS STOP
	endfor
	struc={stat_synthese_unsorted:dblarr(8, n_elements(files)), hists:dblarr(n_elements(files), 2, N/Nb_classes), Nb_classes:0., Nsamples:0.}
	struc.stat_synthese_unsorted=stat_synthese_unsorted
	struc.hists=hist_all
	struc.Nb_classes=Nb_classes
	struc.Nsamples=N
return, struc
end

; Function that evaluates a1 and inc, if the model is not explicitly fitting them
; This is required to ensure that MS_Global_rotinc_correlations works as this program expect a1 (and possibly inc) to be variables
function evaluate_a1_inc, modelname, dir_files_samples, idl_format
	
parameters_length=read_plength(dir_files_samples +'/plength.txt')

Nmax=parameters_length[0]
lmax=parameters_length[1] ; number of visibilities
if idl_format eq 0 then begin
	Nf=total(parameters_length[2:5]) ; The total number of parameter is total(Nf_ls)
	k=parameters_length[6]
	l=parameters_length[7]
	mm=parameters_length[8]
endif else begin
	Nf=parameters_length[2] ; The total number of parameter is total(Nf_ls)
	k=parameters_length[3]
	l=parameters_length[4]
	mm=parameters_length[5]
endelse

i_a1=Nmax+Nf+lmax
restore, dir_files_samples  + format_filename(i_a1, idl_format) +'.sav'
a1=param
varname_a1=variable_name
i_inc=Nmax+Nf+k+l+mm+lmax
restore, dir_files_samples  + format_filename(i_inc, idl_format) +'.sav'
inc=param
varname_inc=variable_name

if n_elements(a1) eq 1 then a1=replicate(a1, 10) ; This is to avoid the crash of the stddev function
if n_elements(inc) eq 1 then inc=replicate(inc, 10) ; This is to avoid the crash of the stddev function

if stddev(a1) eq 0 and stddev(inc) eq 0 and modelname eq 'model_MS_Global_a1etaa3_HarveyLike' then begin ; We need to compute a1 and inc using sqrt(a1).cos(i) and sqrt(a1).sin(i)
	restore, dir_files_samples  + format_filename(i_a1+3, idl_format) +'.sav'	
	projcosi=param
	restore, dir_files_samples  + format_filename(i_a1+4, idl_format) +'.sav'	
	projsini=param
	
	if n_elements(projcosi) eq 1 then a1=replicate(projcosi, 10) ; This is to avoid the crash of the stddev function
	if n_elements(projsini) eq 1 then inc=replicate(projsini, 10) ; This is to avoid the crash of the stddev function

	if stddev(projcosi) ne 0 AND stddev(projsini) ne 0 then begin
		print, '   model=' + modelname + '...'
		print, '      ===> Convert sqrt(a1).cos(i) and sqrt(a1).sin(i) to a1 and inc...'
		a1=projcosi^2 + projsini^2
		inc=atan(projsini/projcosi)*180./!pi
		print, '      ===> Updating sav files for a1 and inc'
		param=a1
		variable_name=varname_a1
		save, filename=dir_files_samples  + format_filename(i_a1, idl_format) +'.sav', param, variable_name
		param=inc
		variable_name=varname_inc
		save, filename=dir_files_samples  + format_filename(i_inc, idl_format) +'.sav', param, variable_name
		
		;stop
	endif else begin
		print, 'Something wrong: There is no free inclination or splitting parameters'
		print, 'Cannot compute anything. Debug required'
		print, 'The program will stop now'
		stop
	endelse
	out=dblarr(2, n_elements(param))
	out[0,*]=a1
	out[1,*]=inc
	return, out;
endif else begin
	return, 0;
endelse

	;stop
end

