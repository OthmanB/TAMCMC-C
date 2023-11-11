@read_ID_list_file
; very simple function that read a one column format that contains the list of kics
function getKIClist, listfile
	openr, 3, listfile
		a=''
		cpt=0
		kiclist=strarr(1)
		while eof(3) eq 0 do begin
			readf, 3, a
			if cpt eq 0 then kiclist[0]=strtrim(a,2) else kiclist=[kiclist, strtrim(a,2)]
			cpt=cpt+1
		endwhile
	close, 3
	return, kiclist
end

; basic function that given a prefix/suffix, identifies the kic number of stars
function guessKIClist, dir_files

	;prefix='kplr'
	;suffix='_COR_PSD_filt_inp'
	prefix=''
	suffix=''
	;prefix='LC_KOI_4.1_'
	;suffix='.tf'
	
	files=file_search(dir_files + prefix + '*' + suffix + '.sav') ; look for all stars with the suffix
	
	bpre=byte(prefix) & Nbpre=n_elements(bpre)
	bsuf=byte(suffix) & Nsuf=n_elements(bsuf)
	
	if suffix eq '' then Nsuf=0
	if prefix eq '' then Nbpre=0
	
	;stop
	KIClist=strarr(n_elements(files))
	for i=0, n_elements(files)-1 do begin
		fileonly=splitfilename(files[i])
		b=byte(fileonly[1])
		k=strtrim(b[Nbpre: n_elements(b)-1-Nsuf],2)
		print, '[' + strtrim(i,2) + ']  ' + k
		KIClist[i]=k
		;stop
	endfor
	;stop
	return,kiclist
end


; Date: 17/07/2015
; Program that Run envelope_measure and then create a proper star_list.txt file
; Such a star_list.txt file can be used as a configuration file for the Preset-Analysis
; Program. The star_list.txt file is created into the dir_files directory.
; dir_files: Directory in which we look for power spectra. ALL FILE THAT MATCH THE 
;			 EXPECTED SYNTAX ARE PROCESSED. The expecte name of the files should be in
; 			 the format: 'PSF*_'+KIC_number + '*.sav'
;			 It is also expected that each file contains two variables:
;			 		- freq : the frequency in microHz
;			 		- spec_reg : the power spectrum in ppm^2/microHz
; KIC_list: KIC number of all files that we want to process
; search_range: range in which we look for modes. Typical values are:
;			     	- [600, 4000] microHz for Main sequence stars
;			     	- [300, 1500] microHz for Subgiants
;			     	- [20 ,  400] microHz for Red Giants  
; suffix: put a suffix after the KIC
; prefix: put a prefix in front of the KIC
; numax_step: parameter that defines the search grid for numax. Not a critical parameter;
;			  But you need to be careful not to put a too large value. It will be faster
;			  but might be problematic for Red Giants. 
;			  Unless specified, numax_step is set to 10 microHz.
pro get_numax_Amax, dir_files, KIC_list, search_range, suffix, prefix, numax_step=numax_step

if n_elements(suffix) eq 0 then suffix=''
if n_elements(prefix) eq 0 then prefix=''
if n_elements(numax_step) eq 0 then numax_step=10

NKIC=n_elements(KIC_list)

print, '---------------------'
print, 'Searching for the files matching the requested KIC numbers...'
cpt=0.
cpt2=0
for i=0, NKIC-1 do begin
	print, ' - Looking for a file with given suffix/prefix (set to blank if not provided)...'
	file0=file_search(dir_files + prefix + '*' + strtrim(KIC_list[i],2) + suffix + '*.sav')
	;stop
	if file0 ne '' then begin ; If the file was found!
		print, '   File ' + file0 + ' found...'
		files=autofill(files, file0, cpt) ; create/increment a table. Creation occurs if cpt=0.
		KIC_found=autofill(KIC_found, strtrim(KIC_list[i],2), cpt)
		cpt=cpt+1.
	endif else begin ; If the file is not found...
		print, 'File with given suffix/prefix not found. Other possibilities will be tried...'
	
		print, ' - Looking for a file with the syntax: ' + 'PSF_' + strtrim(KIC_list[i],2) + '.sav'
		file0=file_search(dir_files + 'PSF_' + strtrim(KIC_list[i],2) + '.sav')
		;file0=file_search(dir_files  + strtrim(KIC_list[i],2) + '*.sav')
		if file0 ne '' then begin ; If the file was found!
			print, '   File ' + file0 + ' found...'
			files=autofill(files, file0, cpt) ; create/increment a table. Creation occurs if cpt=0.
			KIC_found=autofill(KIC_found, strtrim(KIC_list[i],2), cpt)
			cpt=cpt+1.
		endif else begin ; If the file is not found...
			print, '   No file found with the assumed syntax...'
			print, ' - Trying the syntax: ' + strtrim(KIC_list[i],2) + '.sav'
			file0=file_search(dir_files  + strtrim(KIC_list[i],2) + '.sav') ; Test another format... 
			if file0 ne '' then begin ; If the file was found...
				print, '   File ' + file0 + ' found!'
				files=autofill(files, file0, cpt) ; create/increment a table. Creation occurs if cpt=0.
				KIC_found=autofill(KIC_found, strtrim(KIC_list[i],2), cpt)
				cpt=cpt+1.
			endif else begin
				print, '-----------------------'
				print, 'Warning: No file found matching the search criteria for KIC=' + strtrim(KIC_list[i],2)
				print, 'Check whether your filename/KIC_number/search criteria is correct!' 
				print, 'The star will be ignored!'
				;print, '...The program will stop now'
				print, '-----------------------'
				if cpt2 gt 0 then rem_i=[rem_i, i] else rem_i=i
				cpt2=cpt2+1
			endelse
		endelse
	endelse	
endfor
Nfiles=n_elements(files)
print, 'Search finished. A total of ' + strtrim(Nfiles,2) + ' files over a total of ' + $
		strtrim(NKIC,2) + ' requested files'
print, '---------------------'


psfiles=files ; the extension '.ps' will be added by the write_on_ps_on function.
txtfile=dir_files + 'star_list.txt'

sfactor=6

;stop ; for debug

print, 'Proceeding to the numax and Amax determination...'
a=''
openw, 3, txtfile
a='# Next line is the directory where to look for spectra. Next follows the KIC/numax/e(numax)/Amax' 
printf,3, a
a=dir_files
printf,3, a
for i=0, Nfiles-1 do begin
	restore, files[i]
	res=envelope_measure(freq, spec_reg, search_range, numax_step, $
					sfactor=sfactor, psfile=psfiles[i])
	a=string(strtrim(KIC_found[i],2), format='(i15)') + ' ' + $  ; KIC number
	  string(strtrim(res[0],2), format='(f15.2)') + ' ' + $ ; numax
	  string(strtrim(res[1],2), format='(f15.2)') + ' ' + $ ; uncertainty on numax
	  string(strtrim(res[2],2), format='(f15.2)')           ; Amax
	printf, 3, a
endfor
close,3

print, 'Finished. A summary file -star_list.txt- was created in ' + $
		strtrim(dir_files,2)
print, 'The program will exit now'
end

; Date: 13/07/2015
; determine numax by computing a slice-averaged power spectrum and look for 
; the maximum power in each of the slices. The width of the slices is defined using the Dnu
; derived using the scaling relation with numax found in Stello et al. (2009).
; freq, spec_reg: the spectrum
; fold_range=[fmin, fmax]: On which range of frequency we look for a repetitive pattern
; sfactor: Defines the level of smooth applied on the normalised detection criteria
;          (normalised slice-averaged power spectrum) max(s)/smooth(s) for each slice.  
;           This reduce a bit the noise effect.
;          If not set, sfactor=1.
; nu_rebin: if set, the spectrum is rebined at a resolution nu_rebin before proceeding
;           to the analysis. This is useful if the spectra has a very high resolution
;            ie, much higher than the typical precision on Dnu.
;            Thus it is adequate to set to have a faster process. The rebin is done using 
;            congrid and the interpolation (/interp) option.
; psdir: If set, create a ps file with the plot showing the normalised slice-averaged 
;		  power spectrum. 
;		 When set, it must contain the directory on which we wish to create the ps file.
function envelope_measure, freq, spec_reg, fold_range, numax_step, $
					sfactor=sfactor, nu_rebin=nu_rebin, $
					psfile=psfile

if n_elements(sfactor) eq 0 then sfactor=1

; ---------- Constants ---------
Teff_sun=5777.
Dnu_sun=135.
numax_sun=3100.
a0=0.77 ; Stello et al. (2009) for the relation between Dnu and numax
; ---- Values for the relation between Gamma and numax (Benomar 2013, Eq.1)
alpha=double(0.267) & beta=double(0.76) & aa=4.06d-42 & po=double(11)
K=alpha^(8./3.) * numax_sun^2 * Teff_sun / Dnu_sun^(8./3.)
; ----------------- These are used as rough estimates ---------------------
; ---- This is the average M/Teff/L of the observed stellar population ----
M=1.25
Teff=6000
L=1.5
b0=M^(0.5 -a0) * (Teff/Teff_sun)^(3-3.5*a0) / L^(0.75-a0)
; -------------------------------------------------------------------------


min_dpx=4

Ndata=n_elements(freq)

if max(freq) lt 10 then freq=freq*1d6 ; convert from Hz to microHz

if fold_range[0] lt min(freq) then fold_range[0]=min(freq)
if fold_range[1] gt max(freq) then fold_range[1]=max(freq)

resol=freq[1]-freq[0]

x0=reform(freq)
s0=reform(spec_reg)

if n_elements(nu_rebin) ne 0 then begin
	if nu_rebin gt resol then begin
		Nrebin=(max(x)-min(x))/nu_rebin ; new number of points
		x=congrid(x0, Nrebin, /interp)
		s=congrid(s0, Nrebin, /interp)
		print, 'Rebin requested...'
		print, 'Previous resolution:', resol, ' for a total of ' + strtrim(Ndata,2) + ' points'
		print, 'New resolution:', nu_rebin, ' for a total of ' + strtrim(Nrebin,2) + ' points'
		resol=nu_rebin
	endif else begin
		x=x0
		s=s0
		print, 'rebin of the spectrum requested, but the user provided a rebin resolution'
		print, 'nu_rebin='+strtrim(nu_rebin,2) + ' smaller than the current spectrum '
		print, 'resolution resol='+strtrim(resol,2)+ ' !'
		print, '---->  nu_rebin parameter was discarded!'
	endelse	
endif else begin
x=x0
s=s0
endelse

print, 'Selecting data within the requested folding range ' + 'fmin='+strtrim(fold_range[0],2) + $
		' fmax='+strtrim(fold_range[1],2) + ' ...'
zone=where(x ge fold_range[0] and x le fold_range[1])
x=x[zone]
s=s[zone]

; ------- Setting the list of tested numax -----
Nstep=(fold_range[1] - fold_range[0])/numax_step
numax=numax_step*findgen(Nstep+1) + fold_range[0]
Dnu=Dnu_sun * b0 * ( numax/numax_sun )^a0 ; Stello et al. (2009) Defines the Heavy smooth coef 
coef=3. ;1.2 ; Defines the average window size (window=coef * Dnu)
Gamma0=double(aa) * K^po * M^(2.*po/3.) * numax^(-2.*po*(3.-4.*beta)/3.)  ; Benomar 2013. Defines the light smooth coef
coefsmoothG=5 ;3 ; Defines the coeficient of smoothing for the lightly smoothed spectra used for mode detection
; ----------------------------------------------

print, 'Creating a slice-averaged spectrum centered on each possible numax....'
print, '   1. Setting the number of datapoints/slice using the numax_step parameter...'
;pslice=where(x ge (numax[0] - total(Dnu[0]*coef)) AND x le  (numax[0] + total(Dnu[0]*coef)) )
pslice=where(x ge (numax[0] - total(Dnu[0]*coef)) AND x le  (numax[0] + total(Dnu[0]*coef)) )

px0_0=min( pslice ) & px1_0=max( pslice )
dpx=1.*(px1_0-px0_0)
print, '   Initial dpx set to ' + strtrim(dpx, 2)
if dpx lt min_dpx then begin
	print, ' This value is too small... please either increase numax_step or reduce nu_rebin (if used)'
	print, ' the minimum dpx is currently set to:', min_dpx
	print, 'program stopped'
	stop
endif

print, '   2. Computing...'
px0=px0_0 & px1=px1_0
	for cpt=long(0), n_elements(numax)-2 do begin
		;print, px0, px1
		x_tmp=x[px0:px1] 
		s_tmp=s[px0:px1] ; The median of this is used to define the noise background
	    
	    ; smooth the original spectrum over the expected mode width*coefsmoothG
	    ss2=smooth(s0, coefsmoothG*Gamma0[cpt]/resol, /edge_truncate) 
	    ss2=ss2[zone]
	    s2_tmp=ss2[px0:px1] ; we look for the max on this lightly smoothed spectra
	    ;stop
	    x_smooth=autofill(x_smooth, mean(x_tmp), cpt) ; update/initialize x_smooth
		s_smooth=autofill(s_smooth, median(s_tmp), cpt)
	    s_smooth_max=autofill(s_smooth_max, max(s2_tmp), cpt)
	    
	    if px1-px0 lt 10 then scoefmax=2. else scoefmax=1.*(px1-px0)/50 ; basically smooth over 0.1*Dnu
	    if px1-px0 lt 2 then scoefmax=1.
	    ;s_smooth_max_smooth=autofill(s_smooth_max_smooth, max(smooth(s_tmp, 2, /edge_truncate)), cpt) ; update/initialize s_cc
		s_smooth_stddev=autofill(s_smooth, stddev(s_tmp), cpt)
		;stop
		; --- shifting to the next slice ---
		pslice=where(x ge (numax[cpt+1] - total(Dnu[cpt+1]*coef)) AND x le  (numax[cpt+1] + total(Dnu[cpt+1]*coef)) )
		px0=min( pslice )
		px1=max( pslice )
		; ----------------------------------
		;stop  ; debug only
	endfor
print, '   3. Extracting the guess for the maximum position and its significance...'
	x_r=x_smooth
	;s_r=smooth(s_smooth_max_smooth/s_smooth, sfactor, /edge) ; this should mitigate the noise background increase at LF
	s_r=smooth(s_smooth_max/s_smooth, sfactor, /edge_truncate) ; this should mitigate the noise background increase at LF
	s_r0=s_r
	s_r=(s_r-min(s_r))
	maxsr=max(s_r)
	s_r2=s_r/maxsr
    s_r_stddev=s_r + smooth(s_smooth_stddev/s_smooth/maxsr, sfactor, /edge_truncate) 
    ; ------ Output Values -------
    numax_guess0=x_r[where(s_r eq max(s_r))]
    numax_guess0=numax_guess0[0] ; to solve issues when several s_r have same max
    uncertainty0=interpol(Dnu, numax, numax_guess0) * coef ; the uncertainty is the half width of the slice
    
    ; --- refining using a polyfit, order 2 ----
    rzone=where(x_r ge numax_guess0 - 2*uncertainty0 AND x_r le numax_guess0 + 2*uncertainty0)
    peak_coef=poly_fit(x_r[rzone], s_r[rzone], 3, yfit=yfit)
    sig_peak=max(yfit)
    xselect=x_r[rzone]
    numax_guess=xselect[where(yfit eq max(yfit))]
    numax_guess=numax_guess[0]
    uncertainty=interpol(Dnu, numax, numax_guess) * coef
	
    Gamma0_guess=double(aa) * K^po * M^(2.*po/3.) * numax_guess^(-2.*po*(3.-4.*beta)/3.) ; expected with a numax_guess
    s4Amax=smooth(s0, coefsmoothG*Gamma0_guess/resol, /edge_truncate) ; define the optimum smooth a around numax
    Amax_guess=max(s4Amax[where(x0 ge numax_guess[0] - uncertainty[0] AND x0 le numax_guess[0] + uncertainty[0])]) 
    ;Amax_guess=median(s[where(x0 ge numax_guess[0] - uncertainty[0] AND x0 le numax_guess[0] + uncertainty[0])]) ; we look at the spectrum with optimum smooth 
    coefnoise=5
    exit=0
    while coefnoise ge 1 AND exit eq 0 do begin
    	noise_pos=where(x_r le numax_guess - uncertainty*coefnoise OR x_r ge numax_guess + uncertainty*coefnoise)
    	if noise_pos[0] eq -1 then coefnoise=coefnoise/1.2
    	if noise_pos[0] ne -1 then begin
    		nel=n_elements(noise_pos)
    		if nel gt 10 then begin
    			noise_mean=mean(s_r0[noise_pos])
    			noise_stddev=stddev(s_r0[noise_pos]) ;*sqrt(n_elements(s_r0[noise_pos]))
    			exit=1 
    		endif else begin
    			print, 'Warning during calculation of the peak significance: Less than 10 points in order to calculate the noise level'
    			if nel lt 2 then begin
    				coefnoise=coefnoise/1.2
    				noise_pos=where(x_r ge numax_guess - uncertainty*coefnoise AND x_r le numax_guess + uncertainty*coefnoise)
    				print, 'Warning during calculation of the peak significance: Only few points to calculate noise level!'
    				print, '    !!!!!!   PROCEED WITH EXTREME CAUTION WHEN INTERPRETING THE PEAK SIGNIFICANCE   !!!!!!'
    				noise_mean=mean(s_r0[noise_pos])
    				noise_stddev=stddev(s_r0[noise_pos]) ;*sqrt(n_elements(s_r0[noise_pos]))
    			endif else begin
    				noise_mean=mean(s_r0[noise_pos])
    				noise_stddev=stddev(s_r0[noise_pos]) ;*sqrt(n_elements(s_r0[noise_pos]))
    			endelse
    			exit=1
    		endelse
    	endif
    endwhile
    
    ms_r=mean(s_r)
    ;sigs_r=stddev(s_r)
    ; ---- old procedure to evaluate the significance (BUGGY) ----
    sig_gauss_sigma=[1,  2   ,  3  ,  4   ,   5    , 6  ,  10,  20 ] ; significance in sigma
    sig_gauss_tab=[68.2, 95.4, 99.7, 99.99, 99.9999, 100,  100, 100] ; significance in percent
    ;sig_peak=max(s_r)/sigs_r ; in unit of sigma
    ;significance=interpol(sig_gauss_tab, sig_gauss_sigma, sig_peak) ; get the significance by interpolation
    ; -------------------------------------------------------------
    
    noise_cdf=gauss_cdf(noise_mean, noise_stddev)
    ;significance=interpol(noise_cdf[1,*], noise_cdf[0,*], sig_peak) ; get the significance by interpolation
    lvl_1sigma=interpol(noise_cdf[0,*], noise_cdf[1,*], 68.2) ; determine the 1 sigma detection level
    significance=interpol(sig_gauss_tab, sig_gauss_sigma, Amax_guess/(lvl_1sigma + min(s_r0)))
    sigs_r=lvl_1sigma
    print, '------------------------------------------------------------------------------'
    print, '       numax_guess              = '+ strtrim(numax_guess,2) + ' microHz'
    print, '       uncertainty (slice size) = '+ strtrim(uncertainty, 2) + ' microHz'
    print, '       Amax_guess               = ' + strtrim(Amax_guess, 2) + ' ppm^2/microHz'
    print, '       Significance probability = '+strtrim(significance,2) + ' %'
    print, '       1-sigma detection level  = '+strtrim(lvl_1sigma,2) + ' ppm^2/microHz'
    print, '------------------------------------------------------------------------------'
    ; ----------------------------
    s_r=(s_r-min(s_r))/max(s_r)
    yfit=yfit/max(yfit)
    ;stop
    if n_elements(psfile) ne 0 then begin
    	print, 'As per requested, saving the image with the first guess on' + strtrim(psfile, 2) + '.ps'
    	e=write_on_ps_on(psfile)
    	plot, x_r, s_r, thick=2, xtitle='Frequency (microHz)', ytitle='No unit', $
    		title='max(s)/average(s) calculated per slice of PSF. Peak significance: '+ $
    					strtrim(significance, 2) + '%', yr=[0,1], /yst
		;oplot, x_r, s_r_stddev, color=fsc_color('Red'), thick=2
		plots, [numax_guess, numax_guess], [0, 1], linestyle=2, color=fsc_color('Blue'), $
			thick=3
		xyouts, numax_guess + (fold_range[1]-fold_range[0])/100, 0.9, $
			textoidl('\nu_{max}='+strtrim(numax_guess, 2))+ '+/-' + strtrim(uncertainty, 2), $
			color=fsc_color('Blue'), charsize=1.75
		plots, [min(x_r), max(x_r)], [ms_r, ms_r], linestyle=3, thick=3
		plots, [min(x_r), max(x_r)], [sigs_r, sigs_r]/max(s_r0), $
			linestyle=2, thick=3, color=fsc_color('Red')
		oplot, x_r[rzone], yfit, color=fsc_color('Green'), thick=3

		;plots, [min(x_r), max(x_r)], [ms_r+sigs_r, ms_r+sigs_r], $
		;	linestyle=2, thick=3
		; --------- Second plot: original spectra in search range ---------
		xmin=fold_range[0]
		xmax=fold_range[1]
		ymin=0.
		ymax=Amax_guess*2	
		plot, x0, s4Amax, thick=2, xtitle='Frequency (microHz)', ytitle='ppm^2/microHz', $
    		title='PSF at optimal smooth for modes enhancement: Local view', /yst, /xst, $
    		xr=[xmin, xmax], yr=[ymin, ymax]
    	plots, [numax_guess, numax_guess], [ymin, ymax], linestyle=2, color=fsc_color('Blue'), $
			thick=3
		xyouts, numax_guess + (fold_range[1]-fold_range[0])/100, ymax- 10.*(ymax-ymin)/100., $
			textoidl('\nu_{max}='+strtrim(numax_guess, 2))+ '+/-' + strtrim(uncertainty, 2), $
			color=fsc_color('Blue'), charsize=1.75
		plots, [numax_guess-uncertainty, numax_guess+uncertainty], [Amax_guess, Amax_guess], $
			linestyle=2, color=fsc_color('Blue'), thick=3
		xyouts, numax_guess + (fold_range[1]-fold_range[0])/100, $
				Amax_guess + (ymax-ymin)/100., $
			textoidl('A_{max}='+strtrim(Amax_guess, 2)), $
			color=fsc_color('Blue'), charsize=1.75
		; --------- Third plot: The whole original spectra ---------
		xmin=0
		xmax=max(x0)
		ymin=0.
		ymax=Amax_guess*10
		plot, x0, s4Amax, thick=2, xtitle='Frequency (microHz)', ytitle='ppm^2/microHz', $
    		title='PSF at optimal smooth for modes enhancement: Global view', /yst, /xst, $
    		xr=[xmin, xmax], yr=[ymin, ymax]
    	plots, [numax_guess, numax_guess], [ymin, ymax], linestyle=2, color=fsc_color('Blue'), $
			thick=3
;		xyouts, numax_guess + (fold_range[1]-fold_range[0])/100, ymax- (ymax-ymin)/100., $
;			textoidl('\nu_{max}='+strtrim(numax_guess, 2))+ '+/-' + strtrim(uncertainty, 2), $
;			color=fsc_color('Blue'), charsize=1.75
		plots, [numax_guess-uncertainty, numax_guess+uncertainty], [Amax_guess, Amax_guess], $
			linestyle=2, color=fsc_color('Blue'), thick=3
;		xyouts, numax_guess + (fold_range[1]-fold_range[0])/100, $
;				Amax_guess + (ymax-ymin)/100., $
;			textoidl('A_{max}='+strtrim(Amax_guess, 2)), $
;			color=fsc_color('Blue'), charsize=1.75
		e=write_on_ps_off('')
	endif

;stop ; for debug
tab=[numax_guess, uncertainty, Amax_guess, significance]

return, tab ; significance is in percent.
end

; simple function that returns a curve of the cdf provided mu and sigma
; if x is not provided, we range it between -10 sigma and + 10 sigma
; this is an approximation... that might be suficient
function gauss_cdf, mu, sigma, x=x
	
	if n_elements(x) eq 0 then begin
		resol=20*sigma/199
		x=findgen(200)*resol  - 10*sigma + mu
	endif else begin
		resol=mean(deriv(x))
	endelse
	logcdf=-alog(sqrt(2. * !pi) * sigma) - 0.5*(x - mu)^2 / sigma^2
	
	maxi=total(exp(logcdf))
	cdf=dblarr(n_elements(x))
	xcdf=cdf
	for i=long(0), n_elements(x)-1 do begin
		cdf[i]=total(exp(logcdf[0:i]))/maxi
		xcdf[i]=x[i] + resol/2
	endfor
	out=dblarr(2, n_elements(x))
	out[0,*]=xcdf
	out[1,*]=cdf*100
return, out
end

; function that declare an array if ind=0, increment it otherwise
function autofill, arr, val, ind

	if ind eq 0 then new_arr=val else new_arr=[arr, val]

return, new_arr
end

; separate the directory from the filename
function splitfilename, namein
	
	b=byte(namein)
	posslash=where((b eq 47) OR (b eq 92))
	posdot=where(b eq 46)
	
	f=b[max(posslash)+1:max(posdot)-1]
	d=b[0:max(posslash)]
	ext=b[max(posdot)+1:*]
	
	return, [strtrim(d,2), strtrim(f,2), strtrim(ext,2)]
end

