pro nice_hist1D, hist, ps=ps, file_out=file_out, xr=xr, $
	title=title, xtitle=xtitle, overplot=overplot, col=col, fclose=fclose, show_stats=show_stats, $
	legendsize=legendsize, legend_precision=legend_precision, pmulti=pmulti

if n_elements(ps) eq 0 then ps=0
if n_elements(file_out) eq 0 then ps=0
if n_elements(title) eq 0 then title=''
if n_elements(xtitle) eq 0 then xtitle=''
if n_elements(overplot) eq 0 then overplot=0
if n_elements(pmulti) eq 0 then pmulti=0
if n_elements(col) eq 0 then col=['Black', 'Red', 'Red']
if n_elements(fclose) eq 0 then fclose=0 ; a variable that forces to close the ps file after the plot... to use in case of overplot
if n_elements(show_stats) eq 0 then show_stats=1 ; defaut, we plot the stats
if n_elements(legendsize) eq 0 then size_char=0.85 else size_char=legendsize
if n_elements(legend_precision) eq 0 then txtformat='(f8.2)' else txtformat='(f8.' + strtrim(floor(legend_precision),2) +  ')'
smooth_coef=1

		xscale=(max(hist[0,*]) - min(hist[0,*]))*0.1

		y_max=max(hist[1,*])*1.1
		if n_elements(xr) eq 0 then begin
			x_min=min(hist[0,*]) - xscale
			x_max=max(hist[0,*]) + xscale
		endif else begin
			x_min=xr[0]
			x_max=xr[1]
		endelse
		frac=0.02

		if pmulti eq 0 then size_c=1.5 else size_c=0.75
		
		if ps eq 1 and overplot eq 0 AND pmulti eq 0 then nimp,name=file_out,/paper,/eps

		if overplot eq 1 AND pmulti eq 0 then oplot, hist[0,*], smooth(hist[1,*], smooth_coef), psym=10, color=fsc_color(col[0])

		if overplot eq 0 OR pmulti eq 1 then plot, hist[0,*], smooth(hist[1,*], smooth_coef), psym=10,$
				xtitle=xtitle, ytitle='Probability density', title=title, charsize=size_c, $
				background=fsc_color('White'), color=fsc_color(col[0]),yr=[0, y_max], xr=[x_min, x_max], /xst


if show_stats eq 1 then begin
		tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

		stats=estimate_1_sigma_error( hist[0,*],hist[1,*],68.3,2,tab_critere)

		b_max=interpol(smooth(hist[1,*],smooth_coef), hist[0,*], stats[6])
		plots, [stats[6],stats[6]], [0, b_max] , color=fsc_color('Green'), thick=5. ; median
		xyouts, stats[6]+(stats[1]*0.01), b_max*(1-0.82), strtrim(string(stats[6], format=txtformat),1), color=fsc_color('Green'),charsize=size_char

		b_max=interpol(smooth(hist[1,*],smooth_coef), hist[0,*], stats[5])
		plots, [stats[5],stats[5]],[0, b_max],  color=fsc_color('Orange'), thick=5.,linestyle=2 ; -1sigma
		xyouts, stats[5]+(stats[1]*0.03	), y_max*(1-0.94), strtrim(string(stats[5], format=txtformat),1), color=fsc_color('Orange'),charsize=size_char

		b_max=interpol(smooth(hist[1,*],smooth_coef), hist[0,*], stats[7])
		plots,  [stats[7],stats[7]],[0, b_max], color=fsc_color('Orange'), thick=5.,linestyle=2 ; +1sigma
		xyouts, stats[7]+(stats[1]*0.03), y_max*(1-0.94), strtrim(string(stats[7], format=txtformat),1), color=fsc_color('Orange'),charsize=size_char
end

	if (ps eq 1 AND overplot eq 0 AND fclose eq 0) OR fclose eq 1 then fimp

end

function build_histogram, y,Nb_classes, normalize=normalize

if n_elements(normalize) eq 0 then normalize=0 ; by defaut, no normalisation

;Nb_classes=100.*N/10000.
N=n_elements(y)
hist=dblarr(N/Nb_classes)
x_hist=dblarr(N/Nb_classes)

hist=HISTOGRAM(y,nbins=1.*N/Nb_classes)

res=(max(y)-min(y))/(N/Nb_classes) ; the resolution of the graph
x_hist=findgen(N/Nb_classes)*res+min(y) + res/2 ; we res/2 because we center the bins... UPDATED THE 02/09/2013

output_table=dblarr(2,N/Nb_classes)
output_table[0,*]=x_hist
output_table[1,*]=hist

;if normalize eq 1 then output_table[1,*]=output_table[1,*]/N ; Not correct normalisation
if normalize eq 1 then begin
	int_y=int_tabulated(double(x_hist), double(hist), /double)
	output_table[1,*]=output_table[1,*]/int_y
endif

return, output_table
end


function build_histogram_2D, y1,y2,Nb_classes

;Nb_classes=100.*N/10000.
N1=n_elements(y1)
N2=n_elements(y2)
hist=dblarr(Nb_classes,Nb_classes)


;hist=HISTOGRAM(y,nbins=1.*N/Nb_classes)
;hist = HIST_2D( y1, y2, BIN1=Nb_classes , BIN2=Nb_classes )
hist = HIST_2D( y1, y2)
x1=findgen(N1/Nb_classes)*(max(y1)-min(y1))/(N1/Nb_classes)+min(y1)
x2=findgen(N2/Nb_classes)*(max(y2)-min(y2))/(N2/Nb_classes)+min(y2)

output_struc={x1:dblarr(N1/Nb_classes),x2:dblarr(N2/Nb_classes), hist:dblarr(N1/Nb_classes,N2/Nb_classes)}
output_struc.x1=x1
output_struc.x2=x2
output_struc.hist=hist
return, output_struc
end

; create an histogram of an y(x) binary variable (0 or 1)
; x : x variable
; y : y(x) variable
; Dx : step in x
function build_histo_2var,x,y,Dx,seuil
count=0l
x0=x[0]
Nt=n_elements(x) & Ntot=0l
while Ntot lt Nt do begin
	Nl=n_elements(where(x ge x0 and x lt x0+Dx)) ; number of elements in the range [x0,x0+xM]
	sum=0l
	for l=0,Nl-1 do begin
		i=count + l
		if y[i] ge seuil then sum=sum+1
	endfor
	sum=1d*sum/Nl ; number in percent by class
	if Ntot ne 0 then begin
		y_result=[y_result,sum]
		x_result=[x_result,x0+Dx/2]
	endif else begin
		y_result=sum
		x_result=x0+Dx/2
	endelse
	count=count+Nl
	x0=x0+Dx

	Ntot=Ntot+Nl ; Sum(Nl) calculated iteratively
endwhile

output=dblarr(2,n_elements(y_result))
output[0,*]=x_result
output[1,*]=y_result
return, output
end

; 2D cdf
function cdf_2D, pdf_2D

a=n_elements(pdf_2D[0,*]) & b=n_elements(pdf_2D[*,0])
cdf=dblarr(a,b)

N=a
for i=0, N-1 do for j=0, N-1 do cdf[i,j]=total(pdf_2D[0:i, 0:j])


return, cdf
end

pro show_correl_V3, dir_out, var1, var2, name_var1, name_var2, unit_var1, unit_var2,$
		 xline=xline, yline=yline, legend_precision=legend_precision

mvsini=0 ; Measure of the vsin i ... or not.

tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma
!p.thick=3

  if n_elements(yline) ne 0 then begin
  	if n_elements(yline) gt 1 then begin
    	ylines=dblarr(3, 4)
    	ylines[0,*]=[yline[0], yline[0], min(var1), max(var1)]
    	ylines[1,*]=[yline[1], yline[1], min(var1), max(var1)]
    	ylines[2,*]=[yline[2], yline[2], min(var1), max(var1)] 
  	endif else begin
		ylines=dblarr(3, 4)
    	ylines[0,*]=[yline[0], yline[0], min(var1), max(var1)]
  	endelse
  endif else begin
    Inc=-1
    errm_Inc=-1
    errp_Inc=-1
    ylines=-1
    range_B=0
  endelse 
  
   if n_elements(xline) ne 0 then begin
  	if xline ge 0 then begin
  		Nxline=n_elements(xline[0,*])
    	xlines=dblarr(Nxline, 4)
    	for kk=0, Nxline-1 do xlines[kk,*]=[min(var2), max(var2), xline[kk], xline[kk]]
  	endif else begin
		xlines=dblarr(3, 4)
    	xlines[0,*]=[min(var2), max(var2), xline[0], xline[0]]
  	endelse
  endif else begin
    xline=-1
    errm_xline=-1
    errp_xline=-1
    xlines=-1
    range_B=0
  endelse 
  
nimp,name=dir_out+ name_var1 + name_var2+'.eps',/paper,/eps ; Savita's code

	N1=n_elements(var1)
	N2=n_elements(var2)

	if N1 eq N2 then Nb_classes=100*N1/10000. else begin
		print, 'BIZARRERIE !!! Tailles des tables diff�rentes. ARRET PREMATRUEE'
		stop
	endelse
	resol1=( max(var1) - min(var1) ) / 20
	resol2=( max(var2) - min(var2) ) / 20

	hist = HIST_2D( var1, var2, bin1=resol1,bin2=resol2,min1=min(var1),min2=min(var2))/(1.*N1)

	x_s=max([n_elements(hist[*,0]), n_elements(hist[0,*])])

	hist=congrid(hist, x_s, x_s)

	x1=findgen(n_elements(hist[*,0]))*(max(var1)-min(var1))/(n_elements(hist[*,0])-1)+min(var1)
	x2=findgen(n_elements(hist[0,*]))*(max(var2)-min(var2))/(n_elements(hist[0,*])-1)+min(var2)

	k=20
	N2=k*n_elements(hist[0,*])
	pdf0=congrid(hist, N2, N2) ; VIRTUALY increase the precision of the map

	pdf=smooth(pdf0, 52,/edge_truncate)
	pdf=pdf
	pdf=pdf/total(pdf)
	pdf_show=smooth(pdf0, 30, /edge_truncate)
	pdf_show=pdf_show
	pdf_show=pdf_show/total(pdf_show)

	x_int=rebin(x1, k*n_elements(x1))
	y_int=rebin(x2, k*n_elements(x2))


	histo1D_1=HISTOGRAM(var1,nbins=0.5*N1/Nb_classes)
    x_hist1D_1=findgen(n_elements(histo1D_1))*(max(var1)-min(var1))/ $
    	(n_elements(histo1D_1)-1)+min(var1)
	err_1=estimate_1_sigma_error( x_hist1D_1,histo1D_1,68.3,2,tab_critere)

	h1D_1=dblarr(2, n_elements(histo1D_1))
	h1D_1[0,*]=x_hist1D_1 & h1D_1[1,*]=histo1D_1

    histo1D_2=HISTOGRAM(var2,nbins=1.*N1/Nb_classes)
    x_hist1D_2=findgen(n_elements(histo1D_2))*(max(var2)-min(var2))/ $
    	(n_elements(histo1D_2)-1)+min(var2)
	err_2=estimate_1_sigma_error( x_hist1D_2,histo1D_2,68.3,2,tab_critere)

	h1D_2=dblarr(2, n_elements(histo1D_2))
	h1D_2[0,*]=x_hist1D_2 & h1D_2[1,*]=histo1D_2

	loadct,27

  s0=min(x_int)
  s1=max(x_int)
  range_A=[s0, s1]
  if xlines[0] ne -1 then begin
    lmin=min(xlines[*,2:3]) ; minimum of the lines in the splitting axis
    lmax=max(xlines[*,2:3]) ; maximum of the lines in the splitting axis
  endif else begin
    lmin=1d6 ; dummy if no lines
    lmax=-1 ; dummy if no lines
  endelse
  range_A[0]=min([s0, lmin])
  range_A[1]=max([s1, lmax*1.1])
  
  
  i0=min(y_int)
  i1=max(y_int)
  range_B=[i0, i1]
  if ylines[0] ne -1 then begin
    lmin=min(ylines[*,2:3]) ; minimum of the lines in the splitting axis
    lmax=max(ylines[*,2:3]) ; maximum of the lines in the splitting axis
  endif else begin
    lmin=1d6 ; dummy if no lines
    lmax=-1 ; dummy if no lines
  endelse
  range_B[0]=min([i0, lmin])
  range_B[1]=max([i1, lmax*1.1])
  
	plot_fullhist_2D, x_int, y_int,transpose(pdf_show)/max(pdf_show), h1D_1,  h1D_2, err_1, err_2, $
		Name_var1 + ' ' + unit_var1, Name_var2 + ' ' + unit_var2, 1, $
		range_A=0, range_B=0, xlines=xlines, ylines=ylines, legend_precision=legend_precision;, projected=hist_vsini, label_proj=Name_proj, range_proj=0, stat_proj=err_vsini

fimp

end

pro plot_lines, line, xvar, yvar, type

  if n_elements(line) ne 0 then begin
  		Nxline=n_elements(line)
    	;xlines=dblarr(Nxline, 4)
    	for kk=0, Nxline-1 do begin
    		;xlines[kk,*]=[min(xvar), max(yvar), xline[kk], xline[kk]]
  			if type eq 'vertical' then $
  				plots, [line[kk], line[kk]], [min(yvar), max(yvar)], color=fsc_color('Dark Gray'), thick=5, linestyle=4
  			if type eq 'horizontal' then $
  				plots, [min(xvar), max(xvar)], [line[kk], line[kk]], color=fsc_color('Dark Gray'), thick=5, linestyle=4
  			if type ne 'horizontal' AND type ne 'vertical' then begin
  				print, 'Warning in plot_xlines <--- build_histogram_v2.pro'
  				print, '    requested lines not plotted because their type is incorrect'
  				print, '    Please set type to horizontal or vertical'
  			endif
  		endfor
  endif else begin
    xline=-1
    xlines=-1
  endelse 

end

function get_ylines, yline, var
	if n_elements(yline) ne 0 then begin
  	if n_elements(yline) ge 5 then begin
  	  	Nyline=n_elements(yline[0,*])
    	ylines=dblarr(Nyline, 4)
    	for kk=0, Nyline-1 do ylines[kk,*]=[yline[kk], yline[kk], min(var), max(var)]
  	endif else begin
		ylines=dblarr(3, 4)
    	ylines[0,*]=[yline[0], yline[0], min(var), max(var)]
  	endelse
  endif else begin
    yline=-1
    ylines=-1
  endelse 

end

; Show a1, a2, a3 and inclination in a matricial form ==> correlations
pro show_hist_a1_a2_a3_inc_matrix, a1, a2, a3, inc, dir_out, name_out, Radius=Radius, err_Radius=err_Radius, spectro_vsini=spectro_vsini, $
		a1_lines=a1_lines, a2_lines=a2_lines, $ ;, seismo_only=seismo_only
		a3_lines=a3_lines, inc_lines=inc_lines
	

  tab_critere=[2.25,16,50,84,97.75]
  N1=n_elements(a1)
  
  if n_elements(spectro_vsini) eq 0 then spectro_inc_rot=-1 else spectro_inc_rot=spectro_vsini
  
  resol1=( max(a1) - min(a1) ) / 30
  resol2=( max(a2) - min(a2) ) / 30
  hist = HIST_2D(a1, a2,bin1=resol1,bin2=resol2,min1=min(a1),min2=min(a2))/(1.*N1)
  x_s=n_elements(hist[*,0])
  hist=congrid(hist, x_s, x_s)
  x_a1_a1a2=findgen(n_elements(hist[*,0]))*(max(a1)-min(a1))/(n_elements(hist[*,0])-1)+min(a1)
  x_a2_a1a2=findgen(n_elements(hist[0,*]))*(max(a2)-min(a2))/(n_elements(hist[0,*])-1)+min(a2)
  hist_a1a2=hist
  
  resol1=( max(a1) - min(a1) ) / 30
  resol2=( max(a3) - min(a3) ) / 30
  hist = HIST_2D(a1, a3,bin1=resol1,bin2=resol2,min1=min(a1),min2=min(a3))/(1.*N1)

  hist=congrid(hist, x_s, x_s)
  x_a1_a1a3=findgen(n_elements(hist[*,0]))*(max(a1)-min(a1))/(n_elements(hist[*,0])-1)+min(a1)
  x_a3_a1a3=findgen(n_elements(hist[0,*]))*(max(a3)-min(a3))/(n_elements(hist[0,*])-1)+min(a3) 
  hist_a1a3=hist
  
  resol1=( max(a1) - min(a1) ) / 30
  resol2=( max(inc) - min(inc) ) / 30
  hist = HIST_2D(a1, inc,bin1=resol1,bin2=resol2,min1=min(a1),min2=min(inc))/(1.*N1)

  hist=congrid(hist, x_s, x_s)
  x_a1_a1inc=findgen(n_elements(hist[*,0]))*(max(a1)-min(a1))/(n_elements(hist[*,0])-1)+min(a1)
  x_inc_a1inc=findgen(n_elements(hist[0,*]))*(max(inc)-min(inc))/(n_elements(hist[0,*])-1)+min(inc)
  hist_a1inc=hist

  resol1=( max(a2) - min(a2) ) / 30
  resol2=( max(a3) - min(a3) ) / 30
  hist = HIST_2D(a2, a3,bin1=resol1,bin2=resol2,min1=min(a2),min2=min(a3))/(1.*N1)
  hist=congrid(hist, x_s, x_s)
  x_a2_a2a3=findgen(n_elements(hist[*,0]))*(max(a2)-min(a2))/(n_elements(hist[*,0])-1)+min(a2)
  x_a3_a2a3=findgen(n_elements(hist[0,*]))*(max(a3)-min(a3))/(n_elements(hist[0,*])-1)+min(a3)
  hist_a2a3=hist

  resol1=( max(a2) - min(a2) ) / 30
  resol2=( max(inc) - min(inc) ) / 30
  hist = HIST_2D(a2, inc,bin1=resol1,bin2=resol2,min1=min(a2),min2=min(inc))/(1.*N1)
  hist=congrid(hist, x_s, x_s)
  x_a2_a2inc=findgen(n_elements(hist[*,0]))*(max(a2)-min(a2))/(n_elements(hist[*,0])-1)+min(a2)
  x_inc_a2inc=findgen(n_elements(hist[0,*]))*(max(inc)-min(inc))/(n_elements(hist[0,*])-1)+min(inc)
  hist_a2inc=hist

  resol1=( max(a3) - min(a3) ) / 30
  resol2=( max(inc) - min(inc) ) / 30
  hist = HIST_2D(a3, inc,bin1=resol1,bin2=resol2,min1=min(a3),min2=min(inc))/(1.*N1)
  hist=congrid(hist, x_s, x_s)
  x_a3_a3inc=findgen(n_elements(hist[*,0]))*(max(a3)-min(a3))/(n_elements(hist[*,0])-1)+min(a3)
  x_inc_a3inc=findgen(n_elements(hist[0,*]))*(max(inc)-min(inc))/(n_elements(hist[0,*])-1)+min(inc)
  hist_a3inc=hist
  
  ;if n_elements(x_a1_a1a2) ne n_elements(x_a2_a1a2) AND n_elements(x_a1_a1a2) ne n_elements(x_a1_a1a3) AND $
  ;		n_elements(x_a1_a1a2) ne n_elements(x_a3_a1a3) AND n_elements(x_a1_a1a2) ne n_elements(x_a1_a1inc) AND $
  ;		n_elements(x_a1_a1a2) ne n_elements(x_inc_a1inc) AND n_elements(x_a1_a1a2) ne n_elements(x_a2_a2a3) AND $
  ;		n_elements(x_a1_a1a2) ne n_elements(x_a3_a2a3) AND n_elements(x_a1_a1a2) ne n_elements(x_a2_a2inc) AND $
  ;		n_elements(x_a1_a1a2) ne n_elements(x_inc_a2inc) AND n_elements(x_a1_a1a2) ne n_elements(x_a3_a3inc) AND $
  ;		n_elements(x_a1_a1a2) ne n_elements(x_inc_a3inc) $
  ;then begin
  ;	print, 'Size of the table are not matching! This procedure can only handle vectors of samples of the same size!'
  ;	print, 'The program will stop now'
  ;	stop		
  ;endif else begin
  	x_a1=x_a1_a1a2
  	x_a2=x_a2_a1a2
  	x_a3=x_a3_a1a3
  	x_inc=x_inc_a1inc
  	N0=n_elements(x_a1)
  ;endelse
  
  k=10
  N2=k*N0
  if n_elements(x_a1) ne n_elements(x_inc) then begin
  	print, 'debug stop'
  	stop
  endif
  struc_a1a2_show=rebin_2dpdf(hist_a1a2, x_a1, x_a2, N2)
  struc_a1a3_show=rebin_2dpdf(hist_a1a3, x_a1, x_a3, N2)
  struc_a1inc_show=rebin_2dpdf(hist_a1inc, x_a1, x_inc, N2)
  struc_a2a3_show=rebin_2dpdf(hist_a2a3, x_a2, x_a3, N2)
  struc_a2inc_show=rebin_2dpdf(hist_a2inc, x_a2, x_inc, N2)
  struc_a3inc_show=rebin_2dpdf(hist_a3inc, x_a3, x_inc, N2)
 
  hist_a1=build_histogram(a1, 4*N2, normalize=1)
  err_a1=estimate_1_sigma_error(hist_a1[0,*],hist_a1[1,*],68.3,2,tab_critere)

  hist_a2=build_histogram(a2, 4*N2, normalize=1)
  err_a2=estimate_1_sigma_error(hist_a2[0,*],hist_a2[1,*],68.3,2,tab_critere)
  
  hist_a3=build_histogram(a3, 4*N2, normalize=1)
  err_a3=estimate_1_sigma_error(hist_a3[0,*],hist_a3[1,*],68.3,2,tab_critere)

  hist_inc=build_histogram(inc, 4*N2, normalize=1)
  err_inc=estimate_1_sigma_error(hist_inc[0,*],hist_inc[1,*],68.3,2,tab_critere)

  s0=min(a1)
  s1=max(a1)
  range_A=[s0, s1]
  if spectro_inc_rot ne -1 then begin
  	smin=spectro_inc_rot[2, min(where(spectro_inc_rot[0,*] ge max(y_int)))] ; find the spliting when vsini curve cross the max of shown inclination
  endif else begin
  	smin=1d6 ; dummy value that will always be rejected
  endelse
  range_A[0]=min([s0, smin])
  range_A[1]=s1 ;max([s1])
  if range_A[0] lt 0 then range_A[0]=0
  
  i0=min(inc)
  i1=max(inc)
  range_D=[i0, i1]
  ;if ylines[0] ne -1 then begin
  ;  lmin=min(ylines[*,2:3]) ; minimum of the lines in the splitting axis
  ;  lmax=max(ylines[*,2:3]) ; maximum of the lines in the splitting axis
  ;endif else begin
  ;  lmin=1d6 ; dummy if no lines
  ;  lmax=-1 ; dummy if no lines
  ;endelse
  range_D[0]=i0 ;min([i0, lmin])
  range_D[1]=i1 ;max([i1, lmax*1.1])
  
  !p.thick=3
  nimp,name=dir_out+name_out+'PDF_a1a2a3inc.eps',/paper,/eps ; Savita's code
    Name_split=textoidl('a1=<\delta\nu_{n,l}> (\mu')+'Hz)'
    Name_a2='!3a2'
    Name_a3=textoidl('a3  (\mu')+'Hz)'
    Name_inc='!3Stellar inclination (degree)' ;+ textoidl('I_s (degree)')
   if n_elements(Radius) ne 0 AND n_elements(err_Radius) ne 0 then begin
		Name_proj = textoidl('v_{seis} sin i (km.s^{-1})')
   endif else begin
		Name_proj = textoidl('f_{seis} sin i (\mu')+'Hz)'
   endelse
;if seismo_only eq 1 then begin ; override the vsini and the Prot plots
;   lines=-1
;   spectro_inc_rot=-1
;endif
 
 loadct, 25 
 scoef=2

 ;print, 'debug point 1'
 ;stop
   ;                          A     B     C      D
	plot_fullhist_2D_MATRIX, struc_a1a2_show.x, struc_a2a3_show.x, struc_a3inc_show.x, struc_a1inc_show.y, $
		struc_a1a2_show.hist, struc_a1a3_show.hist, $
		struc_a1inc_show.hist, struc_a2a3_show.hist, struc_a2inc_show.hist, struc_a3inc_show.hist,$ 
		hist_a1,  hist_a2, hist_a3, hist_inc, $
		err_a1, err_a2, err_a3, err_inc, $
		Name_split, Name_a2, Name_a3, Name_inc, scoef, $
		range_A=range_A, range_D=range_D, $
		;range_B=range_B
		;range_C=range_C, range_D=range_D ;, $
		A_lines=a1_lines, B_lines=a2_lines, $ ;, seismo_only=seismo_only
		C_lines=a3_lines, D_lines=inc_lines
		;spectro_inc_rot=spectro_inc_rot, legend_precision=legend_precision
		
fimp

	
end


pro show_hist_a1_a2_a3_mag_asym_inc_matrix, a1, inc, a2=a2, a3=a3, mag_b=mag_b, mag_alfa=mag_alfa, asym=asym, extra=extra, dir_out, name_out
	

  print, extra

  if n_elements(spectro_vsini) eq 0 then spectro_inc_rot=-1 else spectro_inc_rot=spectro_vsini
  
  tab_critere=[2.25,16,50,84,97.75]
  N1=n_elements(a1)
  
  dim=2 ; mimimum we have a1 and inclination

  if n_elements(asym) ne 0 then if total(asym) ne 0 then dim=dim+1
  if n_elements(a2) ne 0 then if total(a2) ne 0 then  dim=dim+1
  if n_elements(a3) ne 0 then if total(a3) ne 0 then  dim=dim+1
  if n_elements(mag_b) ne 0 then if total(mag_b) ne 0 then  dim=dim+1
  if n_elements(mag_alfa) ne 0 then if total(mag_alfa) ne 0 then  dim=dim+1


  ; put all data into one single array
  data=dblarr(dim, N1)
  names=strarr(dim)  

  data[0,*]=a1
  data[1,*]=inc

  names[0]='a1'
  names[1]='inc'

  j=2
  if n_elements(asym) ne 0 then begin
	if total(asym) ne 0 then begin
		data[j,*]=asym
		names[j]='asym'
		j=j+1
	endif
  endif
  if n_elements(a2) ne 0 then begin
	if total(a2) ne 0 then begin
		data[j,*]=a2
		names[j]='a2'
		j=j+1
	endif
  endif
  if n_elements(a3) ne 0 then begin
	if total(a3) ne 0 then begin
		data[j,*]=1000.*a3
		names[j]='a3'
		j=j+1
	endif
  endif
  if n_elements(mag_b) ne 0 then begin
	if total(mag_b) ne 0 then begin
		data[j,*]=mag_b
		names[j]='sqrt(a1).cos(i)'
		j=j+1
	endif
  endif
  if n_elements(mag_alfa) ne 0 then begin
	if total(mag_alfa) ne 0 then begin
		data[j,*]=mag_alfa
		names[j]='sqrt(a1).sin(i)'
		j=j+1
	endif
  endif

  ;stop
  ; Compute the x-axises and the 2D maps...
  for i=0, n_elements(data[*,0])-1 do begin
  	for j=0, n_elements(data[*,0])-1 do begin
  		;if i le j then begin
  			resol1=( max(data[i,*]) - min(data[i,*]) ) / 30.
  			resol2=( max(data[j,*]) - min(data[j,*]) ) / 30.
  			hist = HIST_2D(data[i,*], data[j,*],bin1=resol1,bin2=resol2,min1=min(data[i,*]),min2=min(data[j,*]))/(1.*N1)
  			if i eq 0 AND j eq 0 then begin
  				 x_s=n_elements(hist[*,0])
  				 x_vals=dblarr(n_elements(data[*,0]), n_elements(data[*,0]), x_s)
  				 hists_raw=dblarr(n_elements(data[*,0]), n_elements(data[*,0]), x_s, x_s)
  			endif
  			hist=congrid(hist, x_s, x_s)
  			x_vals[i,j,*]=findgen(n_elements(hist[*,0]))*(max(data[i,*])-min(data[i,*]))/(n_elements(hist[*,0])-1)+min(data[i,*])
   			hists_raw[i,j,*,*]=hist
  		;endif
  	endfor
  endfor

  ; interpolate the maps....
  N0=n_elements(x_vals[0,0,*])
  k=10
  N2=k*N0


  xhists_2dshow=dblarr(n_elements(data[*,0]), n_elements(data[*,0]), N2)
  hists_2dshow=dblarr(n_elements(data[*,0]), n_elements(data[*,0]), N2, N2)
  errs=dblarr(n_elements(data[*,0]), 10)
  
  for i=0, n_elements(data[*,0])-1 do begin
	hist=build_histogram(data[i,*], 4.1*N2, normalize=1)
	err=estimate_1_sigma_error(hist[0,*],hist[1,*],68.3,2,tab_critere)

	if i eq 0 then hists=dblarr(n_elements(data[*,0]), 2, n_elements(hist[1,*]))
	hists[i,*,*]=hist
	errs[i, *]=err
  	for j=0, n_elements(data[*,0])-1 do begin
  	  ;if i le j then begin
	  	struc_show=rebin_2dpdf(reform(hists_raw[i,j,*,*]), reform(x_vals[i,j,*]), reform(x_vals[j,i,*]), N2)
		xhists_2dshow[i,j,*]=struc_show.x
		hists_2dshow[i,j, *,*]=struc_show.hist/max(struc_show.hist)
	  ;endif
	endfor
  endfor

	if n_elements(extra) ne 0 then begin
    if extra[0] ne 0 then begin
		  hist_extra=build_histogram(extra, 4.1*N2, normalize=1)
		  err_extra=estimate_1_sigma_error(hist[0,*],hist[1,*],68.3,2,tab_critere)
		  label_extra='a2'
	  endif else begin
      hist_extra=dblarr(2, 10)
      err_extra=dblarr(8)
      label_extra=''
    endelse
  endif else begin
		hist_extra=dblarr(2, 10)
		err_extra=dblarr(8)
		label_extra=''
	endelse
   ;stop
  !p.thick=3
  nimp,name=dir_out+name_out+'PDF_multidim.eps',/paper,/eps ; Savita's code

  ; Handling ranges and variables names
  ranges=dblarr(n_elements(data[*,0]), 2)
  labels=strarr(n_elements(data[*,0]))
  for i=0, n_elements(data[*,0])-1 do begin
  	if names[i] eq "a1" then begin
  		labels[i]=textoidl('a1 (\mu')+'Hz)' ; =<\delta\nu_{n,l}>
  		s0=min(data[i,*])
  		s1=max(data[i,*])
  		ranges[i,*]=[s0, s1]
  		if spectro_inc_rot ne -1 then begin
  			smin=spectro_inc_rot[2, min(where(spectro_inc_rot[0,*] ge max(y_int)))] ; find the spliting when vsini curve cross the max of shown inclination
  		endif else begin
  			smin=1d6 ; dummy value that will always be rejected
  		endelse
  		ranges[i,0]=min([s0, smin])
  		ranges[i,1]=s1 ;max([s1])
  		if ranges[i, 0] lt 0 then ranges[i, 0]=0
  	endif
  	if names[i] eq "inc" then begin
  		labels[i]='!3Stellar inc. (degree)'
  		i0=min(data[i,*])
  		i1=max(data[i,*])
  		ranges[i,*]=[i0, i1]
  	endif
  	if names[i] eq "a2" then begin
  		;labels[i]=textoidl('\beta_0 (no unit)')
		labels[i]=textoidl('10^5 \Delta R/R')
		label_extra=labels[i]
  	endif
  	if names[i] eq "dRovR" then begin
  		labels[i]=textoidl('10^5 \Delta R/R')
  	endif
  	if names[i] eq "a3" then begin
  		;labels[i]=textoidl('a3  (\mu')+'Hz)'
		labels[i]='a3  (nHz)'
  	endif
  	if names[i] eq "asym" then begin
  		labels[i]="Mode asym. (no unit)"
  	endif
  	if names[i] eq "sqrt(a1).cos(i)" then begin
  		labels[i]=textoidl('sqrt(a1).cos(i) (\mu')+'Hz)'
  	endif
  	if names[i] eq "sqrt(a1).sin(i)" then begin
  		labels[i]=textoidl('sqrt(a1).sin(i) (\mu')+'Hz)'
  	endif
  endfor
   
 loadct, 25 
 scoef=2

	posplots=dblarr(2, n_elements(data[*,0]))
	if n_elements(data[*,0]) ne 5 AND n_elements(data[*,0]) ne 4 then begin
		print, 'This function is designed to better work with 5 or 4 parameters which are the rotation/asymetry/inclination parameters'
		print, 'It was detected that it was not a 5 or 4 parameter space, thus no specific order of the parameter is used'
		posplots[0,*]=findgen(n_elements(data[*,0]))
		posplots[1,*]=findgen(n_elements(data[*,0]))
	endif
	if n_elements(data[*,0]) eq 5 then begin
		posplots[0,*]=[0,3,2,1,4]
		posplots[1,*]=[0,3,2,1,4]
	endif
	if n_elements(data[*,0]) eq 4 then begin
		posplots[0,*]=[0,3,2,1] ;[0,2,1,3]
		posplots[1,*]=[0,3,2,1];[0,2,1,3]
		;posplots[1,*]=[0,1,2,3]
	endif
  smoothcoef=2
	plot_fullhist_2D_MATRIX_Nd_v2, xhists_2dshow, hists_2dshow, hists, errs, labels, $
		ranges=ranges, smoothcoef, lines=lines, spectro_inc_rot=spectro_inc_rot, $
		legend_precision=legend_precision, posplots=posplots, hist_extra=hist_extra, err_extra=err_extra, label_extra=label_extra
fimp
	
end

function rebin_2dpdf, hist_in, x1, x2, N2

  pdf0=congrid(hist_in, N2, N2) ; VIRTUALY increase the precision of the map
  ;pdf_show=smooth(pdf0, 30, /edge_truncate)
  pdf_show=smooth(pdf0, 10, /edge_truncate)
  pdf_show=pdf_show/total(pdf_show)
  x=rebin(x1, N2)
  y=rebin(x2, N2)

  struc={x:dblarr(N2), y:dblarr(N2), hist:dblarr(N2,N2)}

  struc.x=x
  struc.y=y
  struc.hist=pdf_show
  
  return, struc
end

pro show_hist_a1_inc, splitting, angle, dir_out, name_out, Radius=Radius, err_Radius=err_Radius, spectro_vsini=spectro_vsini, $
		Prot=Prot, err_Prot=err_Prot, xline=xline, yline=yline, seismo_only=seismo_only
		    

tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

if n_elements(seismo_only) eq 0 then seismo_only=0

  !p.multi=0

  Teff_sun= 5777d ; same values as in the function: seismic_vsini
  Dnu_sun=135.1d
  numax_sun=3150d
  R_sun=6.96342d5 ; in km
  M_sun=1.98855d30 ; in kg

  if n_elements(spectro_vsini) ne 0 then begin
  	spectro_inc_rot=nu_rot_from_spec_vsini(Radius, spectro_vsini, err_spectro_vsini)
  	if err_spectro_vsini eq -1 then err_spectro_vsini=0.1*spectro_vsini
  endif else begin
  	spectro_inc_rot=-1
  endelse
  if n_elements(Prot) ne 0 then begin
  	if n_elements(err_Prot) ne 0 then begin
		if err_Prot eq -1 then err_Prot=0.1*Prot
		if err_Prot lt -1 then err_Prot=-1
	endif else begin
		err_Prot=-1
	endelse
  endif
  
  N1=n_elements(splitting)
  N2=n_elements(angle)
  if N1 eq N2 and n_elements(Nb_classes) eq 0 then begin
    Nb_classes=50*N1/10000
  endif
  if N1 ne N2 then begin
    print, 'BIZARRERIE !!! Tailles des tables differentes. ARRET PREMATUREE'
    stop
  endif
  
  if n_elements(Prot) ne 0 then begin
  	if Prot ge 0 then begin
    	Nsamples=N1
    	omega_Prot=gauss_random(1d6/Prot, 1d6*err_Prot/Prot^2, Nsamples)
    	test=where(omega_Prot lt 0) ; reject negative values that may arise from round off errors
    	if test[0] ne -1 then omega_Prot[test]= 1d-5 ; negatives are transformed into near 0 values
		stop
    	histProt=build_histogram(omega_Prot,Nb_classes, normalize=1)

    	s=estimate_1_sigma_error( histProt[0,*],histProt[1,*],68.3,2,tab_critere)
    	omega_Prot_stats=s[3:9] ;10^(s[3:9]) ; [0%, tab_critere, 100%]
    	omega_Prot=omega_Prot_stats[3] ;s[6]
    	errm_omega_Prot=omega_Prot_stats[3] - omega_Prot_stats[2];s[6] - s[5]
    	errp_omega_Prot=omega_Prot_stats[4] - omega_Prot_stats[3] ;s[7] - s[6]
    	lines=dblarr(3, 4)
    	lines[0,*]=[min(angle), max(angle), omega_Prot, omega_Prot]
    	lines[1,*]=[min(angle), max(angle), omega_Prot_stats[2], omega_Prot_stats[2]]
    	lines[2,*]=[min(angle), max(angle), omega_Prot_stats[4], omega_Prot_stats[4]] 
  	endif else begin
		lines=dblarr(3, 4)
    	lines[0,*]=[min(angle), max(angle), 1d6/Prot, 1d6/Prot]
    	lines[1,*]=-1 ; [min(angle), max(angle), omega_Prot_stats[2], omega_Prot_stats[2]]
    	lines[2,*]=-1 ; [min(angle), max(angle), omega_Prot_stats[4], omega_Prot_stats[4]] 
  	endelse
  endif else begin
    omega_Prot=-1
    omega_Prot_stats=-1
    omega_Prot=-1
    errm_omega_Prot=-1
    errp_omega_Prot=-1
    lines=-1
    range_A=0
  endelse 
  if n_elements(yline) ne 0 then begin
  	if n_elements(yline) gt 1 then begin
    	ylines=dblarr(3, 4)
    	ylines[0,*]=[yline[0], yline[0], min(splitting), max(splitting)]
    	ylines[1,*]=[yline[1], yline[1], min(splitting), max(splitting)]
    	ylines[2,*]=[yline[2], yline[2], min(splitting), max(splitting)] 
  	endif else begin
		ylines=dblarr(3, 4)
    	ylines[0,*]=[yline[0], yline[0], min(splitting), max(splitting)]
  	endelse
  endif else begin
    Inc=-1
    errm_Inc=-1
    errp_Inc=-1
    ylines=-1
    range_B=0
  endelse 
  
   if n_elements(xline) ne 0 then begin
  	if xline ge 0 then begin
  		Nxline=n_elements(xline[0,*])
    	xlines=dblarr(Nxline, 4)
    	for kk=0, Nxline-1 do xlines[kk,*]=[min(angle), max(angle), xline[kk], xline[kk]]
  	endif else begin
		xlines=dblarr(4)
    	xlines[0,*]=[min(angle), max(angle), xline[0], xline[0]]
  	endelse
  endif else begin
    xline=-1
    errm_xline=-1
    errp_xline=-1
    xlines=-1
    range_B=0
  endelse 

  if N1 gt 500000 then Nb_classes=1000. 
  if N1 lt 500000 AND N1 gt 200000 then Nb_classes=1000.
  if N1 lt 200000 AND N1 gt 50000  then Nb_classes=200.
  if N1 lt 50000 then Nb_classes=50.

  project=abs(splitting*sin(angle*!pi / 180d)) ; abs is there in case the prior authorize i<90 (cosI prior)

  if n_elements(Radius) ne 0 AND n_elements(err_Radius) ne 0 then begin ; If Radius then calculates vsini PDF and plot it
  	vsini=abs(vsini_pdf_RADIUS(project, Radius, err_Radius))
	hist_vsini=build_histogram(vsini, Nb_classes, normalize=1)
  	err_vsini=estimate_1_sigma_error( hist_vsini[0,*],hist_vsini[1,*],68.3,2,tab_critere)
  	;hist_vsini[1,*]=hist_vsini[1,*]/n_elements(vsini)
  endif else begin ; esle if no Radius then calculates nu.sini PDF and plot it
	hist_vsini=build_histogram(project, Nb_classes, normalize=1)
  	err_vsini=estimate_1_sigma_error( hist_vsini[0,*],hist_vsini[1,*],68.3,2,tab_critere)  
  endelse
  
  resol1=( max(splitting) - min(splitting) ) / 30
  resol2=( max(angle) - min(angle) ) / 30
  hist = HIST_2D( splitting, angle,bin1=resol1,bin2=resol2,min1=min(splitting),min2=min(angle))/(1.*N1)
  x_s=max([n_elements(hist[*,0]), n_elements(hist[0,*])])
  hist=congrid(hist, x_s, x_s)
  x1=findgen(n_elements(hist[*,0]))*(max(splitting)-min(splitting))/(n_elements(hist[*,0])-1)+min(splitting)
  x2=findgen(n_elements(hist[0,*]))*(max(angle)-min(angle))/(n_elements(hist[0,*])-1)+min(angle)

  k=20
  N2=k*n_elements(hist[0,*])
  pdf0=congrid(hist, N2, N2) ; VIRTUALY increase the precision of the map

  pdf=smooth(pdf0, 52,/edge_truncate)
  pdf=pdf
  pdf=pdf/total(pdf)
  pdf_show=smooth(pdf0, 20, /edge_truncate)
  pdf_show=pdf_show
  pdf_show=pdf_show/total(pdf_show)
  ;x_int=rebin(x1, k*n_elements(x1))
  ;y_int=rebin(x2, k*n_elements(x2))
  x_int=rebin(x1, N2)
  y_int=rebin(x2, N2)

  histo1D_1=HISTOGRAM(splitting,nbins=0.5*N1/Nb_classes)
  x_hist1D_1=findgen(n_elements(histo1D_1))*(max(splitting)-min(splitting))/ (n_elements(histo1D_1)-1) + min(splitting)
  err_1=estimate_1_sigma_error( x_hist1D_1,histo1D_1,68.3,2,tab_critere)

  histo1D_2=HISTOGRAM(angle,nbins=1.*N2/Nb_classes)
  x_hist1D_2=findgen(n_elements(histo1D_2))*(max(angle)-min(angle))/(n_elements(histo1D_2)-1) + min(angle)
  err_2=estimate_1_sigma_error( x_hist1D_2,histo1D_2,68.3,2,tab_critere)

  h1D_1=dblarr(2, n_elements(histo1D_1))
  h1D_1[0,*]=x_hist1D_1 & h1D_1[1,*]=histo1D_1
  h1D_2=dblarr(2, n_elements(histo1D_2))
  h1D_2[0,*]=x_hist1D_2 & h1D_2[1,*]=histo1D_2
 
  lines0=lines
  s0=min(x_int)
  s1=max(x_int)
  range_A=[s0, s1]
  if spectro_inc_rot ne -1 then begin
  	smin=spectro_inc_rot[2, min(where(spectro_inc_rot[0,*] ge max(y_int)))] ; find the spliting when vsini curve cross the max of shown inclination
  endif else begin
  	smin=1d6 ; dummy value that will always be rejected
  endelse
  if lines[0] ne -1 then begin
    lmin=min(lines[*,2:3]) ; minimum of the lines in the splitting axis
    lmax=max(lines[*,2:3]) ; maximum of the lines in the splitting axis
  endif else begin
    lmin=1d6 ; dummy if no lines
    lmax=-1 ; dummy if no lines
  endelse
  range_A[0]=min([s0, smin, lmin])
  range_A[1]=max([s1, lmax*1.1])
  if range_A[0] lt 0 then range_A[0]=0
  
  i0=min(y_int)
  i1=max(y_int)
  range_B=[i0, i1]
  if ylines[0] ne -1 then begin
    lmin=min(ylines[*,2:3]) ; minimum of the lines in the splitting axis
    lmax=max(ylines[*,2:3]) ; maximum of the lines in the splitting axis
  endif else begin
    lmin=1d6 ; dummy if no lines
    lmax=-1 ; dummy if no lines
  endelse
  range_B[0]=min([i0, lmin])
  range_B[1]=max([i1, lmax*1.1])
  
  !p.thick=3
  nimp,name=dir_out+name_out+'a1inc.eps',/paper,/eps ; Savita's code
    Name_split=textoidl('f_{seis}=<\delta\nu_{n,l}> (\mu')+'Hz)'
    Name_inc='!3Stellar inclination (degree)' ;+ textoidl('I_s (degree)')
   if n_elements(Radius) ne 0 AND n_elements(err_Radius) ne 0 then begin
		Name_proj = textoidl('v_{seis} sin i (km.s^{-1})')
   endif else begin
		Name_proj = textoidl('f_{seis} sin i (\mu')+'Hz)'
   endelse
if seismo_only eq 1 then begin ; override the vsini and the Prot plots
   lines=-1
   spectro_inc_rot=-1
endif
loadct,25
    plot_fullhist_2D, x_int, y_int,transpose(pdf_show)/max(pdf_show), h1D_1,  h1D_2, err_1, err_2, Name_split, Name_inc, 2, $
      range_A=range_A, range_B=range_B, lines=lines, xlines=xlines, ylines=ylines, projected=hist_vsini, label_proj=Name_proj, range_proj=0, stat_proj=err_vsini, $
      spectro_inc_rot=spectro_inc_rot
  fimp
end

; usefull to add the median/1sig/2sig intervals in a 1D pdf
pro replot_hist_1D_wlegends
tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

;; ***** HD16392-Rafa *****
;	;dir='C:\Work_dir\Analysis_Results\HD19392-Rafa_3_spline\files\' ; involving spline for l=0 & l=2
;	dir='C:\Work_dir\Analysis_Results\HD19392-Rafa_3_spline\RELAXED_priorVl2\4.15M_SAMPLES\files\' ; involving spline for l=0 & l=2
;	output_dir=dir
;	file1='individual_param_11.sav' ; V_l1
;	file2='individual_param_12.sav'	; V_l2
;
;; ************************

;; ***** HD16392-Rafa *****
;	dir='C:\Work_dir\Analysis_Results\HD19392-Rafa_l3free\files\Rafa_l3free_6M\'
;	output_dir=dir
;	file1='individual_param_11.sav' ; V_l1
;	file2='individual_param_12.sav'	; V_l2
;
;; ************************

;; ********* HAT P7 ********
;dir='C:\Work_dir\Analysis_Results\Hat-P7\NEW_ANALYSISES\Hat-P7-WeakP\files\'
;;dir='C:\Work_dir\Analysis_Results\Hat-P7\NEW_ANALYSISES\Hat-P7-WeakP2\files\'
;;	file1='individual_param_64.sav' ; SPLITTING
;;	file2='individual_param_93.sav' ; INCLINATION
;; *************************

;; **** Tim's Twin Star *****
;	dir='C:\Work_dir\Analysis_Results\PE15\TWINS_LUKELEIA\10124866A\Files\'
;	output_dir=dir
;	file1='individual_param_14.sav' ; V_l1
;	file2='individual_param_15.sav'	; V_l2
;; ***************************

;; ***** PAT STAR *****
	;dir='C:\Work_dir\Analysis_Results\Pat_Star\New_SET_smooth\1.4M\files\' ; My spectrum
	dir='C:\Work_dir\SYDNEY_WORK\Analysis_Results\Pat_Star\Oct2012_smooth_PatSpec\files\' ; Pat spectrum
	output_dir=dir
	file1='individual_param_13.sav' ; V_l1
	file2='individual_param_14.sav'	; V_l2
	;file2='individual_param_15.sav'	; V_l3
;; ************************


	variable_name=''
	restore, dir+file1
	Vl1=param & name_Vl1=variable_name & N2=n_elements(Vl1)
	restore, dir+file2
	Vl2=param ;/sin(angle*!pi/180d)
	name_Vl2=variable_name & N1=n_elements(Vl2)

	resol1=( max(Vl1) - min(Vl1) ) / 30
	resol2=( max(Vl2) - min(Vl2) ) / 30

	hist_1=build_histogram(Vl1, 35000)
	err_1=estimate_1_sigma_error( hist_1[0,*],hist_1[1,*],68.3,2,tab_critere)

	hist_2=build_histogram(Vl2, 35000)
	err_2=estimate_1_sigma_error( hist_2[0,*],hist_2[1,*],68.3,2,tab_critere)

	print, ' ****** Vl1 *******'
	print, 'Median : '+ strtrim(err_1[6],1)
	print, 'Median (par function) : ' + strtrim(median(Vl1),1)
	print, '-1 sigma : '+ strtrim(err_1[6]-err_1[5],1)
	print, '+1 sigma : '+ strtrim(err_1[7]-err_1[6],1)
	print, '-2 sigma : '+ strtrim(err_1[6]-err_1[4],1)
	print, '+2 sigma : '+ strtrim(err_1[8]-err_1[6],1)
	print, ' ************'
		print, ' ****** Vl2 *******'
	print, 'Median : '+ strtrim(err_2[6],1)
	print, 'Median (par function) : ' + strtrim(median(Vl2),1)
	print, '-1 sigma : '+ strtrim(err_2[6]-err_2[5],1)
	print, '+1 sigma : '+ strtrim(err_2[7]-err_2[6],1)
	print, '-2 sigma : '+ strtrim(err_2[6]-err_2[4],1)
	print, '+2 sigma : '+ strtrim(err_2[8]-err_2[6],1)
	print, ' ************'

	smooth_coef=2
	size_char=1.5
	y_max=max(hist_1[1,*]/n_elements(Vl1))*1.1
	;e=write_on_ps_on(output_dir+'PDF_Vl1')
	nimp,name=output_dir+'PDF_Vl1.eps',/paper,/eps ; Savita's code
		plot, hist_1[0,*], smooth(hist_1[1,*]/n_elements(Vl1),smooth_coef), psym=10,charsize=2,$
			xtitle=textoidl('V^2_{l=1}'), ytitle='Probability density',$
			;xtitle=textoidl('V^2_{l=1}'), ytitle='Probability density',$
			background=fsc_color('White'), color=fsc_color('Black'),yr=[0, y_max],/yst

		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[6])/n_elements(Vl1)
		plots, [err_1[6], err_1[6]], [0, bmax],color=fsc_color('Green'), thick=8
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[7])/n_elements(Vl1)
		plots, [err_1[7], err_1[7]], [0, bmax],color=fsc_color('Orange'), linestyle=2, thick=6
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[5])/n_elements(Vl1)
		plots, [err_1[5], err_1[5]], [0, bmax],color=fsc_color('Orange'), linestyle=2, thick=6

		xyouts, err_1[6]+(err_1[1]*0.01), bmax*(1-0.82), strtrim(string(err_1[6], format='(f8.2)'),1), color=fsc_color('Green'),charsize=size_char

		xyouts, err_1[5]+(err_1[1]*0.03	), y_max*(1-0.94), strtrim(string(err_1[5], format='(f8.2)'),1), color=fsc_color('Orange'),charsize=size_char

		xyouts, err_1[7]+(err_1[1]*0.03), y_max*(1-0.94), strtrim(string(err_1[7], format='(f8.2)'),1), color=fsc_color('Orange'),charsize=size_char
	fimp


	;e=write_on_ps_off('')
	y_max=max(hist_2[1,*]/n_elements(Vl2))*1.2
	;e=write_on_ps_on(output_dir+'PDF_Vl2')
	nimp,name=output_dir+'PDF_Vl2.eps',/paper,/eps ; Savita's code
		plot, hist_2[0,*], smooth(hist_2[1,*]/n_elements(Vl2),smooth_coef), psym=10,charsize=2,$
			xtitle=textoidl('V^2_{l=2}'), ytitle='Probability density',$
			background=fsc_color('White'), color=fsc_color('Black'),yr=[0, y_max],/yst
		bmax=interpol(smooth(hist_2[1,*],smooth_coef), hist_2[0,*], err_2[6])/n_elements(Vl2)
		plots, [err_2[6], err_2[6]], [0, bmax],color=fsc_color('Green'), thick=8

		bmax=interpol(smooth(hist_2[1,*],smooth_coef), hist_2[0,*], err_2[7])/n_elements(Vl2)
		plots, [err_2[7], err_2[7]], [0, bmax],color=fsc_color('Orange'), linestyle=2, thick=6
		bmax=interpol(smooth(hist_2[1,*],smooth_coef), hist_2[0,*], err_2[5])/n_elements(Vl2)
		plots, [err_2[5], err_2[5]], [0, bmax],color=fsc_color('Orange'), linestyle=2, thick=6

		xyouts, err_2[6]+(err_2[1]*0.01), bmax*(1-0.82), strtrim(string(err_2[6], format='(f8.2)'),1), color=fsc_color('Green'),charsize=size_char

		xyouts, err_2[5]+(err_2[1]*0.03	), y_max*(1-0.94), strtrim(string(err_2[5], format='(f8.2)'),1), color=fsc_color('Orange'),charsize=size_char

		xyouts, err_2[7]+(err_2[1]*0.03), y_max*(1-0.94), strtrim(string(err_2[7], format='(f8.2)'),1), color=fsc_color('Orange'),charsize=size_char

	fimp
	;e=write_on_ps_off('')

print, 'End'
end


; usefull to add the median/1sig/2sig intervals in a 1D pdf
pro replot_hist_1D_wlegends_v2, dir_files, file_synthese
tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma

smooth_coef=4

restore, file_synthese ; to get parameters_length


Nhyper=parameters_length[8]/4
Nmax=parameters_length[0]

dir=dir_files

	output_dir=dir

	files=file_search(dir+'*.sav')

; ************************
e=write_on_ps_on('a')
Height_unit=textoidl('ppm^2/ \mu')+'Hz'
Freq_unit=textoidl('\mu')+'Hz'
angle=textoidl('\circ')
e=write_on_ps_off('')

cpt=0d
for i=0, n_elements(files)-1 do begin
	variable_name=''
	restore, files[i]
	N2=n_elements(param)

	resol1=( max(param) - min(param) ) / 30

	Nb_c=35000
	Nb_c=12000
	hist_1=build_histogram(param, Nb_c)
	err_1=estimate_1_sigma_error( hist_1[0,*],hist_1[1,*],68.3,2,tab_critere)

	if i eq total(parameters_length[0:5]) then variable_name[0,*]=textoidl('Inclination (degree)') ; sometime buggy... ensure that the symbol is ok.
 ;***** HYPER PARAMATERS LABELS TUNES MANUALY !!!! ****
	;if variable_name[0,*] eq 'Label a definir...'  OR variable_name[0,*] eq '' then begin
	if i ge total(parameters_length[0:7]) then begin
		if cpt lt Nhyper then begin
			variable_name[0,*]='Height  mixed mode '+', l=1 ('+freq_unit+')
			;stop
		endif
		if cpt ge Nhyper AND cpt lt 2*Nhyper then variable_name[0,*]='Frequency  mixed mode  '+', l=1 ('+freq_unit+')
		if cpt ge 2*Nhyper AND cpt lt 3*Nhyper then variable_name[0,*]='Width  mixed mode  '+', l=1 ('+freq_unit+')
		if cpt ge 3*Nhyper AND cpt lt 4*Nhyper then variable_name[0,*]='Splitting mixed mode  '+', l=1 ('+freq_unit+')
		if cpt ge 4*Nhyper AND cpt lt Nmax + 4*Nhyper then  variable_name[0,*]='Width l=2  '+'('+freq_unit+')
		if cpt ge Nmax + 4*Nhyper  AND cpt lt 2*Nmax + 4*Nhyper then  variable_name[0,*]='Splitting l=2  '+'('+freq_unit+')
		if cpt ge 2*Nmax + 4*Nhyper then  variable_name[0,*]='Splitting l=3  '+'('+freq_unit+')
		cpt=cpt+1
	endif

 ;******************************************************

 	if i lt 10 then begin
 		name_out=output_dir+'PDF_param_00'+strtrim(i,1)+'.eps'
 		output_ascii=output_dir+'Val_param_00'+strtrim(i,1)+'.ASCII'
 	endif
	if i ge 10 AND i lt 100 then begin
		name_out=output_dir+'PDF_param_0'+strtrim(i,1)+'.eps'
		output_ascii=output_dir+'Val_param_0'+strtrim(i,1)+'.ASCII'
	endif
	if i ge 100 then begin
		name_out=output_dir+'PDF_param_'+strtrim(i,1)+'.eps'
		output_ascii=output_dir+'Val_param_'+strtrim(i,1)+'.ASCII'
	endif

	openw,3, output_ascii
	printf,3, ' ****** param *******'
	;printf,3, 'Median : '+ strtrim(err_1[6],1)
	printf,3, 'Median : ' + strtrim(median(param),1)
	printf,3, '-1 sigma : '+ strtrim(err_1[6]-err_1[5],1)
	printf,3, '+1 sigma : '+ strtrim(err_1[7]-err_1[6],1)
	printf,3, '-2 sigma : '+ strtrim(err_1[6]-err_1[4],1)
	printf,3, '+2 sigma : '+ strtrim(err_1[8]-err_1[6],1)
	printf,3, ' ************'
	close,3

	y_max=max(hist_1[1,*]/N2)*1.1
	nimp,name=name_out,/paper,/eps ; Savita's code

		xtitle_v=variable_name[0,*]

		plot, hist_1[0,*], smooth(hist_1[1,*],smooth_coef)/N2, psym=10,charsize=1.5,$
			xtitle=xtitle_v, ytitle='Probability density',$
			;xtitle=textoidl('V^2_{l=1}'), ytitle='Probability density',$
			background=fsc_color('White'), color=fsc_color('Black'),yr=[0, y_max], thick=5
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[6])/N2
		plots, [err_1[6], err_1[6]], [0, bmax],color=fsc_color('Dark Gray'), linestyle=2, thick=7
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[7])/N2
		plots, [err_1[7], err_1[7]], [0, bmax],color=fsc_color('Red'), linestyle=3, thick=7
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[5])/N2
		plots, [err_1[5], err_1[5]], [0, bmax],color=fsc_color('Red'), linestyle=3, thick=7
	fimp
endfor

print, 'End'
end


; err format: [0%, 2.25%, 16%, 50%, 66.7%, 97.75%, 100%, max]
pro simple_hist_plot, x_hist, y_hist, err, file_out, smooth_coef, x_axis_range=x_axis_range

	if n_elements(x_axis_range) eq 0 then xr=[min(x_hist), max(x_hist)] else xr=x_axis_range

	size_char=1.5

	nimp,name=file_out,/paper,/eps ; Savita's code

	y_max=max(y_hist)*1.1

	xlabel=textoidl('\nu_{max} (\mu') + 'Hz)'

		plot, x_hist, smooth(y_hist ,smooth_coef), psym=10,charsize=2,$
			xtitle=xlabel, ytitle='Probability density',$
			background=fsc_color('White'), color=fsc_color('Black'),yr=[0, y_max],/yst, xr=xr, /xst
		bmax=interpol(smooth(y_hist ,smooth_coef), x_hist, err[3])
		plots, [err[3], err[3]], [0, bmax],color=fsc_color('Green'), thick=8, linestyle=2

		bmax=interpol(smooth(y_hist,smooth_coef), x_hist, err[4])
		plots, [err[4], err[4]], [0, bmax],color=fsc_color('Orange'), linestyle=2, thick=6
		bmax=interpol(smooth(y_hist,smooth_coef), x_hist, err[2])
		plots, [err[2], err[2]], [0, bmax],color=fsc_color('Orange'), linestyle=2, thick=6

deltaerr=xr[1] - xr[0]

		xyouts, err[3]-(deltaerr*0.11), bmax*(1-0.75), strtrim(string(err[3], format='(f8.1)'),1), color=fsc_color('Green'),charsize=size_char

		xyouts, err[2]-(deltaerr*0.11), y_max*(1-0.97), strtrim(string(err[2], format='(f8.1)'),1), color=fsc_color('Orange'),charsize=size_char

		xyouts, err[4]+(deltaerr*0.01), y_max*(1-0.97), strtrim(string(err[4], format='(f8.1)'),1), color=fsc_color('Orange'),charsize=size_char

	fimp

end

pro avg_splitting, dir, i0, imax, txt_ext
tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma


Nb_c=35000
Nb_c=12000

smooth_coef=4

;dir='C:\Work_dir\Analysis_Results\8561221_nol1_splitting\Files\'
output_dir=dir

param_t=0d
param_o=0d
cpt=0d

;i0=113
;imax=119
;txt_ext=' / Averaged 5 higher'

;i0=114
;imax=118
for i=i0, imax do begin

	variable_name=''
 	if i lt 10 then begin
 		file=file_search(dir+'param_00'+strtrim(imax,1) + '.sav')
 		;name_out=output_dir+'PDF_param_00'+strtrim(i,1)+'.eps'
 		;output_ascii=output_dir+'Val_param_00'+strtrim(i,1)+'.ASCII'
 	endif
	if i ge 10 AND i lt 100 then begin
		file=file_search(dir+'param_0'+strtrim(imax,1) + '.sav')
		;name_out=output_dir+'PDF_param_0'+strtrim(i,1)+'.eps'
		;output_ascii=output_dir+'Val_param_0'+strtrim(i,1)+'.ASCII'
	endif
	if i ge 100 then begin
		file=file_search(dir+'param_'+strtrim(imax,1) + '.sav')
		name_out=output_dir+'PDF_l1_split_avgALL'+'.eps'
		output_ascii=output_dir+'Val_param_l1_split_avgALL'+'.ASCII'
	endif

	restore, file
	N2=n_elements(param)

	resol1=( max(param) - min(param) ) / 30

	param_t=param_t + param/ ( 1d*imax - 1d*i0 )
	param_o=param_o + param

	print, i
endfor
	hist_1=build_histogram(param_t, Nb_c)
	err_1=estimate_1_sigma_error( hist_1[0,*],hist_1[1,*],68.3,2,tab_critere)

	openw,3, output_ascii
		printf,3, ' ****** param *******'
		;printf,3, 'Median : '+ strtrim(err_1[6],1)
		printf,3, 'Median : ' + strtrim(median(param),1)
		printf,3, '-1 sigma : '+ strtrim(err_1[6]-err_1[5],1)
		printf,3, '+1 sigma : '+ strtrim(err_1[7]-err_1[6],1)
		printf,3, '-2 sigma : '+ strtrim(err_1[6]-err_1[4],1)
		printf,3, '+2 sigma : '+ strtrim(err_1[8]-err_1[6],1)
		printf,3, ' ************'
	close,3

	y_max=max(hist_1[1,*]/N2)*1.1
	nimp,name=name_out,/paper,/eps ; Savita's code

	xtitle_v='l=1 rotational splitting ('+textoidl('\mu')+'Hz)' + txt_ext

		plot, hist_1[0,*], smooth(hist_1[1,*],smooth_coef)/N2, psym=10,charsize=1.5,$
			xtitle=xtitle_v, ytitle='Probability density',$
			;xtitle=textoidl('V^2_{l=1}'), ytitle='Probability density',$
			background=fsc_color('White'), color=fsc_color('Black'),yr=[0, y_max], thick=5
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[6])/N2
		plots, [err_1[6], err_1[6]], [0, bmax],color=fsc_color('Dark Gray'), linestyle=2, thick=7
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[7])/N2
		plots, [err_1[7], err_1[7]], [0, bmax],color=fsc_color('Red'), linestyle=3, thick=7
		bmax=interpol(smooth(hist_1[1,*],smooth_coef), hist_1[0,*], err_1[5])/N2
		plots, [err_1[5], err_1[5]], [0, bmax],color=fsc_color('Red'), linestyle=3, thick=7
	fimp
	print, 'End'
end


function cdf_sigmas, cdf, x, y, criteria
N=n_elements(x) ; CAREFULL WE ASSUME THAT THE CDF TABLE IS A SQUARE!

k=1
N2=k*n_elements(x)

;c=congrid(cdf, N2,N2, /interp) ; congrid allows to rebin to any dimension... rebin cannot do that!
c=rebin(cdf, N2, N2) ; VIRTUALY increase the precision of the map

;c2=smooth(c, 4) ; smooth effect
;c=c2
;c=cdf

x_int=rebin(x, k*n_elements(x), k*n_elements(x))
y_int=rebin(y, k*n_elements(y), k*n_elements(y))

pass=0 & pass_old=0
for i=0, N2-1 do begin
	pass_old=pass
	j=min(where(c[i,*] gt criteria))
	if j ne -1 then begin
		x_i=x_int[i] & y_j=y_int[j]
		if pass eq 0 then begin
			x_pos=x_i
			y_pos=y_j
		endif else begin
			x_pos=[x_pos, x_i]
			y_pos=[y_pos, y_j]
		endelse
		pass=pass+1 ; variable to check if we already pass the condition once Or not
		if pass ne 1 then if pass-pass_old ne 1 then stop ; we stop the process if we passed already once but still didn't pass iteratively then
	endif
endfor

positions=dblarr(2,n_elements(x_pos))
positions[0,*]=x_pos
positions[1,*]=y_pos

return, positions
end

function gaussian_2D, m1,m2,sig1,sig2, x, y

norm=1d/( 2*!pi*sig1*sig2 )

gauss=dblarr(n_elements(x), n_elements(y))
for i=0, n_elements(x)-1 do for j=0, n_elements(y)-1 do begin
	gauss[i,j]=exp(-0.5d*( (x[i]-m1)^2/sig1^2 + (y[j]-m2)^2/sig2^2) )
endfor
g=exp(-((x-m1)^2/sig1^2 + (y-m2)^2/sig2^2)/2d)
no=total(gauss)
result=gauss/ no

return, result
end

pro test_f

x=findgen(300)/30
y=findgen(300)/30

m1=5d & m2=5d
sig1=1d & sig2=0.5d

pdf=gaussian_2D(m1,m2,sig1,sig2, x, y)
cdf=cdf_2D(pdf)
pos=pdf_2D_sigmas(pdf, x,y, 0.667, 500) ; the pdf must be normalized
sig1x=[x[min(pos[0,*])] , x[max(pos[0,*])]] ; 1 sigma positions in x
sig1y=[y[min(pos[1,*])] , y[max(pos[1,*])]] ; 1 sigma positions in y

med=cdf_sigmas( cdf, x, y, 0.5)
sig_m=cdf_sigmas( cdf, y, x, 0.16)
sig_p=cdf_sigmas( cdf, y, x, 0.84)

contour, pdf,x,y,/fill,yr=[min(y),max(y)],/xst,/yst,$
	xtitle='x',ytitle='y',nlevel=20, $
	background=fsc_color('White'),color=fsc_color('Black') ,/nodata
contour, pdf,x,y,/fill,yr=[min(y),max(y)],/xst,/yst,$
	xtitle='x',ytitle='y', nlevel=20,$
	background=fsc_color('Blue'),color=fsc_color('Black'),/overplot
oplot, x[pos[0,*]], y[pos[1,*]], psym=3, color=fsc_color('red'),thick=4


oplot, med[0,*], med[1,*],color=fsc_color('blue');*2
oplot, sig_m[0,*], sig_m[1,*];*2
oplot, sig_p[0,*], sig_p[1,*];*2
print, 'End'
stop
end

function pdf_2D_sigmas, pdf, x,y, criteria, N_slices

N=n_elements(x) ; we assume that the df is a square!

Delta=max(pdf) ; we may instead use the 2D median coordinates
med=median(pdf, dimension=2) ; just to check if it is like this that we can find the median

inc=Delta/N_slices ; define the width of each slices... related to the precision

i=0 & int=0
while (i le N-1) AND (int lt criteria) do begin
	condi=where(pdf ge Delta - i*inc)
	int=total(pdf[condi])
	i=i+1
endwhile
print, 'Error on the integration: ', abs(int-criteria), '   Error In %: ', 100*abs(int/criteria - 1)

pos=array_indices(pdf, condi)
sigx=[x[min(pos[0,*])] , x[max(pos[0,*])]]
sigy=[y[min(pos[1,*])] , y[max(pos[1,*])]]

return, pos
end


pro compute_average_height
tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma


dir_l0='C:\Work_dir\Analysis_Results\HD19392-Rafa_l3free\files\Rafa_l3free_6M\l0\'
dir_l3='C:\Work_dir\Analysis_Results\HD19392-Rafa_l3free\files\Rafa_l3free_6M\l3\'
output_dir=dir_l3
f_core='individual_param_'

files_l0=file_search(dir_l0+f_core+'*.sav')
files_l3=file_search(dir_l3+f_core+'*.sav')

restore, files_l0[0] ; retrieve the size

Vl3=dblarr(n_elements(files_l0), n_elements(param)) ; define the table that will contain all the ratios, order by order.
Vl3_m=0d ; define the initial value for the sum
for i=0, n_elements(files_l0)-1 do begin
	restore, files_l0[i]
	Hl0=param
	restore, files_l3[i]
	Hl3=param

	Vl3[i,*]=Hl3/Hl0 ; compute the ratio, order by order
	Vl3_m=Vl3_m + Vl3[i,*]/n_elements(files_l0) ; compute the mean
endfor

	resol1=( max(Vl3_m) - min(Vl3_m) ) / 30

	hist_1=build_histogram(Vl3_m, 25)
	err_1=estimate_1_sigma_error( hist_1[0,*],hist_1[1,*],68.3,2,tab_critere)

	print, ' ****** Vl1 *******'
	print, 'Median : '+ strtrim(err_1[6],1)
	print, 'Median (par function) : ' + strtrim(median(Vl3_m),1)
	print, '-1 sigma : '+ strtrim(err_1[6]-err_1[5],1)
	print, '+1 sigma : '+ strtrim(err_1[7]-err_1[6],1)
	print, '-2 sigma : '+ strtrim(err_1[6]-err_1[4],1)
	print, '+2 sigma : '+ strtrim(err_1[8]-err_1[6],1)
	print, ' ************'

	y_max=max(hist_1[1,*]/n_elements(Vl3_m))*1.1
	;e=write_on_ps_on(output_dir+'PDF_Vl3_m')
	nimp,name=output_dir+'PDF_Vl3.eps',/paper,/eps ; Savita's code
		plot, hist_1[0,*], smooth(hist_1[1,*]/n_elements(Vl3_m),2), psym=10,charsize=2,$
			xtitle=textoidl('V^2_{l=3}'), ytitle='Probability density',$
			;xtitle=textoidl('V^2_{l=1}'), ytitle='Probability density',$
			background=fsc_color('White'), color=fsc_color('Black'),yr=[0, y_max],/yst, xr=[0,1]
		;bmax=interpol(hist_1[1,*], hist_1[0,*], err_1[6])/n_elements(Vl3_m)
		;plots, [err_1[6], err_1[6]], [0, bmax],color=fsc_color('Dark Gray'), linestyle=2, thick=4
		bmax=interpol(hist_1[1,*], smooth(hist_1[0,*],2), err_1[7])/n_elements(Vl3_m)
		plots, [err_1[7], err_1[7]], [0, bmax],color=fsc_color('Red'), linestyle=2, thick=6
		bmax=interpol(hist_1[1,*], smooth(hist_1[0,*],2), err_1[5])/n_elements(Vl3_m)
		plots, [err_1[5], err_1[5]], [0, bmax],color=fsc_color('Red'), linestyle=2, thick=6
	fimp


end

pro plot_fullhist_2D, x2D_A, x2D_B,hist2D, hist_A,  hist_B, stat_A, stat_B, label_A, label_B, smooth_coef, $
	range_A=range_A, range_B=range_B, lines=lines, xlines=xlines, ylines=ylines, projected=projected, label_proj=label_proj, range_proj=range_proj, stat_proj=stat_proj, $
	spectro_inc_rot=spectro_inc_rot, legend_precision=legend_precision

if n_elements(legend_precision) eq 0 then legend_precision=['(f8.2)', '(f8.1)', '(f8.2)'] ; a1, inc, a1*sin(inc) optimized by default

if n_elements(range_A) lt 2 then range_A=[min(hist_A[0,*]),max(hist_A[0,*])]
if n_elements(range_B) lt 2 then range_B=[min(hist_B[0,*]),max(hist_B[0,*])]
if n_elements(lines) eq 0 then begin
    line=0
    coordlines=-1
endif else begin
    if lines[0] eq -1 then begin
      line=0
      coordlines=-1
    endif else begin
      coordlines=lines
      line=1
    endelse
endelse
if n_elements(xlines) eq 0 then begin
    xxline=0
    xcoordlines=-1
endif else begin
    if xlines[0] eq -1 then begin
      xxline=0
      xcoordlines=-1
    endif else begin
      xcoordlines=xlines
      xxline=1
    endelse
endelse
if n_elements(ylines) eq 0 then begin
    yline=0
    ycoordlines=-1
endif else begin
    if ylines[0] eq -1 then begin
      yline=0
      ycoordlines=-1
    endif else begin
      ycoordlines=ylines
      yline=1
    endelse
endelse
if n_elements(projected) eq 0 then begin
	projected=-1
	label_proj=-1
	range_proj=-1
	stat_proj=-1
endif
if n_elements(label_proj) eq 0 then label_proj='Label not provided!'
if n_elements(range_proj) lt 2 then range_proj=[min(projected[0,*]),max(projected[0,*])]
if n_elements(stat_proj) eq 0 then stat_proj=dblarr(10)
if n_elements(spectro_inc_rot) eq 0 then spectro_inc_rot=-1

!p.multi=[0,2,2]
;smooth_coef=4
size_char=0.75
size_char_axes=1.1

	b_max=max(hist_B[1,*])*1.1 & b_max2=b_max

	plot, hist_B[0,*], smooth(hist_B[1,*],smooth_coef), color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		xtitle=label_B,$
		psym=10, yr=[0, b_max], xr=range_B, /xst, /yst, ytickname=replicate(' ', 24), $
		pos=[0.45, 0.175, 0.9, 0.425]

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[6])
	plots, [stat_B[6],stat_B[6]], [0, b_max] , color=fsc_color('Green'), thick=5. ; median
	xyouts, stat_B[6]+(stat_B[1]*0.01), b_max*(1-0.82), strtrim(string(stat_B[6], format=legend_precision[1]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[5])
	plots, [stat_B[5],stat_B[5]],[0, b_max],  color=fsc_color('Orange'), thick=5.,linestyle=2 ; -1sigma
	xyouts, stat_B[5]+(stat_B[1]*0.03	), b_max2*(1-0.90), strtrim(string(stat_B[5], format=legend_precision[1]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[7])
	plots,  [stat_B[7],stat_B[7]],[0, b_max], color=fsc_color('Orange'), thick=5.,linestyle=2 ; +1sigma
	xyouts, stat_B[7]+(stat_B[1]*0.03), b_max2*(1-0.90), strtrim(string(stat_B[7], format=legend_precision[1]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=max(hist_A[1,*])*1.1 & b_max2=b_max
	;b_max=max(hist_A[1,*])*0.20 & b_max2=b_max ; FOR AD HOC MODIFICATIONS, USE THIS LINE

	plot, smooth(hist_A[1,*],smooth_coef), hist_A[0,*], color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		ytitle=label_A,$
		psym=10, xr=[ b_max, 0], yr=range_A,/xst, /yst, xtickname=replicate(' ', 24), $
		pos=[0.175, 0.45, 0.425, 0.9]

	b_max=interpol(smooth(hist_A[1,*],smooth_coef), hist_A[0,*], stat_A[6])
	plots, [0, b_max], [stat_A[6],stat_A[6]], color=fsc_color('Green'), thick=5. ; median
	xyouts, b_max*(1-0.55), stat_A[6]+(stat_A[1]*0.04), strtrim(string(stat_A[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_A[1,*],smooth_coef), hist_A[0,*], stat_A[5])
	plots, [0, b_max], [stat_A[5],stat_A[5]], color=fsc_color('Orange'), thick=5.,linestyle=2 ; -1sigma
	xyouts, b_max2*(1-0.87), stat_A[5]+(stat_A[1]*0.05), strtrim(string(stat_A[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_A[1,*],smooth_coef), hist_A[0,*], stat_A[7])
	plots, [0, b_max], [stat_A[7],stat_A[7]], color=fsc_color('Orange'), thick=5.,linestyle=2 ; +1sigma
	;xyouts, stat_A[7]+(stat_A[1]*0.01), b_max*(1+0.15), strtrim(string(stat_A[7], format='(f8.3)'),1), color=fsc_color('Orange'),charsize=size_char
	xyouts, b_max2*(1-0.87), stat_A[7]+(stat_A[1]*0.05), strtrim(string(stat_A[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	contour, smooth(hist2D,smooth_coef)/max(smooth(hist2D,smooth_coef)),x2D_B,x2D_A,nlevel=50,/fill, /nodata,$
	  yr=range_A,/xst,/yst,$
	  xr=range_B,$
	  xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  background=fsc_color('White'),color=fsc_color('Black'),charsize=size_char_axes, pos=[0.45, 0.45, 0.9, 0.9]
  polyfill, [range_B[0], range_B[1], range_B[1], range_B[0]], [range_A[0],range_A[0], range_A[1],range_A[1]], color=3, /data ; color=3 should be ok for a bckgrd at 0.05 ;fsc_color('blue'),/data
	contour, smooth(hist2D,smooth_coef)/max(smooth(hist2D,smooth_coef)) +0.05,x2D_B,x2D_A,nlevel=50,/fill,$
		yr=range_A,/xst,/yst,$
		xr=range_B,$
		xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
		background=fsc_color('White'),color=fsc_color('Black'),charsize=size_char_axes, pos=[0.45, 0.45, 0.9, 0.9], overplot=1
 
	if spectro_inc_rot[0] ne -1 then begin
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[3,*], color=fsc_color('Dark Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[2,*], color=fsc_color('Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[4,*], color=fsc_color('Gray'), thick=5, linestyle=2
	endif

  if line ne 0 then begin
    plots, [coordlines[0,0], coordlines[0,1]], [coordlines[0,2], coordlines[0,3]], color=fsc_color('Dark Gray'), thick=5, linestyle=4
    plots, [coordlines[1,0], coordlines[1,1]], [coordlines[1,2], coordlines[1,3]], color=fsc_color('Gray'), thick=5, linestyle=4
    plots, [coordlines[2,0], coordlines[2,1]], [coordlines[2,2], coordlines[2,3]], color=fsc_color('Gray'), thick=5, linestyle=4
  endif
  if yline ne 0 then begin
    plots, [ycoordlines[0,0], ycoordlines[0,1]], [ycoordlines[0,2], ycoordlines[0,3]], color=fsc_color('Dark Gray'), thick=5, linestyle=4
    plots, [ycoordlines[1,0], ycoordlines[1,1]], [ycoordlines[1,2], ycoordlines[1,3]], color=fsc_color('Gray'), thick=5, linestyle=4
    plots, [ycoordlines[2,0], ycoordlines[2,1]], [ycoordlines[2,2], ycoordlines[2,3]], color=fsc_color('Gray'), thick=5, linestyle=4
  endif
  if xxline ne 0 then begin
  	for kk=0, n_elements(xcoordlines[*,0])-1 do begin
  		plots, [xcoordlines[kk,0], xcoordlines[kk,1]], [xcoordlines[kk,2], xcoordlines[kk,3]], color=fsc_color('Dark Gray'), thick=5, linestyle=4
  	endfor
  endif

if projected[0] ne -1 then begin
		b_max=max(smooth(projected[1,*], smooth_coef))*1.02  & b_max2=b_max
		;b_max=max(smooth(projected[1,*], smooth_coef))*0.085  & b_max2=b_max  ; FOR AD HOC MODIFICATIONS, USE THIS LINE
		range_proj=[0,0.7]

		plot, projected[0,*], smooth(projected[1,*],smooth_coef), color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		xtitle=label_proj,$
		psym=10, yr=[0, b_max], xr=range_proj, /xst, /yst, ytickname=replicate(' ', 24), $
		pos=[0.175, 0.175, 0.425, 0.425]

		b_max=interpol(smooth(projected[1,*],smooth_coef), projected[0,*], stat_proj[6])
		plots, [stat_proj[6],stat_proj[6]], [0, b_max] , color=fsc_color('Green'), thick=5. ; median
		xyouts, stat_proj[6]+(stat_proj[1]*0.01), b_max*(1-0.82), strtrim(string(stat_proj[6], format=legend_precision[2]),1), color=fsc_color('Green'),charsize=size_char

		b_max=interpol(smooth(projected[1,*],smooth_coef), projected[0,*], stat_proj[5])
		plots, [stat_proj[5],stat_proj[5]],[0, b_max],  color=fsc_color('Orange'), thick=5.,linestyle=2 ; -1sigma
		xyouts, stat_proj[5]+(stat_proj[1]*0.03	), b_max2*(1-0.94), strtrim(string(stat_proj[5], format=legend_precision[2]),1), color=fsc_color('Orange'),charsize=size_char

		b_max=interpol(smooth(projected[1,*],smooth_coef), projected[0,*], stat_proj[7])
		plots,  [stat_proj[7],stat_proj[7]],[0, b_max], color=fsc_color('Orange'), thick=5.,linestyle=2 ; +1sigma
		xyouts, stat_proj[7]+(stat_proj[1]*0.03), b_max2*(1-0.94), strtrim(string(stat_proj[7], format=legend_precision[2]),1), color=fsc_color('Orange'),charsize=size_char

endif

!p.multi=0

;range_A=0 & range_B=0

end


; used to show 4 correlation map in the same plot. Specifically designed for
; the showing correlations between a1, a2, a3, inclination
pro plot_fullhist_2D_MATRIX, x2D_A, x2D_B, x2D_C, x2D_D, hist2D_AB, hist2D_AC, hist2D_AD, $
		hist2D_BC, hist2D_BD, hist2D_CD, hist_A,  hist_B, hist_C, hist_D, $
		stat_A, stat_B, stat_C, stat_D, label_A, label_B, label_C, label_D, smooth_coef, $
		range_A=range_A, range_B=range_B, range_C=range_C, range_D=range_D, $
		A_lines=A_lines,  B_lines=B_lines,  $ ;, seismo_only=seismo_only
		C_lines=C_lines,  D_lines=D_lines, $
		spectro_inc_rot=spectro_inc_rot, legend_precision=legend_precision

if n_elements(legend_precision) eq 0 then legend_precision=['(f8.2)', '(f8.1)', '(f8.2)'] ; a1, inc, a1*sin(inc) optimized by default

if n_elements(range_A) lt 2 then range_A=[min(hist_A[0,*]),max(hist_A[0,*])]
if n_elements(range_B) lt 2 then range_B=[min(hist_B[0,*]),max(hist_B[0,*])]
if n_elements(range_C) lt 2 then range_C=[min(hist_C[0,*]),max(hist_C[0,*])]
if n_elements(range_D) lt 2 then range_D=[min(hist_D[0,*]),max(hist_D[0,*])]

if n_elements(spectro_inc_rot) eq 0 then spectro_inc_rot=-1

!p.multi=[0,2,2]
;smooth_coef=4
size_char=0.75
size_char_axes=1.25

N=4 ; number of variables plotted, here a1,a2,a3,inc
xymin=0.05 ; 0.175
xymax=0.98 ; 0.9
interspace=0.01 ; 0.025
Delta=xymax-xymin - (N-1)*interspace
Width=Delta/N

!p.multi=[0,N,N]
; --------------------------------------------
; -------- Filling lower line of PDF ---------
; -------------------------------------------- 
	b_max=max(hist_A[1,*])*1.1 & b_max2=b_max
	plot, hist_A[0,*], smooth(hist_A[1,*],smooth_coef), color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		xtitle=label_A,$
		psym=10, xr=range_A, yr=[0, b_max],/xst, /yst, $ ;xtickname=replicate(' ', 24), $
		pos=[xymin+Width+interspace, xymin, xymin+2*Width+interspace, xymin+Width]

	b_max=interpol(smooth(hist_A[1,*],smooth_coef), hist_A[0,*], stat_A[6])
	plots, [stat_A[6],stat_A[6]], [0, b_max], color=fsc_color('Green'), thick=4. ; median
	;xyouts, b_max*(1-0.55), stat_A[6]+(stat_A[1]*0.04), strtrim(string(stat_A[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_A[1,*],smooth_coef), hist_A[0,*], stat_A[5])
	plots, [stat_A[5],stat_A[5]],[0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, b_max2*(1-0.87), stat_A[5]+(stat_A[1]*0.05), strtrim(string(stat_A[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_A[1,*],smooth_coef), hist_A[0,*], stat_A[7])
	plots, [stat_A[7],stat_A[7]], [0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, b_max2*(1-0.87), stat_A[7]+(stat_A[1]*0.05), strtrim(string(stat_A[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	; ----------------

	b_max=max(hist_B[1,*])*1.1 & b_max2=b_max
	plot, hist_B[0,*], smooth(hist_B[1,*],smooth_coef), color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		xtitle=label_B,$
		psym=10, yr=[0, b_max], xr=range_B, /xst, /yst, ytickname=replicate(' ', 24), $
		pos=[xymin+2*Width+2*interspace, xymin, xymin+3*Width+2*interspace, xymin+Width]
		
	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[6])
	plots, [stat_B[6],stat_B[6]], [0, b_max] , color=fsc_color('Green'), thick=4. ; median
	;xyouts, stat_B[6]+(stat_B[1]*0.01), b_max*(1-0.82), strtrim(string(stat_B[6], format=legend_precision[1]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[5])
	plots, [stat_B[5],stat_B[5]],[0, b_max],  color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, stat_B[5]+(stat_B[1]*0.03	), b_max2*(1-0.90), strtrim(string(stat_B[5], format=legend_precision[1]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[7])
	plots,  [stat_B[7],stat_B[7]],[0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, stat_B[7]+(stat_B[1]*0.03), b_max2*(1-0.90), strtrim(string(stat_B[7], format=legend_precision[1]),1), color=fsc_color('Orange'),charsize=size_char

	; ----------------

	b_max=max(hist_C[1,*])*1.1 & b_max2=b_max
	plot, hist_C[0,*], smooth(hist_C[1,*],smooth_coef), color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		xtitle=label_C,$
		psym=10, yr=[0, b_max], xr=range_C, /xst, /yst, ytickname=replicate(' ', 24), $
		pos=[xymin+3*Width+3*interspace, xymin, xymin+4*Width+3*interspace, xymin+Width]
		
	b_max=interpol(smooth(hist_C[1,*],smooth_coef), hist_C[0,*], stat_C[6])
	plots, [stat_C[6],stat_C[6]], [0, b_max] , color=fsc_color('Green'), thick=4. ; median
	;xyouts, stat_C[6]+(stat_C[1]*0.01), b_max*(1-0.82), strtrim(string(stat_C[6], format=legend_precision[1]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_C[1,*],smooth_coef), hist_C[0,*], stat_C[5])
	plots, [stat_C[5],stat_C[5]],[0, b_max],  color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, stat_C[5]+(stat_C[1]*0.03	), b_max2*(1-0.90), strtrim(string(stat_C[5], format=legend_precision[1]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_C[1,*],smooth_coef), hist_C[0,*], stat_C[7])
	plots,  [stat_C[7],stat_C[7]],[0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, stat_C[7]+(stat_C[1]*0.03), b_max2*(1-0.90), strtrim(string(stat_C[7], format=legend_precision[1]),1), color=fsc_color('Orange'),charsize=size_char

; ---------------------------------------------
; -------- Filling left column of PDF ---------
; ---------------------------------------------
	b_max=max(hist_D[1,*])*1.1 & b_max2=b_max
	plot, smooth(hist_D[1,*],smooth_coef), hist_D[0,*], color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		ytitle=label_D,$
		psym=10, xr=[b_max, 0], yr=range_D,/xst, /yst, $ ;, xtickname=replicate(' ', 24), $
		pos=[xymin, xymin+Width+interspace, xymin+Width, xymin+2*Width+interspace]
		
	b_max=interpol(smooth(hist_D[1,*],smooth_coef), hist_D[0,*], stat_D[6])
	plots, [0, b_max], [stat_D[6],stat_D[6]], color=fsc_color('Green'), thick=4. ; median
	;xyouts, b_max*(1-0.55), stat_D[6]+(stat_D[1]*0.04), strtrim(string(stat_D[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_D[1,*],smooth_coef), hist_D[0,*], stat_D[5])
	plots, [0, b_max], [stat_D[5],stat_D[5]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, b_max2*(1-0.87), stat_D[5]+(stat_D[1]*0.05), strtrim(string(stat_D[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_D[1,*],smooth_coef), hist_D[0,*], stat_D[7])
	plots, [0, b_max], [stat_D[7],stat_D[7]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, b_max2*(1-0.87), stat_D[7]+(stat_D[1]*0.05), strtrim(string(stat_D[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	; ----------------
	b_max=max(hist_C[1,*])*1.1 & b_max2=b_max
	plot, smooth(hist_C[1,*],smooth_coef), hist_C[0,*], color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		ytitle=label_C,$
		psym=10, xr=[b_max, 0], yr=range_C,/xst, /yst, xtickname=replicate(' ', 24), $
		pos=[xymin, xymin+2*Width+2*interspace, xymin+Width, xymin+3*Width+2*interspace]

	b_max=interpol(smooth(hist_C[1,*],smooth_coef), hist_C[0,*], stat_C[6])
	plots, [0, b_max], [stat_C[6],stat_C[6]], color=fsc_color('Green'), thick=4. ; median
	;xyouts, b_max*(1-0.55), stat_C[6]+(stat_D[1]*0.04), strtrim(string(stat_D[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_C[1,*],smooth_coef), hist_C[0,*], stat_C[5])
	plots, [0, b_max], [stat_C[5],stat_C[5]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, b_max2*(1-0.87), stat_D[5]+(stat_D[1]*0.05), strtrim(string(stat_D[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_C[1,*],smooth_coef), hist_C[0,*], stat_C[7])
	plots, [0, b_max], [stat_C[7],stat_C[7]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, b_max2*(1-0.87), stat_D[7]+(stat_D[1]*0.05), strtrim(string(stat_D[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	; ----------------
	b_max=max(hist_B[1,*])*1.1 & b_max2=b_max
	plot, smooth(hist_B[1,*],smooth_coef), hist_B[0,*], color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		ytitle=label_B,$
		psym=10, xr=[b_max, 0], yr=range_B,/xst, /yst, xtickname=replicate(' ', 24), $
		pos=[xymin, xymin+3*Width+3*interspace, xymin+Width, xymin+4*Width+ 3*interspace]

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[6])
	plots, [0, b_max], [stat_B[6],stat_B[6]], color=fsc_color('Green'), thick=4. ; median
	;xyouts, b_max*(1-0.55), stat_C[6]+(stat_D[1]*0.04), strtrim(string(stat_D[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[5])
	plots, [0, b_max], [stat_B[5],stat_B[5]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, b_max2*(1-0.87), stat_D[5]+(stat_D[1]*0.05), strtrim(string(stat_D[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hist_B[1,*],smooth_coef), hist_B[0,*], stat_B[7])
	plots, [0, b_max], [stat_B[7],stat_B[7]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, b_max2*(1-0.87), stat_D[7]+(stat_D[1]*0.05), strtrim(string(stat_D[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

    ; ------------------

; -----------------------------------------------
; -------- Filling first line of 2D PDF ---------
; -----------------------------------------------
	contour, smooth(hist2D_AD,smooth_coef)/max(smooth(hist2D_AD,smooth_coef)),x2D_A,x2D_D,nlevel=50,/fill, /nodata,$
	  yr=range_D,/xst,/yst,$
	  xr=range_A,$
	  xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
	  	pos=[xymin+Width+interspace, xymin+Width+interspace, xymin+2*Width+interspace, xymin+2*Width+ interspace]
	polyfill, [range_A[0],range_A[0], range_A[1],range_A[1]], [range_D[0], range_D[1], range_D[1], range_D[0]], color=3 ;fsc_color('Blue'), /data ; color=3 should be ok for a bckgrd at 0.05 ;fsc_color('blue'),/data
	contour, smooth(hist2D_AD,smooth_coef)/max(smooth(hist2D_AD,smooth_coef)) ,x2D_A,x2D_D,nlevel=50,/fill,$
		yr=range_D,/xst,/yst,$
		xr=range_A,$
		xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
		background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
		pos=[xymin+Width+interspace, xymin+Width+interspace, xymin+2*Width+interspace, xymin+2*Width+ interspace], overplot=1
	plot_lines, A_lines, x2D_A, x2D_D, 'vertical'
	plot_lines, D_lines, x2D_A, x2D_D, 'horizontal'

	contour, smooth(hist2D_BD,smooth_coef)/max(smooth(hist2D_AD,smooth_coef)),x2D_B,x2D_D,nlevel=50,/fill, /nodata,$
	  yr=range_D,/xst,/yst,$
	  xr=range_B,$
	  xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
	  	pos=[xymin+2*Width+2*interspace, xymin+Width+interspace, xymin+3*Width+2*interspace, xymin+2*Width+ interspace]
	polyfill, [range_B[0],range_B[0], range_B[1],range_B[1]], [range_D[0], range_D[1], range_D[1], range_D[0]], color=3 ;fsc_color('Blue'), /data ; color=3 should be ok for a bckgrd at 0.05 ;fsc_color('blue'),/data
	contour, smooth(hist2D_BD,smooth_coef)/max(smooth(hist2D_BD,smooth_coef)) ,x2D_B,x2D_D,nlevel=50,/fill,$
		yr=range_D,/xst,/yst,$
		xr=range_B,$
		xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
		background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
		pos=[xymin+2*Width+2*interspace, xymin+Width+interspace, xymin+3*Width+2*interspace, xymin+2*Width + interspace], overplot=1
	plot_lines, B_lines, x2D_B, x2D_D, 'vertical'
	plot_lines, D_lines, x2D_B, x2D_D, 'horizontal'

 	
	contour, smooth(hist2D_CD,smooth_coef)/max(smooth(hist2D_CD,smooth_coef)),x2D_C,x2D_D,nlevel=50,/fill, /nodata,$
	  yr=range_D,/xst,/yst,$
	  xr=range_C,$
	  xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
	  	pos=[xymin+3*Width+3*interspace, xymin+Width+interspace, xymin+4*Width+3*interspace, xymin+2*Width + interspace]
	polyfill, [range_C[0],range_C[0], range_C[1],range_C[1]], [range_D[0], range_D[1], range_D[1], range_D[0]], color=3 ;fsc_color('Blue'), /data ; color=3 should be ok for a bckgrd at 0.05 ;fsc_color('blue'),/data
	contour, smooth(hist2D_CD,smooth_coef)/max(smooth(hist2D_CD,smooth_coef)) ,x2D_C,x2D_D,nlevel=50,/fill,$
		yr=range_D,/xst,/yst,$
		xr=range_C,$
		xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
		background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
		pos=[xymin+3*Width+3*interspace, xymin+Width+interspace, xymin+4*Width+3*interspace, xymin+2*Width + interspace], overplot=1
	plot_lines, C_lines, x2D_C, x2D_D, 'vertical'
	plot_lines, D_lines, x2D_C, x2D_D, 'horizontal'


; -----------------------------------------------
; -------- Filling second line of 2D PDF ---------
; -----------------------------------------------
	contour, smooth(hist2D_AC,smooth_coef)/max(smooth(hist2D_AC,smooth_coef)),x2D_A,x2D_C,nlevel=50,/fill, /nodata,$
	  yr=range_C,/xst,/yst,$
	  xr=range_A,$
	  xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
	  	pos=[xymin+Width+interspace, xymin+2*Width+2*interspace, xymin+2*Width+interspace, xymin+3*Width+2*interspace]
	polyfill, [range_A[0],range_A[0], range_A[1],range_A[1]], [range_C[0], range_C[1], range_C[1], range_C[0]], color=3;fsc_color('blue') ;fsc_color('Blue'), /data ; color=3 should be ok for a bckgrd at 0.05 ;fsc_color('blue'),/data
	contour, smooth(hist2D_AC,smooth_coef)/max(smooth(hist2D_AC,smooth_coef)) ,x2D_A,x2D_C,nlevel=50,/fill,$
		yr=range_C,/xst,/yst,$
		xr=range_A,$
		xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
		background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
		pos=[xymin+Width+interspace, xymin+2*Width+2*interspace, xymin+2*Width+interspace, xymin+3*Width+2*interspace], overplot=1
 	plot_lines, A_lines, x2D_A, x2D_C, 'vertical'
	plot_lines, C_lines, x2D_A, x2D_C, 'horizontal'
 		  	    
    contour, smooth(hist2D_BC,smooth_coef)/max(smooth(hist2D_AD,smooth_coef)),x2D_B,x2D_C,nlevel=50,/fill, /nodata,$
	  yr=range_C,/xst,/yst,$
	  xr=range_B,$
	  xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
	  	pos=[xymin+2*Width+2*interspace, xymin+2*Width+2*interspace, xymin+3*Width+2*interspace, xymin+3*Width+ 2*interspace]
	polyfill, [range_B[0],range_B[0], range_B[1],range_B[1]], [range_C[0], range_C[1], range_C[1], range_C[0]], color=3;fsc_color('Blue'), /data ; color=3 should be ok for a bckgrd at 0.05 ;fsc_color('blue'),/data
	contour, smooth(hist2D_BC,smooth_coef)/max(smooth(hist2D_BC,smooth_coef)),x2D_B,x2D_C,nlevel=50,/fill,$
		yr=range_C,/xst,/yst,$
		xr=range_B,$
		xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
		background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
		pos=[xymin+2*Width+2*interspace, xymin+2*Width+2*interspace, xymin+3*Width+2*interspace, xymin+3*Width + 2*interspace], overplot=1
	plot_lines, B_lines, x2D_B, x2D_C, 'vertical'
	plot_lines, C_lines, x2D_B, x2D_C, 'horizontal'
 
; -----------------------------------------------
; -------- Filling third line of 2D PDF ---------
; -----------------------------------------------
	contour, smooth(hist2D_AB,smooth_coef)/max(smooth(hist2D_AB,smooth_coef)),x2D_A,x2D_B,nlevel=50,/fill, /nodata,$
	  yr=range_B,/xst,/yst,$
	  xr=range_A,$
	  xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
	  	pos=[xymin+Width+interspace, xymin+3*Width+3*interspace, xymin+2*Width+interspace, xymin+4*Width+3*interspace]
	polyfill, [range_A[0],range_A[0], range_A[1],range_A[1]], [range_B[0], range_B[1], range_B[1], range_B[0]], color=3; fsc_color('Blue'), /data ; color=3 should be ok for a bckgrd at 0.05 ;fsc_color('blue'),/data
	contour, smooth(hist2D_AB,smooth_coef)/max(smooth(hist2D_AB,smooth_coef)) ,x2D_A,x2D_B, nlevel=50,/fill,$
		yr=range_B,/xst,/yst,$
		xr=range_A,$
		xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
		background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
		pos=[xymin+Width+interspace, xymin+3*Width+3*interspace, xymin+2*Width+interspace, xymin+4*Width+3*interspace], overplot=1
	plot_lines, A_lines, x2D_A, x2D_B, 'vertical'
	plot_lines, B_lines, x2D_A, x2D_B, 'horizontal'
 
; ------------------------------------------------

	if spectro_inc_rot[0] ne -1 then begin
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[3,*], color=fsc_color('Dark Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[2,*], color=fsc_color('Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[4,*], color=fsc_color('Gray'), thick=5, linestyle=2
	endif

!p.multi=0

end


pro plot_fullhist_2D_Matrix_Nd, x2ds, hist2ds, hists, stats, labels, ranges=ranges, smooth_coef, $
		lines=lines, spectro_inc_rot=spectro_inc_rot, legend_precision=legend_precision
		
Nvars=n_elements(x2ds[*,0])

!p.multi=[0,Nvars+1,Nvars+1]

size_char=0.75
size_char_axes=1.25

xymin=0.05 ; 0.175
xymax=0.95 ; 0.9
interspace=0.005 ; 0.025
Delta=xymax-xymin - (Nvars-1)*interspace
Width=Delta/(Nvars)

;stop

for i=long(0), Nvars-1 do begin
	if n_elements(ranges) gt 0 then begin
		if ranges[i,0] eq ranges[i,1] then ranges[i,*]=[min(hists[i,0,*]),max(hists[i,0,*])]
	endif else begin
			if i eq 0 then ranges=dblarr(Nvars, 2)
			ranges[i,*]=[min(hists[i,0,*]),max(hists[i,0,*])]
	endelse
endfor

if n_elements(spectro_inc_rot) eq 0 then spectro_inc_rot=-1
;fimp

; --------------------------------------------
; -------- Filling upper line of PDF ---------
; -------------------------------------------- 
shift=1 ; used to correct from the fact that we do not need redondant plots...
xy0=xymin + width + interspace
;stop
for i=long(0), Nvars-2 do begin
			b_max=max(hists[i, 1,*])*1.1 & b_max2=b_max
			plot, hists[i, 0,*], smooth(hists[i,1,*],smooth_coef), color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
					;xtitle=labels[i,*], thick=2, $
					thick=2, $
					psym=10, xr=reform(ranges[i,*]), yr=[0, b_max],/xst, /yst, ytickname=replicate(' ', 24), $
					;pos=[xy0+i*(Width+interspace), xymin, xy0+(i+1)*Width+i*interspace, xymin+Width]
					pos=[xy0+i*(Width+interspace), xy0+(Nvars-1-shift)*(Width+interspace), xy0+(i+1)*Width+i*interspace, xy0+(Nvars-shift)*Width+(Nvars-1-shift)*interspace]
					Axis, xaxis=1, xtitle=labels[i,*], charsize=size_char_axes, color=fsc_color('Black'), /xst

			b_max=interpol(smooth(hists[i,1,*],smooth_coef), hists[i,0,*], stats[i, 6])
			plots, [stats[i,6],stats[i,6]], [0, b_max], color=fsc_color('Green'), thick=4. ; median
			;xyouts, b_max*(1-0.55), stat_A[6]+(stat_A[1]*0.04), strtrim(string(stat_A[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

			b_max=interpol(smooth(hists[i,1,*],smooth_coef), hists[i,0,*], stats[i, 5])
			plots, [stats[i,5],stats[i,5]],[0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
			;xyouts, b_max2*(1-0.87), stat_A[5]+(stat_A[1]*0.05), strtrim(string(stat_A[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

			b_max=interpol(smooth(hists[i,1,*],smooth_coef), hists[i,0,*], stats[i,7])
			plots, [stats[i,7],stats[i,7]], [0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
			;xyouts, b_max2*(1-0.87), stat_A[7]+(stat_A[1]*0.05), strtrim(string(stat_A[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char
endfor
;stop
; ---------------------------------------------
; -------- Filling left column of PDF ---------
; ---------------------------------------------

for j=long(shift), Nvars-1 do begin
	b_max=max(hists[j, 1,*])*1.1 & b_max2=b_max
	plot, smooth(hists[j,1,*],smooth_coef), hists[j, 0,*], color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		ytitle=labels[j,*], thick=2, $
		psym=10, xr=[b_max, 0], yr=ranges[j,*],/xst, /yst, xtickname=replicate(' ', 24), $
		;pos=[xymin, xy0+j*(Width+interspace), xymin+Width, xy0+(j+1)*Width+j*interspace]
		pos=[xymin, xymin+(j-shift)*(Width+interspace), xymin+Width, xymin+(j-shift+1)*Width+(j-shift)*interspace]

	b_max=interpol(smooth(hists[j, 1,*],smooth_coef), hists[j, 0,*], stats[j,6])
	plots, [0, b_max], [stats[j,6],stats[j,6]], color=fsc_color('Green'), thick=4. ; median
	;xyouts, b_max*(1-0.55), stat_D[6]+(stat_D[1]*0.04), strtrim(string(stat_D[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hists[j,1,*],smooth_coef), hists[j, 0,*], stats[j, 5])
	plots, [0, b_max], [stats[j, 5],stats[j, 5]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, b_max2*(1-0.87), stat_D[5]+(stat_D[1]*0.05), strtrim(string(stat_D[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hists[j, 1,*],smooth_coef), hists[j, 0,*], stats[j, 7])
	plots, [0, b_max], [stats[j, 7],stats[j, 7]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, b_max2*(1-0.87), stat_D[7]+(stat_D[1]*0.05), strtrim(string(stat_D[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char
endfor
;stop
; -----------------------------------------------
; ------ Filling lines/columns with 2D PDF ------
; -----------------------------------------------
for i=long(0), Nvars-1 do begin
	for j=long(0), Nvars-1 do begin
		if (i lt j) then begin
			contour, smooth(hist2Ds[i,j,*,*],smooth_coef)/max(smooth(hist2Ds[i,j,*,*],smooth_coef)),x2Ds[i,j,*],x2Ds[j,i,*],nlevel=50,/fill, /nodata,$
	  			yr=ranges[j,*],/xst,/yst,$
				xr=ranges[i,*],$
	  			xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  			background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
	  			;pos=[xy0+i*(Width+interspace), xy0+ j*(Width+interspace), xy0+(i+1)*Width+i*interspace, xy0+(j+1)*Width+ j*interspace]
				pos=[xy0+i*(Width+interspace), xymin+ (j-shift)*(Width+interspace), xy0+(i+1)*Width+i*interspace, xymin+(j-shift+1)*Width+ (j-shift)*interspace]	

			polyfill, [ranges[i,0],ranges[i,0], ranges[i,1],ranges[i,1]], [ranges[j,0], ranges[j,1], ranges[j,1], ranges[j,0]], color=3
			contour, smooth(hist2Ds[i,j,*,*],smooth_coef)/max(smooth(hist2Ds[i,j,*,*],smooth_coef)) ,x2Ds[i,j,*],x2Ds[j,i,*],nlevel=50,/fill,$
					yr=ranges[j,*],/xst,/yst,$
					xr=ranges[i,*],$
					xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
					background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
					pos=[xymin+i*(Width+interspace), xymin+ j*(Width+interspace), xymin+(i+1)*Width+i*interspace, xymin+(j+1)*Width+ j*interspace],$
					overplot=1
			;plot_lines, lines[i,*], x2D[i,j,*], x2D[j,i,*], 'vertical'
			;plot_lines, lines[j,*], x2D[i,j,*], x2D[j,i,*], 'horizontal'
		endif
	endfor
endfor

; ------------------------------------------------

	if spectro_inc_rot[0] ne -1 then begin
		print, "spectro_inc_rot is not implemented"
		print, "You need to think about how to put it here...."
		print, "The program will stop now"
		stop
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[3,*], color=fsc_color('Dark Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[2,*], color=fsc_color('Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[4,*], color=fsc_color('Gray'), thick=5, linestyle=2
	endif

!p.multi=0

;stop
end


pro plot_fullhist_2D_Matrix_Nd_v2, x2ds, hist2ds, hists, stats, labels, ranges=ranges, smooth_coef, $
		lines=lines, spectro_inc_rot=spectro_inc_rot, legend_precision=legend_precision, posplots=posplots, $
		hist_extra=hist_extra, err_extra=err_extra, label_extra=label_extra
		
Nvars=n_elements(x2ds[*,0])

!p.multi=[0,Nvars+1,Nvars+1]

;size_char=0.75
size_char_axes=1.9

ymin=0.01
xymin=0.10 ; 0.175
xymax=1.0 ; 0.9
interspace=0.01 ; 0.025
dpdf=0.05
Delta=xymax-xymin - (Nvars-1)*interspace
width=Delta/(Nvars) ; in y-axis

ratio=aspect(1.) ; gives coordinates centered on the window which correspond to a square plot... used to rescale the widths such that I have square plots
r=(ratio[2]-ratio[0])/(ratio[3]-ratio[1])
;r=(ratio[3]-ratio[1])/(ratio[2]-ratio[0])
length=width*r ; rescaled length in order to have squares
;stop

if n_elements(posplots) eq 0 then begin
	posplots=dblarr(2,Nvars)
	posplots[0,*]=findgen(Nvars)
	posplots[1,*]=findgen(Nvars)
endif

for i=long(0), Nvars-1 do begin
	if n_elements(ranges) gt 0 then begin
		if ranges[i,0] eq ranges[i,1] then ranges[i,*]=[min(hists[i,0,*]),max(hists[i,0,*])]
	endif else begin
			if i eq 0 then ranges=dblarr(Nvars, 2)
			ranges[i,*]=[min(hists[i,0,*]),max(hists[i,0,*])]
	endelse
endfor

if n_elements(spectro_inc_rot) eq 0 then spectro_inc_rot=-1
;fimp

; --------------------------------------------
; -------- Filling upper line of PDF ---------
; -------------------------------------------- 
shift=1 ;1 ; used to correct from the fact that we do not need redondant plots...
xy0=xymin + width + interspace
y0=ymin + width + interspace
;stop
for i=long(0), Nvars-2 do begin
			b_max=max(hists[posplots[0,i], 1,*])*1.1 & b_max2=b_max
			ticku='Numeric'
			plot, hists[posplots[0,i], 0,*], smooth(hists[posplots[0,i],1,*],smooth_coef), color=fsc_color('black'), background=fsc_color('white'), $
				charsize=size_char_axes, thick=2, xtickunits=ticku, $
				psym=10, xr=reform(ranges[posplots[0,i],*]), yr=[0, b_max],/xst, /yst, ytickname=replicate(' ', 24), xtickname=replicate(' ', 24), $
				pos=[xy0+i*(length+interspace), y0+(Nvars-1-shift)*(width+interspace), xy0+(i+1)*length+i*interspace, y0+(Nvars-shift)*width+(Nvars-1-shift)*interspace]
				Axis, xaxis=1, xtitle=labels[posplots[0,i],*], charsize=size_char_axes, color=fsc_color('Black'), /xst, xtickunits=ticku

			b_max=interpol(smooth(hists[posplots[0,i],1,*],smooth_coef), hists[posplots[0,i],0,*], stats[posplots[0,i], 6])
			plots, [stats[posplots[0,i],6],stats[posplots[0,i],6]], [0, b_max], color=fsc_color('Cyan'), thick=4. ; median
			;xyouts, b_max*(1-0.55), stat_A[6]+(stat_A[1]*0.04), strtrim(string(stat_A[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

			b_max=interpol(smooth(hists[posplots[0,i],1,*],smooth_coef), hists[posplots[0,i],0,*], stats[posplots[0,i], 5])
			plots, [stats[posplots[0,i],5],stats[posplots[0,i],5]],[0, b_max], color=fsc_color('Cyan'), thick=4.,linestyle=2 ; -1sigma
			;xyouts, b_max2*(1-0.87), stat_A[5]+(stat_A[1]*0.05), strtrim(string(stat_A[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

			b_max=interpol(smooth(hists[posplots[0,i],1,*],smooth_coef), hists[posplots[0,i],0,*], stats[posplots[0,i],7])
			plots, [stats[posplots[0,i],7],stats[posplots[0,i],7]], [0, b_max], color=fsc_color('Cyan'), thick=4.,linestyle=2 ; +1sigma
			;xyouts, b_max2*(1-0.87), stat_A[7]+(stat_A[1]*0.05), strtrim(string(stat_A[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

			if labels[posplots[0,i]] eq label_extra then begin
				oplot, hist_extra[0,*],  smooth(hist_extra[1,*]*max(hists[posplots[0,i],1,*])/max(hist_extra[1,*]),smooth_coef), color=fsc_color('Red')

				;b_max=interpol(smooth(hist_extra[1,*],smooth_coef), hist_extra[0,*], err_extra[6])
				;plots, [stats[posplots[0,i],6],stats[posplots[0,i],6]], [0, b_max], color=fsc_color('Green'), thick=4.,linestyle=2 ; median

				;b_max=interpol(smooth(hist_extra[1,*],smooth_coef), hist_extra[0,*], err_extra[5])
				;plots, [stats[posplots[0,i],5],stats[posplots[0,i],5]], [0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma

				;b_max=interpol(smooth(hist_extra[1,*],smooth_coef), hist_extra[0,*], err_extra[7])
				;plots, [stats[posplots[0,i],7],stats[posplots[0,i],7]], [0, b_max], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
			endif

endfor
;stop
; ---------------------------------------------
; -------- Filling left column of PDF ---------
; ---------------------------------------------
for j=long(shift), Nvars-1 do begin
	b_max=max(hists[posplots[1,j], 1,*])*1.1 & b_max2=b_max
	ticku='Numeric'
	plot, smooth(hists[posplots[1,j],1,*],smooth_coef), hists[posplots[1,j], 0,*], color=fsc_color('black'), background=fsc_color('white'),charsize=size_char_axes, $
		ytitle=labels[posplots[1,j],*], thick=2, xtickunits=ticku, $
		psym=10, xr=[b_max, 0], yr=ranges[posplots[1,j],*],/xst, /yst, xtickname=replicate(' ', 24), $
		pos=[xymin+dpdf, ymin+(j-shift)*(width+interspace), xymin+dpdf+length, ymin+(j-shift+1)*width+(j-shift)*interspace]
	
	b_max=interpol(smooth(hists[posplots[1,j], 1,*],smooth_coef), hists[posplots[1,j], 0,*], stats[posplots[1,j],6])
	plots, [0, b_max], [stats[posplots[1,j],6],stats[posplots[1,j],6]], color=fsc_color('Cyan'), thick=4. ; median
	;xyouts, b_max*(1-0.55), stat_D[6]+(stat_D[1]*0.04), strtrim(string(stat_D[6], format=legend_precision[0]),1), color=fsc_color('Green'),charsize=size_char

	b_max=interpol(smooth(hists[posplots[1,j],1,*],smooth_coef), hists[posplots[1,j], 0,*], stats[posplots[1,j], 5])
	plots, [0, b_max], [stats[posplots[1,j], 5],stats[posplots[1,j], 5]], color=fsc_color('Cyan'), thick=4.,linestyle=2 ; -1sigma
	;xyouts, b_max2*(1-0.87), stat_D[5]+(stat_D[1]*0.05), strtrim(string(stat_D[5], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	b_max=interpol(smooth(hists[posplots[1,j], 1,*],smooth_coef), hists[posplots[1,j], 0,*], stats[posplots[1,j], 7])
	plots, [0, b_max], [stats[posplots[1,j], 7],stats[posplots[1,j], 7]], color=fsc_color('Cyan'), thick=4.,linestyle=2 ; +1sigma
	;xyouts, b_max2*(1-0.87), stat_D[7]+(stat_D[1]*0.05), strtrim(string(stat_D[7], format=legend_precision[0]),1), color=fsc_color('Orange'),charsize=size_char

	if labels[posplots[1,j]] eq label_extra then begin
		oplot, smooth(hist_extra[1,*]*max(hists[posplots[0,j],1,*])/max(hist_extra[1,*]),smooth_coef), hist_extra[0,*], color=fsc_color('Red')
		;b_max=interpol(smooth(hist_extra[0,*],smooth_coef), hist_extra[1,*], err_extra[6])
		;plots, [0, b_max], [stats[posplots[0,i],6],stats[posplots[0,i],6]], color=fsc_color('Green'), thick=4.,linestyle=2 ; median

		;b_max=interpol(smooth(hist_extra[0,*],smooth_coef), hist_extra[1,*], err_extra[5])
		;plots, [0, b_max], [stats[posplots[0,i],5],stats[posplots[0,i],5]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; -1sigma

		;b_max=interpol(smooth(hist_extra[0,*],smooth_coef), hist_extra[1,*], err_extra[7])
		;plots, [0, b_max], [stats[posplots[0,i],7],stats[posplots[0,i],7]], color=fsc_color('Orange'), thick=4.,linestyle=2 ; +1sigma
	endif

endfor
;stop
; -----------------------------------------------
; ------ Filling lines/columns with 2D PDF ------
; -----------------------------------------------
for i=long(0), Nvars-1 do begin
	for j=long(0), Nvars-1 do begin
		if (i lt j) then begin
			contour, smooth(hist2Ds[posplots[0,i],posplots[1,j],*,*],smooth_coef)/max(smooth(hist2Ds[posplots[0,i],posplots[1,j],*,*],smooth_coef)),x2Ds[posplots[0,i],posplots[1,j],*],x2Ds[posplots[1,j],posplots[0,i],*],nlevel=50,/fill, /nodata,$
	  			yr=ranges[posplots[1,j],*],/xst,/yst,$
				xr=ranges[posplots[0,i],*],$
	  			xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
	  			background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
				pos=[xy0+i*(length+interspace), ymin+ (j-shift)*(width+interspace), xy0+(i+1)*length+i*interspace, ymin+(j-shift+1)*width+ (j-shift)*interspace]	

			polyfill, [ranges[posplots[0,i],0],ranges[posplots[0,i],0], ranges[posplots[0,i],1],ranges[posplots[0,i],1]], [ranges[posplots[1,j],0], ranges[posplots[1,j],1], ranges[posplots[1,j],1], ranges[posplots[1,j],0]], color=3
			contour, smooth(hist2Ds[posplots[0,i],posplots[1,j],*,*],smooth_coef)/max(smooth(hist2Ds[posplots[0,i],posplots[1,j],*,*],smooth_coef)) ,x2Ds[posplots[0,i],posplots[1,j],*],x2Ds[posplots[1,j],posplots[0,i],*],nlevel=50,/fill,$
					yr=ranges[posplots[1,j],*],/xst,/yst,$
					xr=ranges[posplots[0,i],*],$
					xtickname=replicate(' ', 24) ,ytickname=replicate(' ', 24),$
					background=fsc_color('Blue'),color=fsc_color('Black'),charsize=size_char_axes, $
					pos=[xymin+i*(length+interspace), ymin+ j*(width+interspace), xymin+(i+1)*length+i*interspace, ymin+(j+1)*width+ j*interspace],$
					overplot=1
			;plot_lines, lines[i,*], x2D[i,j,*], x2D[j,i,*], 'vertical'
			;plot_lines, lines[j,*], x2D[i,j,*], x2D[j,i,*], 'horizontal'
		endif
	endfor
endfor

; ------------------------------------------------

	if spectro_inc_rot[0] ne -1 then begin
		print, "spectro_inc_rot is not implemented"
		print, "You need to think about how to put it here...."
		print, "The program will stop now"
		stop
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[3,*], color=fsc_color('Dark Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[2,*], color=fsc_color('Gray'), thick=5, linestyle=2
		oplot, spectro_inc_rot[0,*], spectro_inc_rot[4,*], color=fsc_color('Gray'), thick=5, linestyle=2
	endif

!p.multi=0

;stop
end

; either number or samples
function seismic_vsini, nusini, Dnu, Teff, nu_max

	Teff_sun= 5777d
	Dnu_sun=135.1d
	nu_max_sun=3150d
	R_sun=6.96342d5 ; in km

	R= (1d*nu_max/nu_max_sun) * (1d*Dnu/Dnu_sun)^(-2d) * (1d*Teff/Teff_sun)^(1d/2)
	vsini=2d * !pi * R * R_sun * nusini * 1d-6

return, vsini
end

function vsini_pdf, nusini_samples, Dnu, Teff, nu_max, err_Dnu, err_Teff, err_nu_max

	N=n_elements(nusini_samples)

	Dnu_samples=gauss_random(Dnu, err_Dnu, N)
	Teff_samples=gauss_random(Teff, err_Teff, N)
	nu_max_samples=gauss_random(nu_max, err_nu_max, N)

	vsini=seismic_vsini(nusini_samples, Dnu_samples, Teff_samples, nu_max_samples)

return, vsini
end


function vsini_pdf_Radius, nusini_samples, Radius, err_Radius
	R_sun=6.96342d5 ; in km

	N=n_elements(nusini_samples)

	R=gauss_random(Radius, err_Radius, N)
	vsini=2d * !pi * R * R_sun * nusini_samples * 1d-6


return, vsini
end


function nu_rot_from_spec_vsini, Radius, spectro_vsini, err_spectro_vsini
	R_sun=6.96342d5 ; in km

	inc=findgen(91)*!pi/180 ; values from 0 to !pi/2 with step of 1 degree
	Period=2d * !pi * Radius * R_sun * sin(inc)/spectro_vsini
	nu_rot=1d6/period

  if n_elements(err_spectro_vsini) eq 1 then begin
    err_m=err_spectro_vsini
    err_p=err_spectro_vsini
  endif
  if n_elements(err_spectro_vsini) eq 2 then begin
    err_m=err_spectro_vsini[0]
    err_p=err_spectro_vsini[1]
  endif
  if n_elements(err_spectro_vsini) gt 2 then begin
    PRINT, 'err_spectro_vsini has more than 2 elements... cannot proceed!'
    PRINT, 'Emmergency stop'
    stop
  endif
	Period=2d * !pi * Radius * R_sun * sin(inc)/(spectro_vsini - 2d *err_m)
	nu_rot_2sigm=1d6/period

	Period=2d * !pi * Radius * R_sun * sin(inc)/(spectro_vsini - err_m)
	nu_rot_1sigm=1d6/period

	Period=2d * !pi * Radius * R_sun * sin(inc)/(spectro_vsini  + err_p)
	nu_rot_1sigp=1d6/period

	Period=2d * !pi * Radius * R_sun * sin(inc)/(spectro_vsini + 2d * err_p)
	nu_rot_2sigp=1d6/period

out=dblarr(6, n_elements(inc))
out[0,*]=inc * 180d / !pi
out[1,*]=nu_rot_2sigm
out[2,*]=nu_rot_1sigm
out[3,*]=nu_rot
out[4,*]=nu_rot_1sigp
out[5,*]=nu_rot_2sigp

return, out
end

; generate random numbers of std dev err_x and of mean x
function gauss_random, x, err_x, N

	r= x + randomn(seed, N)*err_x

return, r
end

; use individual frequencies and amplitudes to compute nu_max
function weighted_numax, file

 restore, file ; synthese file

frequencies=reform(stat_synthese_freq[3, 0, *])
frequencies=[frequencies, reform(stat_synthese_freq[3, 1, *])]

 amplitudes=reform(stat_synthese_amplitude[3, 0, *])
 amplitudes=[amplitudes, reform(stat_synthese_amplitude[3, 1, *])]

 nu_max=WTD_MEAN(frequencies, amplitudes)

return, nu_max
end


pro Hist_2D_NICE, param1, param2, Name_1, Name_2, ps_file_out

	tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma
	!p.thick=3

	nu=param1
	nus=param2
	Name_nu=Name_1
	Name_nus=Name_2

	file_out=ps_file_out

	resol1=( max(nu) - min(nu) ) / 20
	resol2=( max(nus) - min(nus) ) / 20

	N1=n_elements(nu)
	N2=n_elements(nus)

	if N1 eq N2 then Nb_classes=100*N1/10000. else begin
		print, 'BIZARRERIE !!! Tailles des tables diff�rentes. ARRET PREMATRUEE'
		stop
	endelse


	hist = HIST_2D( nu, nus,bin1=resol1,bin2=resol2,min1=min(nu),min2=min(nus))/(1.*N1)

	x_nu=max([n_elements(hist[*,0]), n_elements(hist[0,*])])

	hist=congrid(hist, x_nu, x_nu)

	x1=findgen(n_elements(hist[*,0]))*(max(nu)-min(nu))/(n_elements(hist[*,0])-1)+min(nu)
	x2=findgen(n_elements(hist[0,*]))*(max(nus)-min(nus))/(n_elements(hist[0,*])-1)+min(nus)

	k=30
	N2=k*n_elements(hist[0,*])
	pdf0=congrid(hist, N2, N2) ; VIRTUALY increase the precision of the map

	pdf=smooth(pdf0, 52,/edge_truncate)
	pdf=pdf
	pdf=pdf/total(pdf)
	pdf_show=smooth(pdf0, 30, /edge_truncate)
	pdf_show=pdf_show
	pdf_show=pdf_show/total(pdf_show)

	x_int=rebin(x1, k*n_elements(x1))
	y_int=rebin(x2, k*n_elements(x2))


	histo1D_1=HISTOGRAM(nu,nbins=0.5*N1/Nb_classes)
    x_hist1D_1=findgen(n_elements(histo1D_1))*(max(nu)-min(nu))/ $
    	(n_elements(histo1D_1)-1)+min(nu)
	err_1=estimate_1_sigma_error( x_hist1D_1,histo1D_1,68.3,2,tab_critere)

	h1D_1=dblarr(2, n_elements(histo1D_1))
	h1D_1[0,*]=x_hist1D_1 & h1D_1[1,*]=histo1D_1

    histo1D_2=HISTOGRAM(nus,nbins=1.*N1/Nb_classes)
    x_hist1D_2=findgen(n_elements(histo1D_2))*(max(nus)-min(nus))/ $
    	(n_elements(histo1D_2)-1)+min(nus)
	err_2=estimate_1_sigma_error( x_hist1D_2,histo1D_2,68.3,2,tab_critere)

	h1D_2=dblarr(2, n_elements(histo1D_2))
	h1D_2[0,*]=x_hist1D_2 & h1D_2[1,*]=histo1D_2

	loadct,27
	e=write_on_ps_on(file_out)
	plot_fullhist_2D, x_int, y_int,transpose(pdf_show)+ 0.001+ (max(transpose(pdf_show))-min(transpose(pdf_show)))*0.1, h1D_1,  h1D_2, err_1, err_2, Name_nu, Name_nus, 1, $
		range_A=0, range_B=0, lines=lines;, projected=hist_vsini, label_proj=Name_proj, range_proj=0, stat_proj=err_vsini
	e=write_on_ps_off('')
end


pro Corel_nu_SplitHeightGamma, dir, i0, imax, incr, i_incl,i0_l2, imax_l2, i0_l3, imax_l3
	tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma
	!p.thick=3


;	dir= 'C:\Work_dir\Analysis_Results\8751420_Mtmp\8751420_l2free\Files\'
;	incr=18
;	i0=110
;	imax=127

;	dir= 'C:\Work_dir\Analysis_Results\9574283_Mtmp\Files\'
;	incr=17
;	i0=109
;	imax=125



	for i=i0, imax do begin

		name=read_sav(dir, i)
		restore, name
		nu=param
		Name_nu='l=1 Frequency'

		name=read_sav(dir, i+2*incr)
		restore, name
		nus=param
		Name_nus='l=1 Splitting'

		file_out=dir+'SCorel_'+strtrim(i,1)+'-'+strtrim(i+2*incr,1)
		Hist_2D_NICE, nu, nus, Name_nu, Name_nus, file_out
	endfor

	for i=i0, imax do begin

		name=read_sav(dir, i)
		restore, name
		nu=param
		Name_nu='l=1 Frequency'

		name=read_sav(dir, i+incr)
		restore, name
		nus=param
		Name_nus='l=1 Width'

		file_out=dir+'WCorel_'+strtrim(i,1)+'-'+strtrim(i+incr,1)
		Hist_2D_NICE, nu, nus, Name_nu, Name_nus, file_out
	endfor

	for i=i0, imax do begin

		name=read_sav(dir, i)
		restore, name
		nu=param
		Name_nu='l=1 Frequency'

		name=read_sav(dir, i-incr)
		restore, name
		nus=param
		Name_nus='l=1 Height'

		file_out=dir+'HCorel_'+strtrim(i,1)+'-'+strtrim(i-incr,1)
		Hist_2D_NICE, nu, nus, Name_nu, Name_nus, file_out
	endfor

	for i=i0, imax do begin

		name=read_sav(dir, i_incl)
		restore, name
		nu=param
		Name_nu='Inclination (deg)'

		name=read_sav(dir, i+2*incr)
		restore, name
		nus=param
		Name_nus='l=1 Splitting'

		file_out=dir+'IncCorel_l1_'+strtrim(i+2*incr,1)+'-'+strtrim(i_incl,1)
		Hist_2D_NICE, nu, nus, Name_nu, Name_nus, file_out
	endfor

	if n_elements(i0_l2) ne 0 then for i=i0_l2, imax_l2 do begin

		name=read_sav(dir, i_incl)
		restore, name
		nu=param
		Name_nu='Inclination (deg)'

		name=read_sav(dir, i)
		restore, name
		nus=param
		Name_nus='l=2 Splitting'

		file_out=dir+'IncCorel_l2_'+strtrim(i,1)+'-'+strtrim(i_incl,1)
		Hist_2D_NICE, nu, nus, Name_nu, Name_nus, file_out
	endfor

	if n_elements(i0_l3) ne 0 then for i=i0_l3, imax_l3 do begin

		name=read_sav(dir, i_incl)
		restore, name
		nu=param
		Name_nu='Inclination (deg)'

		name=read_sav(dir, i)
		restore, name
		nus=param
		Name_nus='l=3 Splitting'

		file_out=dir+'IncCorel_l3_'+strtrim(i,1)+'-'+strtrim(i_incl,1)
		Hist_2D_NICE, nu, nus, Name_nu, Name_nus, file_out
	endfor
end



function read_sav, dir, i

	if i lt 10 then begin
 		name=dir+'param_00'+strtrim(i,1)+'.sav'
 	endif
	if i ge 10 AND i lt 100 then begin
		name=dir+'param_0'+strtrim(i,1)+'.sav'
	endif
	if i ge 100 then begin
		name=dir+'param_'+strtrim(i,1)+'.sav'
	endif

return, name
end

pro show_correl_standard_fit

; HatP7....
dir='C:\Work_dir\TOKYO_WORK\Exoplanets\Analysed_Stars\HatP7-2014-ID1-short\Files\'
core='param_'

filter_ind='061' ; inclination
thld_min=60d
thld_max=90d

nus_ind='035'

gamma_indi=48
gamma_indf=53

restore, dir+ core + filter_ind + '.sav'
filter=where(param ge thld_min AND param le thld_max)


restore, dir+core + nus_ind + '.sav'
nus=param[filter]

	for i=gamma_indi, gamma_indf do begin

		Name_nu='Splitting (microHz)'

		name=dir+core+ '0'+ strtrim(floor(i),1) + '.sav'
		restore, name
		nu_width=param[filter]
		Name_nus='l=0 Width (microHz)'

		file_out=dir+'SplitCorel_With_0'+strtrim(floor(i),1)+'-'+ nus_ind
		Hist_2D_NICE, nus, nu_width, Name_nu, Name_nus, file_out
	endfor

stop
end
