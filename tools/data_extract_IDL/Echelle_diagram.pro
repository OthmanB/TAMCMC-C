@read_getmodel_file
; updated to be compatible with cpp program version > 1.6
pro show_ech_diag_CPP, getmodel_file, data, ps=ps, shifts=shifts, trunc_spec=trunc_spec

struc=read_getmodel_file(getmodel_file)
els=struc.modeparams[*,0]

lmax=max(els)
Nf_el=dblarr(lmax+1)
for el=0, lmax do begin
	pos=where(els eq el)
	Nf_el[el]=n_elements(pos)
endfor
freq_tab_t0=dblarr(lmax+1, max(Nf_el))
for el=0, lmax do begin
	pos=where(els eq el)
	if pos[0] ne -1 then freq_tab_t0[el, 0:n_elements(pos)-1]=struc.modeparams[pos,1]
endfor

freq=reform(data[0,*])
spec_reg=reform(data[1,*])
noise=reform(data[2,*])

fit_freq=linfit(findgen(n_elements(freq)), freq) ; This fit is way much more precise than the difference of the first points (subject to round-off errors)
resol=fit_freq[1] ; fit_freq[0] is by definition min(freq0) and fit_freq[1] is the exact resolution

if resol gt 0.5 then smooth_coef=median(struc.modeparams[*,3])/resol
if resol le 0.5 then smooth_coef=median(struc.modeparams[*,3])/resol/2.

print, 'Applied smoothing coefficient (nbins): ', smooth_coef ,   '     equivalent to ', smooth_coef*resol, '   (microHz)' 

s0=spec_reg/noise

if n_elements(ps) eq 0 then ps=0
if n_elements(shifts) eq 0 then shifts=0 ; no shift by defaut
if n_elements(trunc_spec) eq 0 then limit=max(s0) else limit=trunc_spec
cor=1
print=0 ; To have Black/white data or Color
bar=0; by defaut we put the graduation bar
extended=0

;residu=continuity_condition(reform(freq_tab_t0[0,0:Nf_el[0]-1]),  n_elements(freq_tab_t0[0,0:Nf_el[0]-1]), n_elements(freq_tab_t0[0,0:Nf_el[0]-1]))
;residu_l1=continuity_condition(reform(freq_tab_t0[1,0:Nf_el[1]-1]),  n_elements(freq_tab_t0[1,0:Nf_el[[1]-1]), n_elements(freq_tab_t0[1,0:Nf_el[[1]-1]))

c=linfit(findgen(Nf_el[0]),reform(freq_tab_t0[0,0:Nf_el[0]-1]) ) 
Dnu=c[1]

r=min(freq_tab_t0[1,where(freq_tab_t0[1,*] ne 0)])/Dnu
n0=floor(r)
epsilon=r-n0
n0=n0-1

overplot=2 ; show diamonds only
color=['Orange', 'Red','Dark Gray', 'Yellow']

n_number=n_elements(freq_tab_t0[0,*])+2
n0=n0-2

if n_elements(freq) eq 0 then begin
	print, 'BEWARE: We didn t find the original power spectrum (files freq and spec_reg)'
	print, 'Trying to use the file_spec option...'
	if n_elements(file_spec) ne 0 then restore, file_spec else begin
		print, 'file_spec has not been specified ! Emmergency stop!'
		stop
	endelse
endif

if max(freq) lt 10 then x0=freq*1d6 else x0=freq
if min(x0) gt resol then begin ; The echelle diagram works only if freq and spec_reg are vectors such that freq[0]=0... ensure this is not the case

	Nadd=x0[0]/resol ; number of points that need to be added
	fadd=findgen(Nadd-1)*resol
	x0=[fadd[1:*], x0]
	mean_at0=mean(s0[0:5./resol]) ; average over 5 microHz of the spectrum values
	s0=[replicate(mean_at0, n_elements(fadd)), s0]
endif

if ps eq 0 then window, 0, xsize=1400, ysize=800

test=where(s0 ge limit)
if test[0] ne -1 then s0[test]=limit

;stop
echellecorot_v3,x0,s0,Dnu,n0+shifts,uu,freq_table=freq_tab_t0,smooth_coef=smooth_coef,ps=ps,overplot=overplot,$
			color=color,n_number=n_number, print=print, bar=bar


;stop
;level=0.999
;if smooth_coef ge 30 then smooth_coef=30
;threshold=proba_calc_smooth_signal(1d + 0.05, smooth_coef, level)


end


; procedure that overplot frequencies on the echelle diagram of the data
; freq : the frequencies
; spec_reg : the power spectrum
; delta : the large separation
; nbegin : first n index
; uu : an output of the echelle diagram
; freq_table : input table of frequencies to overplot
; smooth_coef : coeficient of smoothing of the power spectrum (default = 1, ie no smoothing !)
; ps : if ps = 1 then we switch into device, decomposed=0 ... resolve problems when ploting on a ps file
; overplot : nature of the overplot (default overplot=2),
	; - overplot=1 then we show only lines for the input frequency table. Table syntax: [0:lmax, 0:Nmax]
	; - overplot=2 then show lines + diamonds for the input frequency table. Table syntax: [0:lmax, 0:Nmax]
	; - overplot=3 then show only diamonds for the input frequency table. Table syntax: [0:lmax, 0:Nmax]
	; - overpolot=4 then use whisker boxes ... frequency_table MUST BE a table COMPATIBLE WITH THE FUNCTION : Draw_BoxAndWiskers_horizontal
	; - overplot=5 then we overplot 2 echelle diagram. For example, one containing true data and one containing a model of the data (variable spec2)
; color : color code to be used (default color='red')... used by the function fsc_color.

;ADDED 01/04/2014 : freq_scd2... usefull to overplot two tables : eg. model and obs
pro echellecorot_v3,freq,spec_reg,delta,nbegin,uu,freq_table=freq_t,freq_scd2=freq_t2, smooth_coef=smooth_coef,ps=ps,overplot=overplot,$
	color=color, scd_color=scd_color, show_ylabel=show_ylabel, charsize=charsize,$
	n_number=n_number,print=print, doublecolor=doublecolor, doublepos=doublepos, bar=bar, $
	;tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, $
	tag_label=tag_label, tag_pos=tag_pos, tag_size=tag_size, tag_color=tag_color, set_plot=set_plot,$
	extralines=extralines, pextralines_col=pextralines_col, style_extralines=style_extralines, recenter=recenter ;, rota=rota, congridon=congridon


A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A);, /FILL
sym_perso=dblarr(2, n_elements(A)) ; for the legends
sym_perso[0,*]=cos(A) ; for the legends
sym_perso[1,*]=sin(A) ; for the legends

printer='' ; relicat of the appourchaux function
index=''
;bar=1

if n_elements(set_plot) eq 0 then set_plot='X' ; Default, we plot on a linux system
if n_elements(recenter) eq 0 then recenter=0
if n_elements(n_number) eq 0 then n_number= 21
if n_elements(smooth_coef) eq 0 then smooth_coef=1
;if n_elements(tag_label) eq 0 then tag_label=''
if n_elements(tag_pos) eq 0 then tag_pos=[5.1, 90]
if n_elements(tag_size) eq 0 then tag_size=1.
;if n_elements(rota) eq 0 then rota=0
;if n_elements(congridon) eq 0 then congridon=0

if n_elements(show_ylabel) eq 0 then show_ylabel=1
if n_elements(charsize) eq 0 then charsize=2.
a=smooth(spec_reg,smooth_coef)


if n_elements(extralines) eq 0 then begin
	extralines=-1
	extralines_color=-1
endif
if n_elements(extralines) ne n_elements(pextralines_col) then pextralines_col=strarr(n_elements(extralines)) + 'Black'


n_numb=n_number
if nbegin lt 0 then nbegin=0
nbeg=nbegin

if n_elements(freq_t) ne 0  then if  freq_t[0] ne -1 then begin
if overplot ne 4 AND overplot ne 5 then freq_table=freq_t
if overplot eq 4 then begin ; if overplot = 4 we assume that the table is COMPATIBLE with Draw_BoxAndWiskers_horizontal
	stat_synthese_freq=freq_t
	freq_table=dblarr(n_elements(stat_synthese_freq[0,*,0]),n_elements(stat_synthese_freq[0,0,*]))
	freq_table[*,*]=stat_synthese_freq[3,*,*] ; freq_table contain ONLY THE MEDIAN POSITION
endif
endif

if n_elements(freq_t2) ne 0  then if  freq_t2[0] ne -1 then begin
if overplot ne 4 AND overplot ne 5 then freq_table2=freq_t2
if overplot eq 4 then begin ; if overplot = 4 we assume that the table is COMPATIBLE with Draw_BoxAndWiskers_horizontal
	stat_synthese_freq2=freq_t2
	freq_table2=dblarr(n_elements(stat_synthese_freq2[0,*,0]),n_elements(stat_synthese_freq2[0,0,*]))
	freq_table2[*,*]=stat_synthese_freq2[3,*,*] ; freq_table contain ONLY THE MEDIAN POSITION
endif
endif

if n_elements(overplot) eq 0 then overplot=-1

	nn=N_params()
	s={echelleplot,windowi:0.,delta:0.,echelletitle:''}
	if (nn LT N_tags(s)) then begin
		nn=nn-1
		read_structure,s,nn
		names=tag_names(s)
		for i=nn,N_tags(s)-1 do begin
			zzz=names(i)+'=s.'+names(i)
			;print,zzz
			r=execute(names(i)+'=s.'+names(i))
		endfor
	endif

	if (printer ne '') then begin
		set_plot,'ps',/interpolate
		file='corotvg'+index+'.ps'
		device,/landscape,/color,bits=8,filename=file
	endif
	if (n_elements(ps) ne 0) then begin
		set_plot,'ps',/interpolate
		device,/times
	endif

	del=1d*delta

	windowi=del

	window=windowi/2.

	Na=N_elements(a)
	resol=1d*(freq[10]-freq[9])

	;print,'Na=',Na
	;print,'resolution',resol
	Nz=1d*del/resol+2
	b=fltarr(Nz,n_numb)
	if overplot eq 5 then b2=fltarr(Nz,n_numb)
	middlefreq=1d*del/2d/resol

	nmax=floor(max(freq)/del)-1
	;print,"nmax=",nmax
;	nn=min([nmax,nbeg+20])
	nn=min([nmax,nbeg+n_numb-1])

	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
		nu=dblarr(n_elements(freq_table[*,0]),n_elements(freq_table[0,*])*3)
		freq_y=dblarr(n_elements(freq_table[*,0]),n_elements(freq_table[0,*])*3)
	endif
	if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
		nu2=dblarr(n_elements(freq_table2[*,0]),n_elements(freq_table2[0,*])*3)
		freq_y2=dblarr(n_elements(freq_table2[*,0]),n_elements(freq_table2[0,*])*3)
	endif

	if nbeg lt 0 then begin
		nbeg=0
		nn=10
	endif
	iend_max=round(1d*(nn+1)*del/resol)
	if iend_max gt n_elements(a) then nn=fix(n_elements(a)*resol/del)-1 ; if we are at upper edge of the spectrum!
	corr_t=0d
	corr_t2=0d
	for i=nbeg,nn do begin
		ibegin=round(1d*i*del/resol)
		iend=round(1d*(i+1)*del/resol)
		b(0:(iend-ibegin),i-nbeg)=a(ibegin:iend)
		if overplot eq 5 then b2(0:(iend-ibegin),i-nbeg)=a2(ibegin:iend)

		;print,i,min(b(*,i-nbeg)),max(b(*,i-nbeg))
		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
			f1=freq[floor(ibegin)] & f2=freq[ceil(ibegin)]
			f= (f2 - f1) * (ibegin - floor(ibegin)) + f1
			corr_t=corr_t+ f/(resol*ibegin)
		endif
		if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
			f1=freq[floor(ibegin)] & f2=freq[ceil(ibegin)]
			f= (f2 - f1) * (ibegin - floor(ibegin)) + f1
			corr_t2=corr_t2+ f/(resol*ibegin)
		endif
	endfor

	if n_elements(freq_table) ne 0 then begin

		lmax=n_elements(freq_table[*,0])-1

	corr_t=corr_t/(nn-nbeg +1)
	kk=0d & phase=0d
	for i=nbeg,nn do begin
		loop=0d
		ibegin=1d*i*del/resol
		iend=1d*(i+1)*del/resol ;-1

		l_tmpQ=dblarr(lmax+1)
		tmpQ=dblarr(lmax+1, 20)-1
		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
			for l=0, lmax do begin
				ta=where(freq_table[l,*] ge 1d*freq[ibegin] AND freq_table[l,*] le 1d*freq[iend])
				if ta[0] ne -1 then begin
					tmpQ[l, 0:n_elements(ta)-1]=ta
					l_tmpQ[l]=n_elements(ta)
				endif
			endfor
			dup_max=max(l_tmpQ)
			for j=0, dup_max-1 do begin
				for l=0, lmax do begin
					if tmpQ[l,j] ne -1 then begin
						nu[l, kk]=freq_table[l, tmpQ[l,j]]- corr_t*resol*ibegin
						freq_y[l, kk]=freq_table[l, tmpQ[l,j]]
					endif
				endfor
				kk=kk+1d
			endfor
		endif
	endfor

	endif
; --------- table 2 ---------
if n_elements(freq_table2) ne 0 then begin

	lmax=n_elements(freq_table2[*,0])-1

	corr_t2=corr_t2/(nn-nbeg +1)
	kk=0d & phase=0d
	for i=nbeg,nn do begin
		loop=0d
		ibegin=1d*i*del/resol
		iend=1d*(i+1)*del/resol ;-1

		l_tmpQ2=dblarr(lmax+1)
		tmpQ2=dblarr(lmax+1, 20)-1
		if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
			for l=0, lmax do begin
				ta=where(freq_table2[l,*] ge 1d*freq[ibegin] AND freq_table2[l,*] le 1d*freq[iend])
				if ta[0] ne -1 then begin
					tmpQ2[l, 0:n_elements(ta)-1]=ta
					l_tmpQ2[l]=n_elements(ta)
				endif
			endfor
			dup_max=max(l_tmpQ2)
			for j=0, dup_max-1 do begin
				for l=0, lmax do begin
					if tmpQ2[l,j] ne -1 then begin
						nu2[l, kk]=freq_table2[l, tmpQ2[l,j]]- corr_t2*resol*ibegin
						freq_y2[l, kk]=freq_table2[l, tmpQ2[l,j]]
					endif
				endfor
				kk=kk+1d
			endfor
		endif
	endfor

	endif

;	help,b
	;if congridon eq 0 then zz=b else zz=congrid(b, Nz, Nz)
	zz=b
	if recenter eq 1 then uu=zz(middlefreq-window/resol:middlefreq+window/resol,*)
	if recenter eq 0 then begin
		uu=zz ;(middlefreq-window/resol:middlefreq+window/resol,*)
		;middlefreq=0
	endif

	if overplot eq 5 AND recenter eq 1 then begin
		zz2=congrid(b2,Nz,100)
		uu2=zz2(middlefreq-window/resol:middlefreq+window/resol,*)
	endif
	if overplot eq 5 AND recenter eq 0 then begin
		zz2=congrid(b2,Nz,100)
		uu2=zz2 ;(middlefreq-window/resol:middlefreq+window/resol,*)
	endif

	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin

		if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then nu_scale=nu-middlefreq*resol*(1d +(corr_t-1)*2)

		pos_dummy=where(nu_scale eq -middlefreq*resol*(1d +(corr_t-1)*2)) ; these positions had 0... we need to let the 0
		nu_scale[pos_dummy]=0d

		if overplot eq 4 then begin
			param_stat=dblarr(n_elements(nu_scale[0,*]),n_elements(nu_scale[*,0]),n_elements(stat_synthese_freq[*,0,0]))

			for i=0,n_elements(nu_scale[0,*])-1 do begin
				for j=0,n_elements(stat_synthese_freq[*,0,0])-1 do begin
					for l=0,n_elements(nu_scale[*,0])-1 do $
						param_stat[i,l,j]=nu_scale[l,i]+(stat_synthese_freq[j,l,i]-stat_synthese_freq[3,l,i]) ; using nu_scale, we compute the error boxes
				endfor
			endfor
		endif

		if n_elements(overplot) eq 0 then overplot=2 ; definee the defaut value of overplot : we plot lines + diamonds !
		if n_elements(color) eq 0 then color=['Red','Red','Red','Red','Red']
		if n_elements(color) eq 1 then color=[color,color,color] ; 1 color code by degree
		if n_elements(scd_color) eq 0 then color_scd=color
		if n_elements(scd_color) eq 1 then color_scd=[scd_color, scd_color,scd_color]
		if n_elements(scd_color) gt 1 then color_scd=scd_color

		x0=-delta/2
		no_zero=where(freq_y ne 0)

		if n_elements(doublepos) ne 0 then $
			posy_dble=where(freq_y eq doublepos) ; find the position of the specially tagged color !!

		freq_y[no_zero]=freq_y[no_zero] - (nu_scale[no_zero]- x0) ;+ window ; the last term is to center on the y axis the frequency plot

	endif

; ----- Table 2 ------
	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, COS(A), SIN(A) ;, /FILL
	if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin

		if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then nu_scale2=nu2-middlefreq*resol*(1d +(corr_t2-1)*2)

		pos_dummy2=where(nu_scale2 eq -middlefreq*resol*(1d +(corr_t2-1)*2)) ; these positions had 0... we need to let the 0
		nu_scale2[pos_dummy2]=0d

		if overplot eq 4 then begin
			param_stat2=dblarr(n_elements(nu_scale2[0,*]),n_elements(nu_scale2[*,0]),n_elements(stat_synthese_freq2[*,0,0]))

			for i=0,n_elements(nu_scale2[0,*])-1 do begin
				for j=0,n_elements(stat_synthese_freq2[*,0,0])-1 do begin
					for l=0,n_elements(nu_scale2[*,0])-1 do $
						param_stat2[i,l,j]=nu_scale2[l,i]+(stat_synthese_freq2[j,l,i]-stat_synthese_freq2[3,l,i]) ; using nu_scale, we compute the error boxes
				endfor
			endfor
		endif

		if n_elements(overplot) eq 0 then overplot=2 ; definee the defaut value of overplot : we plot lines + diamonds !
		if n_elements(color) eq 0 then color=['Red','Red','Red','Red','Red']
		if n_elements(color) eq 1 then color=[color,color,color] ; 1 color code by degree

		x0=-delta/2
		no_zero=where(freq_y2 ne 0)

		if n_elements(doublepos2) ne 0 then $
			posy_dble=where(freq_y2 eq doublepos2) ; find the position of the specially tagged color !!

		freq_y2[no_zero]=freq_y2[no_zero] - (nu_scale2[no_zero]- x0) ;+ window ; the last term is to center on the y axis the frequency plot

	endif

	help,uu
	freq_u=TeXtoIDL('\mu')+'Hz'
	print,max(uu)


	fbegin=freq[round(1d*(nbeg)*del/resol)] ;+delta/2
	fend=freq[round(1d*nn*del/resol)] ;+ delta/2
	print, 'Frequency range:', fbegin, fend

if n_elements(ps) eq 0 then begin
	set_plot,set_plot
	device, decomposed=0
endif
if n_elements(ps) eq 1 then if ps eq 0 then begin
	set_plot,set_plot
	device, decomposed=0
endif
!p.font=0 ;& device,set_font='times' ; police de caractere times
;loadct,39  ; met une table de couleur avec 255 = Noir et 1 = Blanc
loadct,39
if n_elements(print) eq 0 then print=0
if print eq 1 then loadct, 0 ; table with a white background

!EDIT_INPUT=50

uu0=uu

xlabel='!3Frequency ('+freq_u+')'
if show_ylabel eq 1 then begin
	ylabel='!3Frequency ('+freq_u+')'
endif else begin
	ylabel=''
endelse

;if recenter eq 0 then begin
;	bmin=0
;	bmax=delta
;endif
if recenter eq 1 OR recenter eq 0 then begin
	bmin=-window
	bmax=window
endif

if n_elements(bar) eq 0 then bar=1
	if print eq 0 then begin
		if bar eq 1 then tvframe,uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel,/bar else $ ;,/center
		tvframe,uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel;,/center
	endif
	if print ne 0 then begin
		if bar eq 1 then tvframe,-uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel,/bar else $ ;,/center
		tvframe,-uu,xrange=[bmin,bmax],yrange=[fbegin,fend],charsize=charsize,$
			ytitle=ylabel,xtitle=xlabel ;,/center
	endif
	if extralines[0] ne -1 then begin
		for i=0, n_elements(extralines)-1 do $
			if extralines[i]+ bmin ge fbegin AND extralines[i]+bmax le fend then begin
;				norm0=sqrt( ( -window + window)^2 + (extralines[i] + window - nu_center)^2 )
;				norm1=sqrt( ( window + window)^2 + (extralines[i] - window - nu_center)^2 )
;				x0= -window + norm0 * cos( rota*!pi/180d + atan( (extralines[i] + window - nu_center)/(-window + window)) )
;				x1= -window + norm1 * cos( rota*!pi/180d + atan( (extralines[i] - window - nu_center)/( window + window)) )
;				y0=nu_center + norm0 * sin( rota*!pi/180d + atan( (extralines[i] + window - nu_center)/(-window + window)) )
;				y1=nu_center + norm1 * sin( rota*!pi/180d + atan( (extralines[i] - window - nu_center)/( window + window)) )
;				plots, [x0 , x1], [y0, y1], color=fsc_color(pextralines_col[i]), linestyle=style_extralines[i]
				plots, [bmin , bmax], [extralines[i]+ bmax, extralines[i]+bmin], color=fsc_color(pextralines_col[i]), linestyle=style_extralines[i]
				;plots, [-window , window], [extralines[i], extralines[i]], color=fsc_color(pextralines_col[i]), linestyle=style_extralines[i]
			endif
	endif


	if n_elements(freq_table) ne 0  then if  freq_t[0] ne -1 then begin
		if n_elements(nu_scale[*,0]) ge 6 then begin

			maxi=3
			test=where(nu_scale[5,*] ne 0) ; the last coloumn contains extra_dots
			if test[0] ne -1 then $
				oplot, nu_scale[5,where(nu_scale[5,*] ne 0)],freq_y[5,where(nu_scale[5,*] ne 0)],color=fsc_color(color[5]),psym=4,thick=3.,symsize=1.5

		endif else maxi=n_elements(nu_scale[*,0])-1

		if overplot eq 1 then for l=0,maxi do begin
				;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color)
				test=where(nu_scale[l,*] ne 0)
				if test[0] ne -1 then $
					oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
				endfor

		if overplot eq 2 then for l=0,maxi do begin
					test=where(nu_scale[l,*] ne 0)
					if test[0] ne -1 then begin
						tab=reform(nu_scale[l,where(nu_scale[l,*] ne 0)])
						freq_tab=reform(freq_y[l,where(nu_scale[l,*] ne 0)])
					endif
					if test[0] ne -1 then begin
						secu=0d & jj=0d & secu2=0d & ii=0d
						while (ii lt n_elements(tab)-1) AND secu lt 10 do begin
							d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
							d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
							secu2=0d
							while (d1_p le d2_p) AND (d1_m le d2_m) AND ii lt n_elements(tab)-2 AND secu2 lt 100 do begin
								;if ii eq 7 then stop
								plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(color[l])
								ii=ii+1
								d1_p=abs(tab[ii+1] -tab[ii]) & d2_p=abs(tab[ii+1]+2*window -tab[ii])
								d1_m=abs(tab[ii+1] -tab[ii]) & d2_m=abs(tab[ii+1]-2*window -tab[ii])
								secu2=secu2+1
							endwhile
							secu=secu+1
							ii=ii+1
						endwhile
						ii=ii-1
						if d1_p le d2_p then if n_elements(tab) ge 2 then plots, [tab[ii], tab[ii+1]], [freq_tab[ii], freq_tab[ii+1]], color=fsc_color(color[l])
						;oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
						oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l]),psym=4,thick=3.,symsize=1.5
					endif
				endfor
		if overplot eq 3 then for l=0,maxi do begin
					;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color),psym=4,symsize=1.5,thick=3
					test=where(nu_scale[l,*] ne 0)
					if test[0] ne -1 then $
						oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l]),psym=4,thick=3.,symsize=1.5
				endfor
		if overplot eq 4 then for l=0,n_elements(nu_scale[*,0])-1 do begin

					if l eq 0 then begin
						c='Orange' & t=0
					endif
					if l eq 1 then begin
						c='Yellow' & t=0
					endif
					if l eq 2 then begin
						c='red' & t=0
					endif
					if l eq 3 then begin
						c='brown' & t=0
					endif
					x_axis=findgen(nn)+nbeg+0.5
					;width=delta/4
					width=1d/4
					for i=0,n_elements(param_stat[*,0,0])-1 do $
						;Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,c], WIDTH=width, XLOC=stat_synthese_freq[3,l,i],1,t
						Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,'white'], WIDTH=width, XLOC=x_axis[i],1,t
					endfor
	endif

; ------- Table 2 -------
	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, COS(A), SIN(A) ;, /FILL
	if n_elements(freq_table2) ne 0  then if  freq_t2[0] ne -1 then begin
		if n_elements(nu_scale2[*,0]) ge 6 then begin

			maxi=3
			test=where(nu_scale2[5,*] ne 0) ; the last coloumn contains extra_dots
			if test[0] ne -1 then $
				oplot, nu_scale2[5,where(nu_scale2[5,*] ne 0)],freq_y2[5,where(nu_scale2[5,*] ne 0)],color=fsc_color(color_scd[5]),psym=8,thick=5.,symsize=1.5

		endif else maxi=n_elements(nu_scale2[*,0])-1

		if overplot eq 1 then for l=0,maxi do begin
				;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color)
				test=where(nu_scale2[l,*] ne 0)
				if test[0] ne -1 then $
					oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(color_scd[l]),psym=8
				endfor

		if overplot eq 2 then for l=0,maxi do begin
					test=where(nu_scale2[l,*] ne 0)
					if test[0] ne -1 then begin
						tab2=reform(nu_scale2[l,where(nu_scale2[l,*] ne 0)])
						freq_tab2=reform(freq_y2[l,where(nu_scale2[l,*] ne 0)])
					endif
					if test[0] ne -1 then begin
						secu=0d & jj=0d & secu2=0d & ii=0d
						while (ii lt n_elements(tab)-1) AND secu lt 10 do begin
							d1_p=abs(tab2[ii+1] -tab2[ii]) & d2_p=abs(tab2[ii+1]+2*window -tab2[ii])
							d1_m=abs(tab2[ii+1] -tab2[ii]) & d2_m=abs(tab2[ii+1]-2*window -tab2[ii])
							secu2=0d
							while (d1_p le d2_p) AND (d1_m le d2_m) AND ii lt n_elements(tab2)-2 AND secu2 lt 100 do begin
								;if ii eq 7 then stop
								plots, [tab2[ii], tab2[ii+1]], [freq_tab2[ii], freq_tab2[ii+1]], color=fsc_color(color_scd[l])
								ii=ii+1
								d1_p=abs(tab2[ii+1] -tab2[ii]) & d2_p=abs(tab2[ii+1]+2*window -tab2[ii])
								d1_m=abs(tab2[ii+1] -tab2[ii]) & d2_m=abs(tab2[ii+1]-2*window -tab2[ii])
								secu2=secu2+1
							endwhile
							secu=secu+1
							ii=ii+1
						endwhile
						ii=ii-1
						if d1_p le d2_p then if n_elements(tab2) ge 2 then plots, [tab2[ii], tab2[ii+1]], [freq_tab2[ii], freq_tab2[ii+1]], color=fsc_color(color_scd[l])
						;oplot, nu_scale[l,where(nu_scale[l,*] ne 0)],freq_y[l,where(nu_scale[l,*] ne 0)],color=fsc_color(color[l])
						oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(color_scd[l]),psym=8,thick=3.,symsize=1.5
					endif
				endfor
		if overplot eq 3 then for l=0,maxi do begin
					;oplot, nu_scale[l,*],findgen(nn)+nbeg+0.5,color=fsc_color(color),psym=8,symsize=1.5,thick=3
					test=where(nu_scale2[l,*] ne 0)
					if test[0] ne -1 then $
						oplot, nu_scale2[l,where(nu_scale2[l,*] ne 0)],freq_y2[l,where(nu_scale2[l,*] ne 0)],color=fsc_color(color_scd[l]),psym=8,thick=3.,symsize=1.5
				endfor
		if overplot eq 4 then for l=0,n_elements(nu_scale2[*,0])-1 do begin

					if l eq 0 then begin
						c='Orange' & t=0
					endif
					if l eq 1 then begin
						c='Yellow' & t=0
					endif
					if l eq 2 then begin
						c='red' & t=0
					endif
					if l eq 3 then begin
						c='brown' & t=0
					endif
					x_axis=findgen(nn)+nbeg+0.5
					;width=delta/4
					width=1d/4
					for i=0,n_elements(param_stat[*,0,0])-1 do $
						;Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,c], WIDTH=width, XLOC=stat_synthese_freq[3,l,i],1,t
						Draw_BoxAndWiskers_horizontal, param_stat[i,l,*], COLOR=[c,'white'], WIDTH=width, XLOC=x_axis[i],1,t
					endfor
	endif

	if overplot eq 5 then stop

	if n_elements(doublepos) ne 0 then begin
		A = FINDGEN(17) * (!PI*2/16.)
		USERSYM, COS(A), SIN(A), /FILL
		sym_perso=dblarr(2, n_elements(A)) ; for the legends
		sym_perso[0,*]=cos(A) ; for the legends
		sym_perso[1,*]=sin(A) ; for the legends
		if posy_dble[0] ne -1 then $
			plots, nu_scale[posy_dble], freq_y[posy_dble], color=fsc_color(doublecolor), psym=8, symsize=2; else $
				;stop
	endif

if n_elements(tag_label) ne 0 then xyouts, tag_pos[0]*window/100 - window, tag_pos[1]*fend/100, tag_label, $
										charsize=tag_size, color=fsc_color(tag_color)


;	if n_elements(ps) eq 0 then device,/close
;	if n_elements(ps) eq 1 then if ps eq 0 then device,/close
end
