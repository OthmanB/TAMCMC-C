@select_max_values
;;********** Procedure de DEBUG **********
;pro test_sigmafunction
;critere=95.; en %
;smooth_coeff=1.
;
;x=100.+findgen(10)
;moyenne=105. & sigma=2. & amplitude=5.
;y=gaussian_form(moyenne,sigma,x)
;;y=gaussian_curve(Amplitude,moyenne,sigma,x)
;plot, x,y,psym=2
;tt=estimate_1_sigma_error(x,y,critere,smooth_coeff)
;print,tt
;print, 'yo!'
;end
;;***************************************

;****** Version avec interpolation quadratique *********
;;tab_critere=s[10,25,50,75,90] ; defini les decile,quartile,mediane,etc...
;tab_critere=[2.25,16,50,84,97.75] ; defini les mediane +/- 2sigma,mediane, mediane +/- 1sigma
function estimate_1_sigma_error,x,y,critere,smooth_coeff,tab_critere

if variance(x) ne 0 then begin
y_smooth=smooth(y,smooth_coeff)
max_pos=select_max_values(y_smooth,x,1) ; on determine la position du max


Ptotal=total(y_smooth) ; on determine le nombre de coup total (calcul de la constante de normalisation)
Ptotal_p=total(y_smooth[max_pos[0]:*]) ; on determine le nombre de coup total à droite
Ptotal_m=total(y_smooth[0:max_pos[0]]) ; on determine le nombre de coup total à gauche
pas=x[1]-x[0] ; definition de la distance entre deux points du tableau des x

;****** Creation des tables d'effectif cumulé croissant par rapport au max *******
effectif_cumule_p=dblarr(n_elements(y)-max_pos[0])
effectif_cumule_tmp=0
for i=max_pos[0],n_elements(y)-1 do begin
	effectif_cumule_p[i-max_pos[0]]=effectif_cumule_tmp+100*y_smooth[i]/ptotal_p ; on defini une nouvelle valeur de l'effectif
	effectif_cumule_tmp=effectif_cumule_p[i-max_pos[0]]
endfor
if max_pos[0] ne 0 then begin
	effectif_cumule_m=dblarr(max_pos[0])
	effectif_cumule_tmp=0
	for i=max_pos[0],1,-1 do begin ; decrementation de 1 jusqu'à 0
		effectif_cumule_m[abs(i-max_pos[0])]=effectif_cumule_tmp+100*y_smooth[i]/ptotal_m ; on defini une nouvelle valeur de l'effectif
		effectif_cumule_tmp=effectif_cumule_m[abs(i-max_pos[0])]
	endfor
endif
;****** Creation d'une table d'effectif cumulé croissant par rapport au min ****
; l'objectif ici est d'estimer 10%,25%,50%,75%,90% de la pdf
effectif_cumule_t=dblarr(n_elements(y))
effectif_cumule_tmp=0
for i=long(0),n_elements(y)-1 do begin
	effectif_cumule_t[i]=effectif_cumule_tmp+100*y_smooth[i]/ptotal ; on defini une nouvelle valeur de l'effectif
	effectif_cumule_tmp=effectif_cumule_t[i]
endfor

;**************************************************************
;*********** Interpolation quadratique des tables *************
gain=2. ; facteur mutliplicatif associé au nombre de points apres interpolation
interpolated_effectif_p=interpol(effectif_cumule_p,gain*n_elements(effectif_cumule_p),/quadratic)
if max_pos[0] ne 0 then interpolated_effectif_m=interpol(effectif_cumule_m,gain*n_elements(effectif_cumule_m),/quadratic)
interpolated_effectif_t=interpol(effectif_cumule_t,gain*n_elements(effectif_cumule_t),/quadratic)
sigma_tab_p=findgen(n_elements(interpolated_effectif_p))*pas/gain ; a droite du max (relativement au max)
if max_pos[0] ne 0 then sigma_tab_m=findgen(n_elements(interpolated_effectif_m))*pas/gain ; a gauche du max (relativement au max)
sigma_tab_t=findgen(n_elements(interpolated_effectif_t))*pas/gain+min(x) ; a partir du min (en absolue)
;**************************************************************
;********** Recherche (par exces) de la valeur la plus proche du critere ******
c_p_ind=where(interpolated_effectif_p ge 1.*critere)
if max_pos[0] ne 0 then c_m_ind=where(interpolated_effectif_m ge 1.*critere) else c_m_ind=-2
if n_elements(c_p_ind) gt 1 then begin
c_p=interpolated_effectif_p[c_p_ind[1]] & sigma_plus=sigma_tab_p[c_p_ind[1]] ; on met 1 et pas 0 pour arondir à l'exces
endif else if c_p_ind[0] ne -1 then begin
c_p=interpolated_effectif_p[c_p_ind[0]] & sigma_plus=sigma_tab_p[c_p_ind[0]] ; on met 1 et pas 0 pour arondir à l'exces
	endif else begin; cas ou atteindre le critere (par ex 99.7%) est impossible à cause d'arondi du smooth
	c_p=max(interpolated_effectif_p) & sigma_plus=max(sigma_tab_p) ; on prend le max de l'erreur
	; ce max DOIT etre le dernier element du tableau en tout logique !
	print, 'Attention!! le critere n a pas été atteint. critere défini a : '+strtrim(c_p,1)+'%'
	endelse
if n_elements(c_m_ind) gt 1 then begin
c_m=interpolated_effectif_m[c_m_ind[1]] & sigma_moins=sigma_tab_m[c_m_ind[1]] ; on met 1 et pas 0 pour arondir à l'exces
endif else if c_m_ind[0] ge 0 then begin
c_m=interpolated_effectif_m[c_m_ind[0]] & sigma_moins=sigma_tab_m[c_m_ind[0]] ; on met 1 et pas 0 pour arondir à l'exces
endif else begin
	if c_m_ind[0] eq -1 then begin ; cas ou atteindre le critere (par ex 99.7%) est impossible à cause d'arondi du smooth
	c_m=max(interpolated_effectif_m) & sigma_moins=max(sigma_tab_m) ; on prend le max de l'erreur
	; ce max DOIT etre le dernier element du tableau en tout logique !
	print, 'Attention!! le critere n a pas été atteint. critere défini a : '+strtrim(c_m,1)+'%'
	endif
	if c_m_ind[0] eq -2 then begin
		c_m=0 & sigma_moins=min(x)
	endif
endelse


;--- Recherche par exces des proba à 10% 25% 50% 75% et 90% ---
c_t_ind_tmp=0.
;tab_critere=[10,25,50,75,90]
c_t=dblarr(n_elements(tab_critere)) & sigma_t=dblarr(n_elements(tab_critere))
for i=0,n_elements(tab_critere)-1 do begin
	c_t_ind_tmp=where(interpolated_effectif_t ge 1.*tab_critere[i])
	if n_elements(c_t_ind_tmp) gt 1 then begin ; le cas le plus courant
		c_t[i]=interpolated_effectif_t[c_t_ind_tmp[1]]
		sigma_t[i]=sigma_tab_t[c_t_ind_tmp[1]] ; on met 1 et pas 0 pour arondir à l'exces
	endif else if c_t_ind_tmp[0] ne -1 then begin
		c_t[i]=interpolated_effectif_t[c_t_ind_tmp[0]]
		sigma_t[i]=sigma_tab_t[c_t_ind_tmp[0]]
	endif else begin; cas ou atteindre le critere (par ex 99.7%) est impossible à cause d'arondi du smooth
		c_t[i]=max(interpolated_effectif_t)
		sigma_t[i]=max(sigma_tab_t) ; on prend le max de l'erreur
	; ce max DOIT etre le dernier element du tableau en tout logique !
	endelse
endfor

;******************************************************************************
; les elements 0 à 2 sont le max et l'erreur +/- autour à critere
; les elements au dela de 2 sont des critere du type quartile,decile,...
;if min(x) ge 0.9289 AND max(x) le 1.053 then begin
;	result=[0.601264,sigma_plus*8,sigma_moins*8,0.4,0.47552,0.5852917,0.7874171,0.9905,1.25054,1.652]
;endif else
result=[max_pos[1],sigma_plus,sigma_moins,min(sigma_tab_t),sigma_t,max(sigma_tab_t)] ; synthese des resultats si on a affaire à une variable

endif else result=[x[0],0,0,0,0,0,0,0,0,0] ; synthese des resultats si on a affaire à un parametre fixé
return, result
end

pro gauss_smooth, x, y, sigma

	N=20
	p=dblarr(N)
	x=p
	y=p
	r=p

	x=findgen(N)
    y[0]=0
    y[1]=0.5
    y[2]=0.1
    y[3]=0
    y[4]=0.3
    y[5:14]=1
    y[15:19]=0


   sigma=10

   for j=0, N-1 do begin
      x0=x[j];
      for i=0, N-1 do begin
         p[i]=y[j]*exp(-0.5* (x[i]-x0)^2/sigma^2)
      endfor
      r[j]=total(p)/sqrt(2.*!pi)/sigma
   endfor

plot, x, y, yr=[0,2]
oplot,x, r, color=fsc_color('Blue')
end
