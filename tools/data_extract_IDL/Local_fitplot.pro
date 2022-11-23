; A procedure that plots the spectrum for MS_Global models with a local zoom
; and 2 level of smoothness. If prior windows are available (aa,bb,cc,dd),
; those are as well plotted
pro Local_fitplot, freq_show, spec_show, x, model, file_ps, $
	fmin, fmax, fmin_loc, fmax_loc, aa=aa,bb=bb,cc=cc,dd=dd

resol=freq_show[1]-freq_show[0]

mini=fmin & maxi=fmax

avg_freq= (fmin + fmax)/2 ; Proxy of numax

if avg_freq ge 800 then begin 
	scoef=0.001*avg_freq/resol ; Defines the level of smooth on the shown power spectrum
endif 
if avg_freq ge 300 and avg_freq lt 800 then begin
	scoef=0.0005*avg_freq/resol
endif
if avg_freq ge 100 and avg_freq lt 300 then begin
	scoef=0.0002*avg_freq/resol
endif
if avg_freq lt 100 then begin
	scoef=0.0001*avg_freq/resol
endif

scoef2=scoef*2
s_show=smooth(spec_show, scoef, /edge_truncate)
s_show2=smooth(spec_show, scoef2, /edge_truncate)

ymax=max(s_show)*1.2
cte=min(s_show[where(freq_show ge fmin AND freq_show le fmax)])

;e=write_on_ps_on(file_ps)
nimp,name=file_ps,/paper,/eps
plot, freq_show, /Nodata,xr=[mini,maxi],yr=[0,ymax],/xst,/yst,background=fsc_color('white'),color=FSC_Color('Black'),thick=1,$
  xtitle='frequency (microHz)',ytitle='power (ppm^2/microHz)'

oplot, freq_show,s_show,color=FSC_Color('Grey')
oplot, freq_show,s_show2,color=FSC_Color('Dark Grey')

if n_elements(aa) gt 1 then oplot, freq_show,aa*100+cte,color=Fsc_color('Red'),thick=4.,linestyle=1
if n_elements(bb) gt 1 then oplot, freq_show,bb*100+cte,color=fsc_color('blue'),thick=4.,linestyle=1
if n_elements(cc) gt 1 then oplot, freq_show,cc*50+cte,color=fsc_color('Brown'),thick=4.,linestyle=1
if n_elements(dd) gt 1 then oplot, freq_show,dd*50+cte,color=fsc_color('Magenta'),thick=4.,linestyle=1
oplot, x,model[0,*],color=FSC_Color('Blue'),thick=4
if n_elements(model[*,0]) gt 1 then $
	for i=1, n_elements(model[1:*,0])-1 do oplot, x,model[i,*],color=FSC_Color('blue'),thick=4, linestyle=2

ymax=max(s_show/model[0,*])*1.2
plot, freq_show, /Nodata,xr=[mini,maxi],yr=[0,ymax],/xst,/yst,background=fsc_color('white'),color=FSC_Color('Black'),thick=1,$
  xtitle='frequency (microHz)',ytitle='residual power (ppm^2/microHz)'
oplot, freq_show,s_show/model[0,*],color=FSC_Color('Grey')
oplot, freq_show,s_show2/model[0,*],color=FSC_Color('Dark Grey')


; -----------------------------
; ---------- LOCAL VIEW -------
ymax=max(s_show[where(freq_show ge fmin_loc AND freq_show le fmax_loc)])*1.2
cte=min(s_show[where(freq_show ge fmin_loc AND freq_show le fmax_loc)])

plot, freq_show, /Nodata,xr=[fmin_loc,fmax_loc],yr=[0,ymax],/xst,/yst,background=fsc_color('white'),color=FSC_Color('Black'),thick=1,$
  xtitle='frequency (microHz)',ytitle='power (ppm^2/microHz)'

oplot, freq_show,s_show,color=FSC_Color('Grey')
oplot, freq_show,s_show2,color=FSC_Color('Dark Grey')

if n_elements(aa) gt 1 then oplot, freq_show,aa*100+cte,color=Fsc_color('Red'),thick=4.,linestyle=1
if n_elements(bb) gt 1 then oplot, freq_show,bb*100+cte,color=fsc_color('blue'),thick=4.,linestyle=1
if n_elements(cc) gt 1 then oplot, freq_show,cc*50+cte,color=fsc_color('Brown'),thick=4.,linestyle=1
if n_elements(dd) gt 1 then oplot, freq_show,dd*50+cte,color=fsc_color('Magenta'),thick=4.,linestyle=1
oplot, x,model[0,*],color=FSC_Color('Blue'),thick=4
if n_elements(model[*,0]) gt 1 then $
	for i=1, n_elements(model[1:*,0])-1 do oplot, x,model[i,*],color=FSC_Color('blue'),thick=4, linestyle=2

ymax=max(s_show/model[0,*])*1.2
plot, freq_show, /Nodata,xr=[fmin_loc, fmax_loc],yr=[0,ymax],/xst,/yst,background=fsc_color('white'),color=FSC_Color('Black'),thick=1,$
  xtitle='frequency (microHz)',ytitle='residual power (ppm^2/microHz)'
oplot, freq_show,s_show/model[0,*],color=FSC_Color('Grey')
oplot, freq_show,s_show2/model[0,*],color=FSC_Color('Dark Grey')

fimp

end
