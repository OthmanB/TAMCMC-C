; A procedure that plots the spectrum for MS_Global models with a local zoom
; and 2 level of smoothness. If prior windows are available (aa,bb,cc,dd),
; those are as well plotted
pro Gaussian_fitplot, freq_show, spec_show, x, model, file_ps, fmin, fmax, aa=aa

resol=freq_show[2]-freq_show[1]

mini=fmin & maxi=fmax

scoef=1./resol ; 1 microHz smoothing
scoef2=scoef*5
s_show=smooth(spec_show, scoef, /edge_truncate)
s_show2=smooth(spec_show, scoef2, /edge_truncate)

ymin=min(model)*0.1
ymax=max(s_show)*1.2
cte=min(s_show[where(freq_show ge fmin AND freq_show le fmax)])


nimp,name=file_ps,/paper,/eps
plot, freq_show, /Nodata,xr=[mini,maxi],yr=[ymin,ymax],/xst,/yst,background=fsc_color('white'),color=FSC_Color('Black'),thick=1,$
  xtitle='frequency (microHz)',ytitle='power (ppm^2/microHz)', /ylog, /xlog

oplot, freq_show,s_show,color=FSC_Color('Grey')
oplot, freq_show,s_show2,color=FSC_Color('Darsk Grey')

if n_elements(aa) gt 1 then oplot, freq_show,aa*100+cte,color=Fsc_color('Red'),thick=4.,linestyle=1
oplot, x,model[0,*],color=FSC_Color('Blue'),thick=4
if n_elements(model[*,0]) gt 1 then $
	for i=1, n_elements(model[1:*,0])-1 do oplot, x,model[i,*],color=FSC_Color('blue'),thick=4, linestyle=2

ymax=max(s_show/model[0,*])*1.2
plot, freq_show, /Nodata,xr=[mini,maxi],yr=[0,ymax],/xst,/yst,background=fsc_color('white'),color=FSC_Color('Black'),thick=1,$
  xtitle='frequency (microHz)',ytitle='residual power (ppm^2/microHz)'
oplot, freq_show,s_show/model[0,*],color=FSC_Color('Grey')
oplot, freq_show,s_show2/model[0,*],color=FSC_Color('Dark Grey')

indice=0
maxi=mini
fimp

end
