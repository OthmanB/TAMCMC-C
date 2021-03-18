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
				xtitle=xtitle, ytitle='!6Probability density', title=title, charsize=size_c, $
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

if normalize eq 1 then output_table[1,*]=output_table[1,*]/N

return, output_table
end