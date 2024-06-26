; A procedure that plots the spectrum for MS_Global models with a local zoom
; and 2 level of smoothness. If prior windows are available (aa,bb,cc,dd),
; those are as well plotted
pro MS_Global_fitplot, freq_show, spec_show, x, model, Dnu_moy, file_ps, $
  fmin, fmax, aa = aa, bb = bb, cc = cc, dd = dd
  compile_opt idl2

  resol = freq_show[1] - freq_show[0]

  mini = fmin
  maxi = fmax

  if Dnu_moy ge 55 then begin
    scoef = 0.01 * Dnu_moy / resol ; Defines the level of smooth on the shown power spectrum
  endif else begin
    if resol gt 1 then scoef = 0.005 * Dnu_moy / resol
    if resol le 1 then scoef = 0.05 * Dnu_moy / resol
  endelse
  if Dnu_moy le 15 then begin
    scoef = 2
  endif

  scoef2 = scoef * 2
  s_show = smooth(spec_show, scoef, /edge_truncate)
  s_show2 = smooth(spec_show, scoef2, /edge_truncate)

  ymax = max(s_show) * 1.2
  cte = min(s_show[where(freq_show ge fmin and freq_show le fmax)])

  ; e=write_on_ps_on(file_ps)
  nimp, name = file_ps, /paper, /eps
  plot, freq_show, /nodata, xr = [mini, maxi], yr = [0, ymax], /xst, /yst, background = fsc_color('white'), color = FSC_Color('Black'), thick = 1, $
    xtitle = 'frequency (microHz)', ytitle = 'power (ppm^2/microHz)'

  oplot, freq_show, s_show, color = FSC_Color('Grey')
  oplot, freq_show, s_show2, color = FSC_Color('Dark Grey')

  if n_elements(aa) gt 1 then oplot, freq_show, aa * 100 + cte, color = Fsc_color('Red'), thick = 4., linestyle = 1
  if n_elements(bb) gt 1 then oplot, freq_show, bb * 100 + cte, color = fsc_color('blue'), thick = 4., linestyle = 1
  if n_elements(cc) gt 1 then oplot, freq_show, cc * 50 + cte, color = fsc_color('Brown'), thick = 4., linestyle = 1
  if n_elements(dd) gt 1 then oplot, freq_show, dd * 50 + cte, color = fsc_color('Magenta'), thick = 4., linestyle = 1
  oplot, x, model[0, *], color = FSC_Color('Blue'), thick = 4
  if n_elements(model[*, 0]) gt 1 then $
    for i = 1, n_elements(model[1 : *, 0]) - 1 do oplot, x, model[i, *], color = FSC_Color('blue'), thick = 4, linestyle = 2

  ymax = max(s_show / model[0, *]) * 1.2
  plot, freq_show, /nodata, xr = [mini, maxi], yr = [0, ymax], /xst, /yst, background = fsc_color('white'), color = FSC_Color('Black'), thick = 1, $
    xtitle = 'frequency (microHz)', ytitle = 'residual power (ppm^2/microHz)'
  oplot, freq_show, s_show / model[0, *], color = FSC_Color('Grey')
  oplot, freq_show, s_show2 / model[0, *], color = FSC_Color('Dark Grey')

  indice = 0
  maxi = mini

  scoef = 0.01 * Dnu_moy / resol ; Defines the level of smooth on the shown power spectrum
  s_show = smooth(spec_show, scoef, /edge_truncate)
  ymax = max(s_show) * 1.2
  cte = min(s_show[where(freq_show ge fmin and freq_show le fmax)])

  while maxi lt fmax do begin
    ymax = max(s_show) * 1.2
    mini = fmin + indice * Dnu_moy - 0.05 * Dnu_moy
    maxi = fmin + (indice + 1) * Dnu_moy + 0.05 * Dnu_moy
    plot, freq_show, /nodata, xr = [mini, maxi], yr = [0, ymax], /xst, /yst, background = fsc_color('white'), color = FSC_Color('Black'), thick = 1, $
      xtitle = 'frequency (microHz)', ytitle = 'power (ppm^2/microHz)'

    oplot, freq_show, s_show, color = FSC_Color('Grey')
    oplot, freq_show, s_show2, color = FSC_Color('Dark Grey')
    if n_elements(aa) gt 1 then oplot, freq_show, aa * 100 + cte, color = Fsc_color('Red'), thick = 4., linestyle = 1
    if n_elements(bb) gt 1 then oplot, freq_show, bb * 100 + cte, color = fsc_color('blue'), thick = 4., linestyle = 1
    if n_elements(cc) gt 1 then oplot, freq_show, cc * 50 + cte, color = fsc_color('Brown'), thick = 4., linestyle = 1
    if n_elements(dd) gt 1 then oplot, freq_show, dd * 50 + cte, color = fsc_color('Magenta'), thick = 4., linestyle = 1
    oplot, x, model[0, *], color = FSC_Color('Blue'), thick = 4
    if n_elements(model[*, 0]) gt 1 then $
      for i = 1, n_elements(model[1 : *, 0]) - 1 do oplot, x, model[i, *], color = FSC_Color('Blue'), thick = 4, linestyle = 2

    ymax = max(s_show / model[0, *]) * 1.2
    plot, freq_show, /nodata, xr = [mini, maxi], yr = [0, ymax], /xst, /yst, background = fsc_color('white'), color = FSC_Color('Black'), thick = 1, $
      xtitle = 'frequency (microHz)', ytitle = 'residual power (ppm^2/microHz)'
    oplot, freq_show, s_show / model[0, *], color = FSC_Color('Grey')
    oplot, freq_show, s_show2 / model[0, *], color = FSC_Color('Dark Grey')

    indice = indice + 1.
  endwhile
  ; e=write_on_ps_off('')
  fimp
end