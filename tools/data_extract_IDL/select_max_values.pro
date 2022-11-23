; Cette fonction donne la positino du maximum d'un ensemble histogramme

;function select_max_values, histogram_values,x_hist_values,smooth_coeff
;
;nb_m=n_elements(histogram_values[*,0])
;ind_tab=dblarr(nb_m)
;max_pos=dblarr(nb_m)
;
;for i=0,nb_m-1 do begin
;	ind_tab[i]=where(smooth(histogram_values[i,*],smooth_coeff) eq max(smooth(histogram_values[i,*],smooth_coeff)))
;	max_pos[i]=x_hist_values[i,ind_tab[i]]
;endfor
;
;return, max_pos
;end


; Cette fonction donne la position du maximum pour un histogramme.
; Elle retourne aussi l'indice associé dans la table des x
function select_max_values, histogram_values,x_hist_values,smooth_coeff

	ind_tab=where(smooth(histogram_values,smooth_coeff) eq max(smooth(histogram_values,smooth_coeff)))
	max_pos=x_hist_values[ind_tab]

return, [ind_tab,max_pos]
end


