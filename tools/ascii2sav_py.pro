; Function specifically designed to convert a temporarily created ascii file by prepare_lc_kepler.read_py
; into a idl sav file. This is why the function has no argument
pro ascii2sav_pyspectrum, file_in

file_in='tmp_in.ascii'
file_out='tmp_out.sav'

skip=1 ; defaut we skip one line only
ref_N=1 ; defaut we identify the zero non-used tab elements with column 1
K=1000000.
N=2
openr, 3, file_in
	param=dblarr(K,N)-1
	a=''
	i=0d
      while EOF(3) ne 1 do begin
      	if i lt skip then readf,3,format='(q)'
      	if i eq skip+1 then begin
      		readf,3, a
        	b=byte(strtrim(a,2)) ; remove any space before and after the text
        	id_number=strtrim(b[3:n_elements(b)-1],2) ; Remove the first symbol used as comment indicator
        endif
        if i ge skip+1 then begin
          	readf,3,a ; read data
          	uu=strsplit(a)
          	N_uu=N_elements(uu)-1
          	for j=0,N_uu-1 do begin
          		param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
          	endfor
			param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		endif
		i=i+1
      endwhile
close,3
param0=param
test=where(param[*,ref_N] ne -1)

param=param[test,*]

freq=dblarr(1, n_elements(param[*,0]))
spec_reg=dblarr(1, n_elements(param[*,0]))
freq[0,*]=param[*,0]
spec_reg[0,*]=param[*,1]
save, freq, spec_reg, id_number, filename=file_out

; Once the file was read properly and saved into a save file, the original file is erased
;file_delete, file_in

end


pro ascii2sav_pylightcurve, file_in

file_in='tmp_in.ascii'
file_out='tmp_out.sav'

skip=1 ; defaut we skip one line only
ref_N=1 ; defaut we identify the zero non-used tab elements with column 1
K=1000000.
N=2
openr, 3, file_in
  param=dblarr(K,N)-1
  a=''
  i=0d
      while EOF(3) ne 1 do begin
        if i lt skip then readf,3,format='(q)'
        if i eq skip+1 then begin
          readf,3, a
          b=byte(strtrim(a,2)) ; remove any space before and after the text
          id_number=strtrim(b[3:n_elements(b)-1],2) ; Remove the first symbol used as comment indicator
        endif
        if i ge skip+1 then begin
            readf,3,a ; read data
            uu=strsplit(a)
            N_uu=N_elements(uu)-1
            for j=0,N_uu-1 do begin
              param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
            endfor
      param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
    endif
    i=i+1
      endwhile
close,3
param0=param
test=where(param[*,ref_N] ne -1)

param=param[test,*]

freq=dblarr(1, n_elements(param[*,0]))
spec_reg=dblarr(1, n_elements(param[*,0]))
time=param[*,0]
flux=param[*,1]
save, time, flux, id_number, filename=file_out

; Once the file was read properly and saved into a save file, the original file is erased
;file_delete, file_in

end