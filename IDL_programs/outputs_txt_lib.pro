; Ensemble of small functions that allows us to read ASCII outputs

; Iteratively remove from the string values until the chain is blank
function str_to_arrdbl, str
	
	str2=strtrim(str,2);
	cpt=0.;
	while str2 ne '' do begin
		b_str=byte(str2)
		pos_blank=min(where(b_str eq 32));
		if pos_blank[0] ne -1 then begin
			if cpt eq 0 then vals=double(string(b_str[0: pos_blank-1]))
			if cpt gt 0 then vals=[vals, double(string(b_str[0: pos_blank-1]))]
			str2=strtrim(b_str[pos_blank+1:*],2);
		endif else begin ; last remaining value
			if cpt eq 0 then vals=double(string(b_str))
			if cpt gt 0 then vals=[vals, double(string(b_str))]
			str2=''
		endelse
		cpt=cpt+1
	endwhile
	
	return, vals
end

; Iteratively remove from the string values until the chain is blank
; THE SEPARATOR IS ASSUMED TO BE WHITE SPACE
function str_to_arrstr, str
	
	separator=' '
	bsep=byte(separator)
	str2=strtrim(str,2);
	cpt=0.;
	while str2 ne '' do begin
		b_str=byte(str2)
		pos_blank=min(where(b_str eq bsep[0]));
		if pos_blank[0] ne -1 then begin
			if cpt eq 0 then vals=string(b_str[0: pos_blank-1]);
			if cpt gt 0 then vals=[vals, string(b_str[0: pos_blank-1])]
			str2=strtrim(b_str[pos_blank+1:*],2);
		endif else begin
			if cpt eq 0 then vals=string(b_str)
			if cpt gt 0 then vals=[vals, string(b_str)]
			str2=''
		endelse
		cpt=cpt+1
	endwhile
	return, vals
end

; N: number of columns
; K: maximum number of lines
; This function gets the session of a file, not the file itself
function get_Ncolumns, file_session,N, K

	param=dblarr(K,N)
	a=''
	i=0d
      while EOF(file_session) ne 1 do begin
        readf, file_session, a ; read data
        uu=strsplit(a)
        N_uu=N_elements(uu)-1
        for j=0,N_uu-1 do begin
          	param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
        endfor
		param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		i=i+1
      endwhile

param0=param

param=param[0:i-1,*]

return,param
end

; Ncols: Exact number of columns
; Nrows: Exact number of lines to read
; This function gets the session of a file, not the file itself
function get_Ncolumns_Nrows, file_session, Ncols, Nrows

	param=dblarr(Nrows,Ncols)
	a=''
	i=0d
	
    while i lt Nrows do begin
        readf, file_session, a ; read data
        uu=strsplit(a)
        N_uu=N_elements(uu)-1
        for j=0,N_uu-1 do begin
          	param(i,j)=float(strmid(a,uu(j),uu(j+1)-uu(j)-1))
        endfor
		param(i, N_uu)= float(strmid(a,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 10))
		i=i+1
    endwhile
	;stop
return,param
end

