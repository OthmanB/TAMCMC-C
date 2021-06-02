function read_getmodel_file, model_file

	Nmax=100
	a=''
	c=''
	skip=''
	openr, 3, model_file
		readf, 3, a
		s=split_char(a, '=')
		modelname=s[1]
		readf, 3, a
		s=split_char(a, '=')
		plength=fix(split_line(s[1]))
		readf, 3, a
		s=split_char(a, '=')
		raw_params=double(split_line(s[1]))
		readf,3, c ; comment abour spectrum parameters
		specparams_comment=c
		readf, 3, a
		specparams=double(split_line(a)) ; 
		
		readf,3, skip ; not so usefull comment about file version
		readf,3, c ; comment
		modeparams_comment=c
		i=0
		b=-1
		while b ne 35 AND EOF(3) ne 1 do begin ; Unless we meet a comment section do...
			readf, 3, a
			b=byte(strtrim(a,2))
			b=b[0]
			if b ne 35 then begin
				line=double(split_line(a)) 
				if i eq 0 then modeparams=dblarr(Nmax, n_elements(line))
				modeparams[i, *]=line 
				i=i+1
			endif
		endwhile
		tmp=modeparams
		modeparams=modeparams[0:i-1, *]

		readf,3,c
		noiseparams_comment=c
		i=0
		b=-1
		while b ne 35 AND EOF(3) ne 1 do begin ; Unless we meet a comment section or an end of file do...
			readf, 3, a
			b=byte(strtrim(a,2))
			b=b[0]
			if b ne 35 then begin
				line=double(split_line(a)) 
				if i eq 0 then noiseparams=dblarr(Nmax, n_elements(line))
				noiseparams[i, *]=line 
				i=i+1
			endif
		endwhile
		tmp=noiseparams
		noiseparams=noiseparams[0:i-1,*]
	close,3

	struc={modelname:'', plength:intarr(n_elements(plength)), raw_params:dblarr(n_elements(raw_params)), $
		   modeparams:dblarr(n_elements(modeparams[*,0]), n_elements(modeparams[0,*])),  noiseparams:dblarr(n_elements(noiseparams[*,0]), n_elements(noiseparams[0,*])), $
		   modeparams_comment:'', noiseparams_comment:'', specparams_comment:''}

	struc.modelname=modelname
	struc.plength=plength
	struc.raw_params=raw_params
	struc.modeparams=modeparams
	struc.noiseparams=noiseparams
	struc.specparams_comment=specparams_comment
	struc.modeparams_comment=modeparams_comment
	struc.noiseparams_comment=noiseparams_comment
	;stop
	return, struc
end

function where_char, str, char
	b=byte(str)
	c=byte(char)
	pos=where(b eq c[0])
	return, pos
end

function split_char, str, char
	pos=where_char(str, char)
	s=byte(str)
	pass=0
	if pos[0] ne -1 then begin
		for i=0, n_elements(pos)-1 do begin
			if i eq 0 then begin
				line=strtrim(s[0:pos[0]-1],2)
			endif else begin
				line=[line, strtrim(s[pos[i-1]:pos[i]-1],2)]
			endelse
		endfor
		pass=1
	endif
	if pass eq 1 then line=[line, strtrim(s[pos[n_elements(pos)-1]+1:*],2)]
	if pass eq 0 then line =''
	return, line
end

function split_line, str
	uu=strsplit(str)
    N_uu=N_elements(uu)-1
    param=strarr(N_uu+1)
    for j=0,N_uu-1 do begin
    	param(j)=strmid(str,uu(j),uu(j+1)-uu(j)-1)
    endfor
	param(N_uu)= strmid(str,uu(N_uu),uu(N_uu)-uu(N_uu-1)+ 50)
 	return, param
end
