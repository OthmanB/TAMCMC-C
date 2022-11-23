function ascii2sav, dir_out, delete_txt_outputs

	variable_name=strarr(4)
	files=file_search(dir_out + '*.ASCII')
	if files[0] ne '' then begin
		for i=0, n_elements(files)-1 do begin
			struc=read_bin2txt_files(files[i])
			if i eq 0 then Nsamples=n_elements(struc.data)
			b=byte(files[i])
			posslash=max(where((b eq 47) OR (b eq 92))) ; detect last '/' OR '\'
			posdot=max(where(b eq 46)) ; detect last dot
			param=struc.data
			index=i
			variable_name=interpret_varnames_Cpp(struc.varname)
			if delete_txt_outputs eq 1 then print, '  - Saving in IDL format the file (the original ASCII will be deleted in the process): ' + files[i]
			if delete_txt_outputs eq 0 then print, '  - Saving in IDL format the file : ' + files[i]
			;stop
			save, index, param, variable_name, filename=dir_out + strtrim(b[posslash:posdot-1],2) + '.sav'
			if delete_txt_outputs eq 1 then spawn, 'rm ' + files[i]
		endfor
	endif else begin
		;if readbinfiles eq 0 then $
			print, 'No ASCII files found... Cannot proceed'
			stop
			;restore, f + '/000.sav' ; load the first file (should aleways exist of course!) ==> Gives Nsamples
			;Nsamples=n_elements(param)
	endelse

return, Nsamples
end
