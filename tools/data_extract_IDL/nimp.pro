
pro nimp,name=name,paper=paper,portrait=portrait,black=black,eps=eps



IF keyword_set(paper) THEN BEGIN
!p.charsize=2.0
!p.charthick=5
!p.thick=5
!x.thick=4
!y.thick=4
;!P.font=6
ENDIF
!P.font=6

set_plot,'ps',/interpolate
if NOT keyword_set(name) THEN BEGIN
name=' '
  print,'Dime el nombre del fichero'
  read,name
ENDIF

IF not keyword_set(portrait) THEN BEGIN
   IF not keyword_set(eps) THEN BEGIN
     device,file=name,/landscape,xsize=25,ysize=18,/color, bits_per_pixel=8
    ENDIF ELSE BEGIN
     device,file=name,/landscape,xsize=25,ysize=18,/color, bits_per_pixel=8;,/encaps
    ENDELSE
ENDIF ELSE BEGIN
   IF not keyword_set(eps) THEN BEGIN
     device,file=name,/portrait,/color,xsize=19,ysize=26,encaps=0,yoffset=1,xoffset=.5
    ENDIF ELSE BEGIN
     device,file=name,/portrait,/color,xsize=19,ysize=26,yoffset=1,xoffset=.5 ;,/encaps

    ENDELSE

ENDELSE
if keyword_set(black) THEN !P.BACKGROUND=0


return
end
