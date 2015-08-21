pro create_cola_script

SPAWN, 'rm -r models'
SPAWN, 'mkdir models'

ellint_filename = 'ellint.txt'
;cube_name = '749237_cube_regrid'
cube_name = '749237_cube'
;nhi_name = '749237_nhi_regrid'
nhi_name = '749237_nhi'
beamsizes = '10 10'
vsys = '377'
center_pos = '* 12 26 23 * 27 44 44.50'

; set parameter ranges
R=10 * FINDGEN(5) + 10.

incl=[49, 59]
incl=indgen(11)*1 + 49
;incl=[90]
;pa=[50, 60, 70]

;Vrot[R] = Vflat*[1 - EXP(-R/lflat)].  
vflat=[20, 40]
vflat=indgen(11)*2 + 20
lflat=[10, 30]
lflat=indgen(11)*2 + 10

;1 arcsec ~ 53 pc, 1.5 arcsec/pix
; so 50 pc to 500 pc --> 1.5 to 6 arcsec
Z0=[1, 21]
Z0=indgen(21)*1 + 1

print, incl
print, vflat
print, lflat
print, Z0

pa=[240]
vdisp=[0]

PRINT, 'Number of models: ', N_ELEMENTS(incl)*N_ELEMENTS(pa)*N_ELEMENTS(vflat)*N_ELEMENTS(lflat)*N_ELEMENTS(vdisp)*N_ELEMENTS(Z0)
count=0

OPENW, outunit, 'gipsy_script.col', /GET_LUN
OPENW, ounit1, 'parameters.txt', /GET_LUN

SPAWN, 'rm -r profiles'
SPAWN, 'mkdir profiles'

PRINTF, outunit, '"RFITS ' + $
                 ' CTYPE4=hide' + $
                 ' FITSFILE=' + cube_name + '.fits' + $
                 ' OUTSET=' + cube_name + $
                 ' OKAY=y"'

PRINTF, outunit, '"RFITS ' + $
                 ' CTYPE3=hide' + $
                 ' FITSFILE=' + nhi_name + '.fits' + $
                 ' OUTSET=' + nhi_name + $
                 ' OKAY=y"'


radii_string = ''
FOR m=0, N_ELEMENTS(R)-1 DO BEGIN
    radii_string = radii_string + STRTRIM(STRING(R[m]),1) + ' '
    print, STRING(R[m])
ENDFOR

print, radii_string

PRINTF, outunit, '"ELLINT' + $
                 ' INSET=' + nhi_name + $
                 ' OPTION=' + $
                 ' RADII=' + radii_string + $
                 ' WIDTH=10' + $
                 ' PA=' + STRING(pa[0]) + $
                 ' INCL=0' + $
                 ' SEGMENTS=' + $
                 ' POS=' + center_pos + $
                 ' SUBPIX=' + $
                 ' MEDIAN=N' + $
                 ' FORMAT=' + $
                 ' PLOTOPT=' + $
                 ' FILENAME=' + ellint_filename + $
                 ' OVERWRITE=y"'

FOR i=0, N_ELEMENTS(incl)-1 DO BEGIN
FOR j=0, N_ELEMENTS(pa)-1 DO BEGIN
FOR k=0, N_ELEMENTS(lflat)-1 DO BEGIN
FOR l=0, N_ELEMENTS(vflat)-1 DO BEGIN
FOR n=0, N_ELEMENTS(vdisp)-1 DO BEGIN
FOR o=0, N_ELEMENTS(Z0)-1 DO BEGIN
	
	count=count+1
	vrot=vflat[l]*(1.-EXP(-R/lflat[k]))	; km/s
    ;nrings = size(R, /N_ELEMENTS)
    ;vrot=vflat[l]*(indgen(nrings) / float(nrings - 1) )
	;cgPLOT, R, vrot, XTITLE='R [arcsec]', YTITLE='V!Drot!N [km/s]', YRANGE=[0,MAX(vflat)]

	OPENW, ounit, 'profiles/profile_'+STRCOMPRESS(STRING(count), /REMOVE_ALL)+'.txt', /GET_LUN
	FOR m=0, N_ELEMENTS(R)-1 DO PRINTF, ounit, R[m], vrot[m]
	FREE_LUN, ounit

	model_name='model_'+STRCOMPRESS(STRING(count), /REMOVE_ALL)

	PRINTF, outunit, '"DELETE INSET=models/model_'+STRCOMPRESS(STRING(count), /REMOVE_ALL)+'; OK=y;"'

	PRINTF, outunit, '"DELETE INSET=models/model_'+STRCOMPRESS(STRING(count), /REMOVE_ALL)+'.smooth; OK=y;"'

	PRINTF, outunit, '"GALMOD BOX= CDENS= CMODE= DENS=file(' + ellint_filename $
    + ', 4, :)*1.e20 '+ $
	        'DRVAL3=1.419351555040E+9 EMPTY=N INCL='+$
                STRCOMPRESS(STRING(incl[i]), /REMOVE_ALL)+$
                ' INSET=' + cube_name + ' V LTYPE=1 NRESTR= NV= OUTSET=models/'+$
                model_name+' PA='+STRCOMPRESS(STRING(pa[j]), /REMOVE_ALL)+' '+$
            'POS=' + center_pos + ' RADII=' + radii_string +' VDISP='+$
                STRCOMPRESS(STRING(vdisp[n]), /REMOVE_ALL)+' VROT=file(profiles/profile_'+$
                STRCOMPRESS(STRING(count), /REMOVE_ALL)+'.txt,2,:)'+' VSYS=' +$
            vsys + ' Z0='+$
                STRCOMPRESS(STRING(Z0[o]), /REMOVE_ALL)+'"'

	PRINTF, outunit, '"SMOOTH AUTO=NO BOX= DECIM= INSET=models/model_'+$
                STRCOMPRESS(STRING(count), /REMOVE_ALL)+$
                ' V NEWBEAM=' + beamsizes + ' NEWPOSANG=0 OLDBEAM= OLDPOSANG= OUTSET=models/model_'+$
                STRCOMPRESS(STRING(count), /REMOVE_ALL)+'.smooth ; scale=;"'

	PRINTF, outunit, '"DELETE INSET=models/model_'+STRCOMPRESS(STRING(count), /REMOVE_ALL)+'; OK=y;"'
	

	PRINTF, outunit, '"WFITS INSET=models/model_'+$
                STRCOMPRESS(count, /REMOVE_ALL)+$
                '.smooth BOX= FITSFILE=models/model_'+$
                STRCOMPRESS(count, /REMOVE_ALL)+$
                '.smooth.FITS BITPIX=-32"'

    PRINTF, outunit, '"DELETE INSET=models/model_'+STRCOMPRESS(STRING(count), /REMOVE_ALL)+'.smooth; OK=y;"'

	PRINTF, ounit1, model_name+' '+STRCOMPRESS(STRING(incl[i]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(pa[j]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(lflat[k]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(vflat[l]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(vdisp[n]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(Z0[o]), /REMOVE_ALL)

	PRINTF, ounit1, model_name+' '+$
                STRCOMPRESS(STRING(incl[i]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(pa), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(lflat[k]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(vflat[l]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(vdisp[n]), /REMOVE_ALL)+' '+$
                STRCOMPRESS(STRING(Z0[o]), /REMOVE_ALL)

ENDFOR; o
ENDFOR; n
ENDFOR; l
ENDFOR; k
ENDFOR; j
ENDFOR; i

FREE_LUN, outunit
FREE_LUN, ounit1
END



