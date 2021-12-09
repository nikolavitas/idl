PRO example_abundances

; This routine has been use to generate figures for the blog.
; It's not meant to be shared.

!p.background = GETCOLOR('white')
!p.color = GETCOLOR('black')
WINDOW, xsi = 500, ysi = 400

a09 = LOAD_ABUNDANCES(s = 'a09')
gs98 = LOAD_ABUNDANCES(s = 'gs98')

mosta = (REVERSE(SORT(a09)))[0:11]

sym = LOAD_LIST_OF_ELEMENTS()

z = INDGEN(92)+1


PLOT, z, a09, xtit = 'Z', ytit = 'A', yr = [-1, 14.5], /ys, xr = [-2, 94], /xs, charsi = 1.4, $
      title = 'Solar chemical composition (A09)', /noda
PSYM10, z, a09, color = GETCOLOR('salmon'), y0 = -1., thick = 2

; WRITE_PNG, 'example_abundances_0.png', tvrd(/true)


PLOT, z[0:39], gs98[0:39], xtit = 'Z', ytit = 'A', yr = [-1, 14.5], /ys, xr = [-2, 42], /xs, charsi = 1.4, $
      title = 'Solar chemical composition', /noda
PSYM10, z[0:39]+0.1, gs98[0:39], color = GETCOLOR('dark gray'), y0 = -1.
PSYM10, z[0:39]-0.1, a09[0:39], color = GETCOLOR('salmon'), y0 = -1.
OPLOT, [30, 33], [10, 10], color = GETCOLOR('dark gray')
XYOUTS, 34, 9.9, 'GS98', charsi = 1.4
OPLOT, [30, 33], [11, 11], color = GETCOLOR('salmon')
XYOUTS, 34, 10.9, 'A09', charsi = 1.4

odd = 1
FOR i = 0, 39 DO BEGIN
  IF WHERE(mosta EQ i) NE -1 THEN BEGIN
    IF odd EQ 0 THEN BEGIN
      XYOUTS, z[i]-0.5, 12.5, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [a09[i]+0.3, 12.3], color = GETCOLOR('light gray')
      odd = 1
    ENDIF ELSE BEGIN
      XYOUTS, z[i]-0.5, 13.2, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [a09[i]+0.3, 13.0], color = GETCOLOR('light gray')
      odd = 0
    ENDELSE
    
    
  ENDIF
ENDFOR
; WRITE_PNG, 'example_abundances_1.png', tvrd(/true)


natoms = 84
aw = LOAD_ATOMIC_WEIGHTS()
abundance = gs98
relabund = 10.D^(abundance)/TOTAL(10.D^(abundance))
xx_gs98 = relabund[0]*aw[0]/(TOTAL(relabund*aw))
yy_gs98 = relabund[1]*aw[1]/(TOTAL(relabund*aw))
zzz_gs98 = relabund[2:natoms-1]*aw[2:natoms-1]/(TOTAL(relabund*aw))
zz_gs98 = TOTAL(zzz_gs98)
mix = abundance

abundance = a09
relabund = 10.D^(abundance)/TOTAL(10.D^(abundance))
xx_a09 = relabund[0]*aw[0]/(TOTAL(relabund*aw))
yy_a09 = relabund[1]*aw[1]/(TOTAL(relabund*aw))
zzz_a09 = relabund[2:natoms-1]*aw[2:natoms-1]/(TOTAL(relabund*aw))
zz_a09 = TOTAL(zzz_a09)

mix[7] = abundance[7]
abundance = mix
relabund = 10.D^(abundance)/TOTAL(10.D^(abundance))
xx_mix = relabund[0]*aw[0]/(TOTAL(relabund*aw))
yy_mix = relabund[1]*aw[1]/(TOTAL(relabund*aw))
zzz_mix = relabund[2:natoms-1]*aw[2:natoms-1]/(TOTAL(relabund*aw))
zz_mix = TOTAL(zzz_mix)



stop
mosta = (REVERSE(SORT(a09)))[2:11]
natoms = 40
PLOT, z[2:natoms-1], zzz_a09, xtit = 'Atomic number Z', ytit = 'Z!Delement!N/Z', yr = [0, 0.6], /ys, xr = [-2, 42], /xs, charsi = 1.4, $
      title = 'Elemental contributions to the mass fraction Z', /noda
PSYM10, z[2:natoms-1]+0.1, zzz_gs98/zz_gs98, color = GETCOLOR('dark gray'), y0 = 0.
PSYM10, z[2:natoms-1]-0.1, zzz_a09/zz_a09, color = GETCOLOR('salmon'), y0 = 0.
OPLOT, z[2:natoms-1], zzz_a09/zz_a09, color = GETCOLOR('red'), psym = 1
OPLOT, z[2:natoms-1], zzz_gs98/zz_gs98, psym = 1

OPLOT, [30], [0.4], color = GETCOLOR('dark gray'), psym = 1
XYOUTS, 31, 0.395, 'GS98', charsi = 1.4
OPLOT, [30], [0.35], color = GETCOLOR('salmon'), psym = 1
XYOUTS, 31, 0.345, 'A09', charsi = 1.4

odd = 1

FOR i = 2, 39 DO BEGIN
  IF WHERE(mosta EQ i) NE -1 THEN BEGIN
    IF odd EQ 0 THEN BEGIN
      XYOUTS, z[i]-0.5, 0.5, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [max([zzz_a09[i-2]/zz_a09, zzz_gs98[i-2]/zz_gs98])+0.02, 0.48], col = GETCOLOR('light gray')
      odd = 1
    ENDIF ELSE BEGIN
      XYOUTS, z[i]-0.5, 0.55, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [max([zzz_a09[i-2]/zz_a09, zzz_gs98[i-2]/zz_gs98])+0.02, 0.52], col = GETCOLOR('light gray')
      odd = 0
    ENDELSE
    
    
  ENDIF
ENDFOR

; WRITE_PNG, 'example_abundances_3.png', tvrd(/true)



z = INDGEN(40)+1
a09 = LOAD_ABUNDANCES(z, s = 'a09', /ppm)
gs98 = LOAD_ABUNDANCES(z, s = 'gs98', /ppm)
rel = a09/gs98

PLOT, z, rel, xtit = 'Z', ytit = 'A09/GS98', yr = [0, 2], /ys, xr = [-2, 42], /xs, $
      charsi = 1.4, title = 'A09/GS98 Ratio', /noda

OPLOT, [-2, 42], [1, 1], col = GETCOLOR('light gray')
PSYM10, z, rel,  col = GETCOLOR('gray'), y0 = 1.
OPLOT, z, rel,  psym = 1

odd = 1
even = 1
FOR i = 0, 39 DO BEGIN
  IF WHERE(mosta EQ i) NE -1 THEN BEGIN
    IF odd EQ 0 THEN BEGIN
      XYOUTS, z[i]-0.5, 1.8, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [rel[i]+0.05, 1.75], col = GETCOLOR('light gray')
      odd = 1
    ENDIF ELSE BEGIN
      XYOUTS, z[i]-0.5, 1.6, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [rel[i]+0.05, 1.55], col = GETCOLOR('light gray')
      odd = 0
    ENDELSE    
  ENDIF ELSE BEGIN
    IF even EQ 0 THEN BEGIN
      XYOUTS, z[i]-0.5, 0.2, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [rel[i]-0.05, 0.30], col = GETCOLOR('light gray')
      even = 1
    ENDIF ELSE BEGIN
      XYOUTS, z[i]-0.5, 0.4, sym[i], charsi = 1.2
      OPLOT, [z[i], z[i]], [rel[i]-0.05, 0.50], col = GETCOLOR('light gray')
      even = 0
    ENDELSE    
  ENDELSE
ENDFOR
; WRITE_PNG, 'example_abundances_2.png', tvrd(/true)



stop

END
