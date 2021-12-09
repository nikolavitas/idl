FUNCTION pf_from_eqc, molname, t

RESTORE, 'catalogue_of_molecules.sav'

;------------------------------------------------------------------------------
; Constants
;------------------------------------------------------------------------------
kkev  = 8.61734D-5           ; eV K^-1
kkerg = 1.38066D-16          ; erg K^-1
cc    = 2.99792D10           ; cm s^-1
mel   = 9.10939D-28          ; g
hh    = 6.62607D-27          ; erg s
amu = 1.660538921D-24
;------------------------------------------------------------------------------

  ; Find the index of molecule im in the catalogue
  ind = WHERE(molname EQ cat.name)
  IF ind EQ -1 THEN BEGIN
    PRINT, 'Data for ' + molname + ' is not available.' 
    STOP
  ENDIF ELSE BEGIN
    ; Array of indices
    mols_in_cat = ind
    ; Fill the m container
    name = molname
    charge = cat[ind].charge
    omega = cat[ind].omega
    d0 = cat[ind].d0
    
    lambda32 = ((2*!pi*cat[ind].mab*amu*kkerg*t)/(hh^2))^1.5
    
    
    ; USER SHOULD BE ALLOWED TO CHOOSE. IF S/HE ISN'T, THEN THIS SHOULD 
    ; BE MADE MORE TRANSPARENT
    ; Get the equilibrium constant for pressure at given T (Sauval and Tatum data is used)
    IF TOTAL(TOTAL(omega) LE 2) THEN BEGIN 
      IF ind LE 293 THEN source = 's' ELSE source = 't' 
    ENDIF ELSE BEGIN
      source = 'rossi'
    ENDELSE
 
    k = MOLECULES(t, molname, data = source, type = 'kp')
    k = k/(kkerg*t)
  ENDELSE  
  
  qaqb = 1.D0
  FOR i = 0, 91 DO begin
    u1 = (pf(t, i+1, data = "i")).u1
      qaqb *= u1^cat[ind].omega[i]
  ENDFOR
  
  pfm = k * 0.D0
  pfm = lambda32 * (qaqb / k) * EXP(-D0/kkev/t)   

  stop
  
return, pfm
END
