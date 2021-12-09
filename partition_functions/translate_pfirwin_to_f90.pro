PRO translate_pfirwin_to_f90

RESTORE, 'pf_irwin.sav'

OPENW, 1, 'pf_irwin.F90'

PRINTF, 1, '  SUBROUTINE pf_irwin_data'


FOR iz = 0, 91 DO BEGIN
FOR iu = 0, 2 DO BEGIN


  ius = STRCOMPRESS(STRING(iu+1), /rem)
  izs = STRCOMPRESS(STRING(iz+1), /rem)
  data_string = '  irwin_pf_coeffs('+izs+', '+ius+', :) = '
  
  data_string = data_string + '(/ '
  FOR it = 0, 2 DO BEGIN
    its = STRCOMPRESS(STRING(it+1), /rem)
    data_string = data_string + STRING(pfcoeffs[iz, iu, it], format = '(E15.8)') + ', '
  ENDFOR  
  PRINTF, 1, data_string + ' &'
  data_string = STRING(' ', format = '(A33)')
  FOR it = 3, 4 DO BEGIN
    its = STRCOMPRESS(STRING(it+1), /rem)
    data_string = data_string + STRING(pfcoeffs[iz, iu, it], format = '(E15.8)') + ', '
  ENDFOR  
  it = 5
  data_string = data_string + STRING(pfcoeffs[iz, iu, it], format = '(E15.8)') + '/) '
  PRINTF, 1, data_string
  
  
ENDFOR
PRINTF, 1, ''
ENDFOR

PRINTF, 1, '  END SUBROUTINE pf_irwin_data'

CLOSE, 1
stop

END
