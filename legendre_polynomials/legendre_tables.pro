PRO legendre_tables

FOR n = 2, 33 DO BEGIN
  y = LEGENDRE_POLYNOMIAL(n)
  OPENW, 1, 'legendre_zeros_n'+NUM2STR(n)+'.dat'
  FOR i = 0, n-1 DO BEGIN
    PRINTF, 1, i+1, y.zeros[i], y.weights[i], format = '(I4, 2F32.24)'
  ENDFOR
  CLOSE, 1
  SAVE, y, file = 'legendre_n'+NUM2STR(n)+'.sav'
ENDFOR

END