PRO legendre_figure

!P.BACKGROUND = GETCOLOR('white')
!P.COLOR = GETCOLOR('black')
!P.MULTI = 0
PLOTSYM, 8, 1.4, /fill
WINDOW, 1, xsi = 500, ysi = 350

x = FINDGEN(1001)/1000
n = 6
y = LEGENDRE_POLYNOMIAL(n)
PLOT, x, POLY(x, y.coeffs), /noda, xr = [0, 1], yr  = [-0.5, 1], /xs, /ys, $
      xtit = 'x', ytit  = 'y', tit = 'Legendre polynomials', charsi = 1.2

OPLOT, x, POLY(x, y.coeffs), thick = 1, color = GETCOLOR('red')
FOR i = 0, n-1 DO OPLOT, [y.zeros[i]], [0], psym = 8, color = GETCOLOR('red')

n = 10
y = LEGENDRE_POLYNOMIAL(n)
OPLOT, x, POLY(x, y.coeffs), thick = 1, color = GETCOLOR('green')
FOR i = 0, n-1 DO OPLOT, [y.zeros[i]], [0], psym = 8, color = GETCOLOR('green')

n = 3
y = LEGENDRE_POLYNOMIAL(n)
OPLOT, x, POLY(x, y.coeffs), thick = 1, color = GETCOLOR('blue')
FOR i = 0, n-1 DO OPLOT, [y.zeros[i]], [0], psym = 8, color = GETCOLOR('blue')

n = 4
y = LEGENDRE_POLYNOMIAL(n)
OPLOT, x, POLY(x, y.coeffs), thick = 1, color = GETCOLOR('orange')
FOR i = 0, n-1 DO OPLOT, [y.zeros[i]], [0], psym = 8, color = GETCOLOR('orange')

n = 5
y = LEGENDRE_POLYNOMIAL(n)
OPLOT, x, POLY(x, y.coeffs), thick = 1, color = GETCOLOR('violet')
FOR i = 0, n-1 DO OPLOT, [y.zeros[i]], [0], psym = 8, color = GETCOLOR('violet')


image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'legendre_example.jpeg', image_jpeg, TRUE=1, QUALITY=100


END