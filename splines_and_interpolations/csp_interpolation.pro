;+
; NAME:
;       CSP_INTERPOLATION
;
; PURPOSE:
;
;       This function returns the coefficients needed for the cubic spline 
;       interpolation of a given set function y that is sampled at n+2 
;       discrete value of the variable x. The number of returned indices is
;       4(n+1). The PERIODIC boundaries are imposed. To evaluate the 
;       interpolation function one may use CSP_EVALUATE as in the example
;       below.
;
; AUTHOR:
;
;       Nikola Vitas
;       Instituto de Astrofisica de Canarias (IAC)
;       C/ Via Lactea, s/n
;       E38205 - La Laguna (Tenerife), Espana
;       Email: n.vitas@iac.es, nikola.vitas@gmail.com
;       Homepage: nikolavitas.blogspot.com
;
; CATEGORY:
;
;       Numerics.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
; 
;      ; The exact function for comparison.
;      nxa = 241
;      period = 2.
;      xa = DINDGEN(nxa)/100-1.2
;      ya = SIN(2.*!pi*xa/period)  
;      
;      ; The 'measurement' samples the function at 9 points (the array must be
;      ; monotonically increasing).
;      x = DOUBLE([-1., -0.8, -0.6, -0.45,  0.0, 0.3, 0.5, 0.6, 1.])
;      nx = N_ELEMENTS(x)
;      
;      ; Force the exact periodicity, y[0] = y[nx-1]
;      y = [SIN(2*!pi*x[0:nx-2]/period), SIN(2*!pi*x[0]/period)]
;        
;      ; Compute the coefficients of the approximative function S and evaluate S 
;      ; at the same grid used to compute the exact function.
;      nx0 = nxa
;      cp = CSP_INTERPOLATION(x, y)
;      x0 = FINDGEN(120)/4.-3
;      y0 = CSP_EVALUATE(x0, x, cp)     
;
; DEPENDENCIES:
;
;       THREEDIAGONAL_PERIODIC_SYSTEM
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, September, 2013.
;-
;================================================================================
; CSP_INTERPOLATION, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; This software is provided by NV ''as is'' and any express or implied warranties, 
; including, but not limited to, the implied warranties of merchantability and 
; fitness for a particular purpose are disclaimed. In no event shall NV be liable 
; for any direct, indirect, incidental, special, exemplary, or consequential 
; damages (including, but not limited to, procurement of substitute goods or 
; services; loss of use, data, or profits; loss of use, data, or profits; or 
; business interruption) however caused and on any theory of liability, whether 
; in contract, strict liability, or tort (including negligence or otherwise) 
; arising in any way out of the use of this software, even if advised of the 
; possibility of such damage.
;================================================================================

FUNCTION csp_interpolation, x, y

;--------------------------------------------------------------------------------
; Define the variables
;--------------------------------------------------------------------------------
; n+2 data points
; n+1 segment between them
; x indices go from 0 to n+1
; h - difference x[i]-x[i-1], n+1 value, indices from 0, to n
; w - n+1 values, indices from 0 to n 
; q - n+1 values, indices from 0 to n
; a, b, c, d - n+1 values, indices form 0 to n
;--------------------------------------------------------------------------------
n = SIZE(x, /N_ELEMENTS)-2
h = DBLARR(n+1)
w = DBLARR(n+1)
q = DBLARR(n+1)
d = DBLARR(n+1)
b = DBLARR(n+1)
a = y[0:n]
c = DBLARR(n+1)

; Test periodicity
IF ABS((y[0]-y[n+1])) GE (machar()).eps*ABS(y[0]) THEN $
  MESSAGE, 'ERROR: Function must be periodic in the sense: y[0] == y[n]'

; Compute the differences h
FOR i = 0, n DO h[i] = x[i+1] - x[i]

; Compute the quantities w and q
w[0] = 2*(h[0] + h[n]) 
q[0] = 3*(y[1] - y[0])/h[0] - 3*(y[0] - y[n])/h[n]
FOR i = 1, n DO BEGIN
  w[i] = 2*(h[i] + h[i-1])
  q[i] = 3*(y[i+1] - y[i])/h[i] - 3*(y[i] - y[i-1])/h[i-1]
ENDFOR

lower = SHIFT(h, 1)
upper = h

c = THREEDIAGONAL_PERIODIC_SYSTEM(lower, w, upper, q)

d[n] = (c[0] - c[n])/(3*h[n])
FOR i = 0, n-1 DO d[i] = (c[i+1] - c[i])/(3*h[i])

b[n] = (y[n+1] - y[n])/(h[n]) - c[n]*h[n] - d[n]*h[n]^2
FOR i = 0, n-1 DO b[i] = (y[i+1] - y[i])/(h[i]) - c[i]*h[i] - d[i]*h[i]^2
 
; ; In the output we omit the last c-coefficient (it's 0 anyway).
output = {a:a, b:b, c:c, d:d}

RETURN, output
END