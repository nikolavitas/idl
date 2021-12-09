;+
; NAME:
;       CSN_INTERPOLATION
;
; PURPOSE:
;
;       This function returns the coefficients needed for the cubic spline 
;       interpolation of a given set function y that is sampled at n+1 
;       discrete value of the variable x. The number of returned indices is
;       4n. At the boundary condition the second derivative of the interpolation
;       function is set to 0 (the NATURAL splines). The method is coded after 
;       the book of D.S.G. Pollock (Handbook of Time Series Analysis, Signal 
;       Processing and Dynamics, ISBN: 978-0-12-560990-6, 
;       http://www.sciencedirect.com/science/book/9780125609906 ).
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
;      nxa = 201
;      xa = DINDGEN(nxa)/100-1.
;      ya = CS_FUNCTION(xa)
;      
;      ; The 'measurement' samples the function at 9 points (the array must be
;      ; monotonically increasing).
;      x = DOUBLE([-1., -0.8, -0.6, -0.45,  0.0, 0.3, 0.5, 0.6, 1.])
;      nx = N_ELEMENTS(x)
;      y = 0.5 * x * COS(1.5*!pi*x + 0.5)
;      
;      ; Compute the coefficients of the approximative function S and evaluate S 
;      ; at the same grid used to compute the exact function.
;      x0 = xa
;      nx0 = nxa
;      co = CSN_INTERPOLATION(x, y)
;      
;      y0 = DBLARR(nx0)
;      FOR i = 0, nx0-1 DO BEGIN
;        j = VALUE_LOCATE(x, x0[i]) 
;        IF (j GE 0) AND (j LT nx-1) THEN $
;          y0[i] = co.a[j] + co.b[j]*(x0[i]-x[j]) + co.c[j]*(x0[i]-x[j])^2 + co.d[j]*(x0[i]-x[j])^3
;        IF (j GE 0) AND (j EQ nx-1) THEN $
;          IF (x[j] EQ x0[i]) THEN $
;            y0[i] = y[j]
;      ENDFOR
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, February 2005.
;       - Some improvments, August, 2013.
;-
;================================================================================
; INTERPOLATE_CSN, IDL routine by Nikola Vitas is licensed under a Creative 
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

FUNCTION csn_interpolation, x, y

;--------------------------------------------------------------------------------
; Define the variables
;--------------------------------------------------------------------------------
; n - number of segments between the elements of the x-array
; h - difference x[i]-x[i-1], n values, indices from 0, to n-1
; w - n-1 values, indices from 1 to n-1 (from 2 to n-1 after computing w')
; q - n-1 values, indices from 1 to n-1 (from 2 to n-1 after computing q')
; a, b, c, d - n-1 values, indices form 0 to n-1
;--------------------------------------------------------------------------------
n = SIZE(X, /N_ELEMENTS)-1
h = DBLARR(n)
w = DBLARR(n)
q = DBLARR(n)
d = DBLARR(n)
c = DBLARR(n+1)
b = DBLARR(n)
a = y[0:n-1]
c = DBLARR(n+1)

; Compute the differences h
FOR i = 0, n-1 DO h[i] = x[i+1] - x[i]

; Compute the quantities w and q
FOR i = 1, n-1 DO BEGIN
  w[i] = 2*(x[i+1] - x[i-1])
  q[i] = 3*(y[i+1] - y[i])/h[i] - 3*(y[i] - y[i-1])/h[i-1]
ENDFOR

; Evaluate w' and q': We don't keep the old values as we don't
; need them anywhere else
FOR i = 2, n-1 DO BEGIN
  w[i] = w[i] - h[i-1]^2/w[i-1]
  q[i] = q[i] - q[i-1]*h[i-1]/w[i-1]
ENDFOR

; Compute the c-coefficients
c[n] = 0
c[n-1] = q[n-1]/w[n-1]
FOR i = 2, n-1 DO c[n-i] = (q[n-i] - h[n-i]*c[n-i+1])/w[n-i]

; Set the boundaries
d[0] = c[1]/(3*h[0])
b[0] = (y[1] - y[0])/h[0] - c[1]*h[0]/3

; Compute d's and b's
FOR i = 1, n-1 DO BEGIN
  d[i] = (c[i+1] - c[i])/(3*h[i])
;   c[i] = ctemp[i]
  b[i] = (c[i] + c[i-1])*h[i-1] + b[i-1]
ENDFOR

; In the output we omit the last c-coefficient (it's 0 anyway).
output = {a:a, b:b, c:c[0:n-1], d:d}

RETURN, output
END