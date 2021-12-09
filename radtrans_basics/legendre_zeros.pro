;+
; NAME:
;       LEGENDRE_ZEROS
;
; PURPOSE:
;
;       This function computes the zeros of Legendre polynomial
;       specified by its coefficients. It used the method of chords. 
;      
; AUTHOR:
;
;       Nikola Vitas
;       Instituto de Astrofisica de Canarias (IAC)
;       C/ Via Lactea, s/n
;       E38205 - La Laguna (Tenerife), Espana
;       Email: n.vitas@iac.es
;       Homepage: nikolavitas.blogspot.com
;
; CATEGORY:
;
;       Numerics.
;
; CALLING SEQUENCE:
;
;       x = LEGENDRE_ZEROS(nx = nx, epsilon = epsilon, coeff = coeff)
;
; INPUTS:
; 
; OUTPUTS:
;
;       x     = Array, float. Contains the zeros.
;
; INPUT KEYWORDS:
;
;       nx      = Scalar, integer. Number of intervals in which the initial
;                 search is performed.
;
;       epsilon = Scalar, float. Required precision.
;
;       coeff   = Array, float. Coefficients of the given Legendre polynomial.
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, January 2006.
;       - small changes and fixes, NV, June 2013.
;-
;================================================================================
; LEGENDRE_POLYNOMIALS and other belonging routines by Nikola Vitas are licensed 
; under a Creative Commons  Attribution-NonCommercial-ShareAlike 3.0 Unported 
; License.
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

FUNCTION legendre_zeros, nx = nx, epsilon = epsilon, coeff = coeff

n = SIZE(coeff, /N_ELEMENTS) - 1
x = FINDGEN(nx+1)
x = x/nx
f = DBLARR(nx+1)
zeros = DBLARR(n)
iz = 0
interval = DBLARR(n, 2)

f = POLY(x, coeff)

; We define approximate intervals in which we'll look for the zeros.
k = 0

FOR i = 1, nx-1 DO BEGIN
  IF (f[i-1] * f[i] LT 0.) THEN BEGIN
    interval[k, 0] = x[i-1]
    interval[k, 1] = x[i]
    k = k + 1
  ENDIF

; Check if we found a zero. If yes, we add it to the array.
  IF (f[i] EQ 0.) THEN BEGIN
    zeros[iz] = -x[i]
    zeros[iz+1] = x[i]
    iz = iz + 2
  ENDIF
ENDFOR

; If the order n is odd, one zero is = 0, but that root we cannot find
; because it is on the edge of the interval. So we add it.

IF FLOAT(n)/2 NE FLOOR(FLOAT(n)/2) THEN BEGIN
  zeros[iz] = 0.d
  iz = iz + 1
ENDIF

FOR ik = 0, k-1 DO BEGIN
  alpha = interval[ik, 0]
  betha = interval[ik, 1]

  ix = 0

  x1 = (betha * ABS(POLY(alpha, coeff)) + alpha * ABS(POLY(betha, coeff)))/ $
               (ABS(POLY(alpha, coeff)) +         ABS(POLY(betha, coeff)))
  IF (POLY(alpha, coeff) * POLY(x1, coeff) LT 0) THEN betha = x1 $
                                                 ELSE alpha = x1
                                                   
  x2 = (betha * ABS(POLY(alpha, coeff)) + alpha * ABS(POLY(betha, coeff)))/ $
               (ABS(POLY(alpha, coeff)) +         ABS(POLY(betha, coeff)))
  IF (POLY(alpha, coeff) * POLY(x2, coeff) LT 0) THEN betha = x2 $
                                                 ELSE alpha = x2

  WHILE (ABS(x2 - x1) GT epsilon) DO BEGIN
    x1 = x2
    x2 = (betha * ABS(POLY(alpha, coeff)) + ALPHA * ABS(POLY(betha, coeff)))/ $
                 (ABS(POLY(alpha, coeff)) +         ABS(POLY(betha, coeff)))
    IF (POLY(alpha, coeff) * POLY(x2, coeff) LT 0) THEN betha = x2 $
                                                   ELSE alpha = x2
    ix = ix + 1
  ENDWHILE
  
  zeros[iz] = -x2
  zeros[iz+1] = x2
  iz = iz + 2
ENDFOR

; Sort the array with zeros.
zsort = SORT(zeros)
zeros = zeros[zsort]

; Graphical test
; PLOT, x, f, background = 16777215, color = 0
; yzero = DBLARR(n)
; OPLOT, zeros, yzero, psym = 6, symsize = 1.2, color = 3858760
; PRINT, zeros

RETURN, zeros

END