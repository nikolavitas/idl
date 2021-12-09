;+
; NAME:
;       LEGENDRE_WEIGHTS
;
; PURPOSE:
;
;       This function computes the weights of Legendre polynomial
;       zeros by solving a system of linear equations. The system
;       is solved using Cramer's method.
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
;       x = LEGENDRE_WEIGHTS(zero = zero)
;
; INPUTS:
; 
; OUTPUTS:
;
;       x     = Array, float. Contains the weights.
;
; INPUT KEYWORDS:
;
;       zero  = Array, float. An array containing zeros of a Legendre
;               polynomial.

; OUTPUT KEYWORDS:
;
;  EXAMPLE:
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

FUNCTION legendre_weights, zero = zero

n  = SIZE(zero, /N_ELEMENTS)
yc = DBLARR(n)
c  = DBLARR(n)
u  = DBLARR(n, n)
t  = DBLARR(n, n)
d  = DBLARR(n)
ds = 0.D

FOR i = 0, n-1 DO BEGIN
  u[*, i] = zero^i
  IF FLOAT(i)/2 EQ FLOOR(FLOAT(i)/2) THEN BEGIN
    yc[i] = 2.D0/(i+1)
  ENDIF ELSE BEGIN
    yc[i] = 0.D0
  ENDELSE
ENDFOR

ds = DETERM(u)

FOR i = 0, n-1 DO BEGIN
  t = u
  t[i, *] = yc
  d[i] = DETERM(t)
  c[i] = d[i]/ds
ENDFOR

RETURN, c

END