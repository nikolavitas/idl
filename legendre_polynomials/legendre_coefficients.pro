;+
; NAME:
;       LEGENDRE_COEFFICIENTS
;
; PURPOSE:
;
;       This function computes the coefficients of Legendre polynomial
;       of a given order. It uses recurrent formula. For more details check:
;       http://nikolavitas.blogspot.com.es/2013/06/gauss-quadrature-and-legendre.html
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
;       x = LEGENDRE_COEFFICIENTS(n)
;
; INPUTS:
; 
; OUTPUTS:
;
;       x     = Array, float. Contains the coefficients.
;
; INPUT KEYWORDS:
;
;       n       = Scalar, integer. Number of intervals in which the initial
;                 search is performed.
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

FUNCTION legendre_coefficients, n

coeff0 = DBLARR(n+1)
coeff1 = DBLARR(n+1)
coeff2 = DBLARR(n+1)

; Coefficients of the zero-order Legendre polynomial
coeff0[0] = 1

; Coefficients of the first-order Legendre polynomial
coeff1[1] = 1

; Using recurrent formula to find the coefficient of the n-th order
FOR i = 2, n DO BEGIN
  coeff2 = (1.D0/i) * ((2*i - 1)*SHIFT(coeff1, 1) - (i - 1)*coeff0)
  coeff0 = coeff1
  coeff1 = coeff2
ENDFOR

RETURN, coeff2

END