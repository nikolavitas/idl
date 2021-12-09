;+
; NAME:
;       LEGENDRE_POLYNOMIALS
;
; PURPOSE:
;
;       Purpose of this function is to compute the coefficients of Legendre
;       polynomials of a given order, the zeros of the polynomial and their
;       weights. This is the main routine that wraps together the results of
;       these three separate problems. The algorthm is written mainly following
;       "Matematicka obrada astronomskih posmatranja" (Djurovic, 1979, 
;       "Mathematical analysis of astronomical observations").
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
;       x = LEGENDRE_POLYNOMIALS(n)
;
; INPUTS:
;
;       n     = Scalar, integer. Order of Legendre polynomial.
; 
; OUTPUTS:
;
;       x     = Structure with three components: x.coeffs contain the 
;               coefficients of the polynomial; x.zeros its zeros; and
;               x.weights their weights.
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
;  EXAMPLE:
;
;       x = LEGENDRE_POLYNOMIALS(6)
;       FOR i = 0, 5 DO PRINT, i, x.zeros(i), x.weights(i)
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, January 2006.
;       - Various small fixes, NV, June 2013.
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

FUNCTION legendre_polynomial, n


result = {zeros:FLTARR(n), weights:FLTARR(n), coeffs:FLTARR(n+1)}

result.coeffs = LEGENDRE_COEFFICIENTS(n)
result.zeros = LEGENDRE_ZEROS(nx = 1000, epsilon = 1.d-12, coeff = result.coeffs)
result.weights = LEGENDRE_WEIGHTS(zero = result.zeros)

RETURN, result

END




