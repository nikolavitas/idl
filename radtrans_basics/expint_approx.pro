;+
; NAME:
;       EXPINT_APPROX
;
; PURPOSE:
;
;       This function evaluates exponential-integral function for given set
;       of arguments (x) and specified order of the integral (n). It uses
;       polynomial fit of E1 and the recurrent formula to compute En, n>1.
;       The coefficients of the fit are from Abramowitz and Stegun (1964).
;       IDL has an intrinsic routine that computes these functions directly
;       from their definition what is, of course, more accurate. This routine
;       thus serves only for the comparison purposes.
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
;       y = EXPINT_APPROX(n, x)
;
; INPUTS:
;
;       n     = Scalar, integer. Order of the integral that is computed.
;
;       x     = Array, float. Contains the arguments.
; 
; OUTPUTS:
;
;       y     = Array, float. Contains the function values.
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, March 2008.
;       - minor changes, NV, June 2013.
;-
;================================================================================
; EXPINT_APPROX by Nikola Vitas is licensed under a Creative Commons Attribution-
; NonCommercial-ShareAlike 3.0 Unported License.
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

FUNCTION expint_approx, n, x

; The best fit coefficients from Abramowitz and Stegun (1964).
a = [0.2677737343D0,  8.6347608925D0, 18.0590169730D0, 8.5733287401D0, 1.0000000000D0]
b = [3.9584969228D0, 21.0996530827D0, 25.6329561486D0, 9.5733223454D0, 1.0000000000D0]
c = [-0.57721566D0, 0.99999193D0, -0.24991055D0, $
      0.05519968D0, -0.00976004D0, 0.00107857D0]
 
nx = N_ELEMENTS(x)
eint = DBLARR(nx)

; Separate fits for x < 1 and x > 1.
int = WHERE(x GT 1.)
IF TOTAL(int) GT 0 THEN $
  eint[int] = EXP(-x[int])*POLY(x[int], a)/(x[int] * POLY(x[int], b))

int = WHERE(x LE 1.)
IF TOTAL(int) GT 0 THEN $
  eint[int] = -ALOG(x[int]) + POLY(x[int], c)
  
; Recurrent formula.
IF n GT 1 THEN $
  FOR i = 1, n-1 DO eint = (EXP(-x) - x * eint)/i 

RETURN, eint
END