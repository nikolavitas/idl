;+
; NAME:
;       MU1NR_B
;
; PURPOSE:
;
;       Purpose of this function is to compute the first-order directional
;       cosine of the quadrature Set B (see Carlson, 1963). We solve Eqs.
;       9, 10 and 11 numerically by the Newton-Raphson formula. The result
;       is in agreement with the values tabulated by Carlson. 
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
;       Radiative transfer, numerics.
;
; CALLING SEQUENCE:
;
;       w = MU1NR_B(omega, n)
;
; INPUTS:
;
;       omega = Array, float. These are the weights of directionals that are
;               precomputed using W1NR_B.PRO.
; 
;       n = Scalar, integer. This is the number of mu-levels. 
;
; OUTPUTS:
;
;       mu = Array, float. Directional cosines.
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
;  EXAMPLE:
;
;       omega = W1NR_B(8)
;       mu = MU1NR_B(omega, 8)
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, November 2012.
;-
;================================================================================
; Carlson Quadrature Suite by Nikola Vitas is licensed under a Creative Commons 
; Attribution-NonCommercial-ShareAlike 3.0 Unported License.
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

FUNCTION mu1nr_b, omega, n

lm = n/2.D
mu = DBLARR(lm)
mu1_1 = 0.50D
mu1_0 = 0.51D
c = 1.D/2.D
i = 0

; PRINT, 'Newton-Raphson iteration (i, W1)'
WHILE ABS(mu1_1 - mu1_0) GT 1.D-7 DO BEGIN
  mu1_0 = mu1_1
  i += 1

  ; Function
  f = 0.D
  FOR l = 1, lm DO $
    f = f + omega[l-1]*SQRT(mu1_0^2*(1.-(6./n)*(l-1)) + (2./n)*(l-1))
  f = f - c

  ; Derivative
  fprim = 0.D
  FOR l = 1, lm DO $
    fprim = fprim + omega[l-1]*mu1_0*(1.-(6./n)*(l-1))/$
                    SQRT(mu1_0^2*(1.-(6./n)*(l-1)) + (2./n)*(l-1))

  mu1_1 = mu1_0 - f/fprim
ENDWHILE

mu[0] = mu1_1
FOR l = 2, lm DO $
  mu[l-1] = SQRT(mu1_0^2*(1.-(6./n)*(l-1)) + (2./n)*(l-1))

RETURN, mu
END