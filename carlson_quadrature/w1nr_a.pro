;+
; NAME:
;       W1NR_A
;
; PURPOSE:
;
;       Purpose of this function is to solve the equation for the weight
;       W1 (see Carlson, 1963). Value of W1 is used to calculate other weights
;       of directional cosines (first column of Tab.2 of Carlson). The 
;       equation is solved using Newton-Raphson method. The values of W1
;       agree with those of Carlson up to the numerical precision. Note that 
;       when using this routine, we compute the exact W1 (no need for epsilon_n
;       parameter, cf. Carlson, p.14).
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
;       w = W1NR_A(n)
;
; INPUTS:
;
;       n = Scalar, integer. This is the number of mu-levels. The number of rays
;           is n(n+2)/8.
;
; OUTPUTS:
;
;       wl         = Array(n/2), float. Weights of directional cosines.
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
;  EXAMPLE:
;
;       w = W1NR_A(8)
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October 2012.
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

FUNCTION w1nr_a, n

delta = 2.D/(n-1.D)
lm = n/2.
W1_1 = 4.D/(3.D*(n-1.D))
W1_0 = 1.
c = (n-2.D)/3.D

i = 0

omega = DBLARR(lm)

; PRINT, 'Newton-Raphson iteration (i, W1)'
WHILE ABS(W1_1 - W1_0) GT 1.D-7 DO BEGIN
  W1_0 = W1_1
  i += 1

  ; Function
  f = 0.D
  FOR l = 1, lm-1 DO $
    f = f + SQRT(W1_0^2 + (l-1)*delta)
  f = f - c

  ; Derivative
  fprim = 0.D
  FOR l = 1, lm-1 DO $
    fprim = fprim + W1_0/SQRT(W1_0^2 + (l-1)*delta)

  W1_1 = W1_0 - f/fprim

ENDWHILE

omega[0] = W1_1
FOR l = 2, lm-1 DO $
  omega[l-1] = SQRT(omega[0]^2 + (l-1)*delta)-TOTAL(omega[0:l-2])
omega[lm-1] = 1.D - TOTAL(omega)

RETURN, omega
END