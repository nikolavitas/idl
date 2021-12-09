;+
; NAME:
;       KOURGANOFF_PLOTS
;
; PURPOSE:
;
;       This function plots Kourganoff graphs that show the result of the 
;       lambda operator applied to various elementary source functions.
;       See Kourganoff (1952, pp.50-55) and Rutten (2003, p.83). 
;       
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
;       Figures.
;
; CALLING SEQUENCE:
;
;       KOURGANOFF_PLOTS
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
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, May 2013.
;-
;================================================================================
; KOURGANOFF_PLOTS and other belonging routines by Nikola Vitas are licensed 
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

PRO kourganoff_plots

!P.BACKGROUND = 16777215
!P.COLOR = 0
!P.MULTI = 0
WINDOW, 0, xsi = 500, ysi = 350

ntau = 301
tau = DINDGEN(ntau)/100

;-------------------------------------------------------------------------------
; Constant source function:
;-------------------------------------------------------------------------------
s = tau/tau
j = 1. - 0.5*EXPINT(2, tau)
PLOT, tau, s, yr = [0, 1.2], xr = [0, 3], /ys, /xs, $
  xtit = '!4s!3', ytickint = 0.5, xtickint = 1., charsi = 1.2
OPLOT, tau, j, lines = 2
XYOUTS, 1, 1.02, 'S = 1', charsi = 1.4
XYOUTS, 0.5, 0.77, '!4K!3[S]', charsi = 1.4

image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'kourganoff_constant.jpeg', image_jpeg, TRUE=1, QUALITY=100

;-------------------------------------------------------------------------------
; Linear source function:
;-------------------------------------------------------------------------------
a = 0.5
s = 1 + a * tau
j = (1. - 0.5*EXPINT(2, tau)) + a* (tau + 0.5*EXPINT(3, tau))
PLOT, tau, s, yr = [0, 2.], xr = [0, 1], /ys, /xs, $
  xtit = '!4s!3', ytickint = 0.5, xtickint = 0.5, charsi = 1.2
OPLOT, tau, j, lines = 2
XYOUTS, 0.24, 1.30, 'S = 1 + 0.5 !4s!3', charsi = 1.4
XYOUTS, 0.5, 1.00, '!4K!3[S]', charsi = 1.4

image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'kourganoff_linear1.jpeg', image_jpeg, TRUE=1, QUALITY=100

a = 1.5
s = 1 + a * tau
j = (1. - 0.5*EXPINT(2, tau)) + a* (tau + 0.5*EXPINT(3, tau))
PLOT, tau, s, yr = [0.5, 2.], xr = [0, 1], /ys, /xs, $
  xtit = '!4s!3', ytickint = 0.5, xtickint = 0.5, charsi = 1.2
OPLOT, tau, j, lines = 2
XYOUTS, 0.05, 1.50, 'S = 1 + 1.5 !4s!3', charsi = 1.4
XYOUTS, 0.1, 0.95, '!4K!3[S]', charsi = 1.4

image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'kourganoff_linear2.jpeg', image_jpeg, TRUE=1, QUALITY=100

a = 3.
s = 1 + a * tau
j = (1. - 0.5*EXPINT(2, tau)) + a* (tau + 0.5*EXPINT(3, tau))
PLOT, tau, s, yr = [0.5, 4.], xr = [0, 1], /ys, /xs, $
  xtit = '!4s!3', ytickint = 1, xtickint = 0.5, charsi = 1.2
OPLOT, tau, j, lines = 2
XYOUTS, 0.5, 2.2, 'S = 1 + 3 !4s!3', charsi = 1.4
XYOUTS, 0.1, 2., '!4K!3[S]', charsi = 1.4

image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'kourganoff_linear3.jpeg', image_jpeg, TRUE=1, QUALITY=100

-------------------------------------------------------------------------------
Exponential source function (Kourganoff, eq.14.37, p.46)
-------------------------------------------------------------------------------
a = 0.5
s = EXP(-a * tau)
j = (1./(2.*a)) * EXP(-a * tau) * (ALOG((a+1)/ABS(a-1)) - EXPINT(1, tau-a*tau)) + (1./(2.*a))*EXPINT(1, tau)
PLOT, tau, s, yr = [0, 1.2], xr = [0, 3], /ys, /xs, $
  xtit = '!4s!3', ytickint = 0.5, xtickint = 1., charsi = 1.2
OPLOT, tau, j, lines = 2
XYOUTS, 0.8, 0.7, 'S = EXP(-0.5 !4s!3)', charsi = 1.4
XYOUTS, 0.5, 0.5, '!4K!3[S]', charsi = 1.4

image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'kourganoff_exponential1.jpeg', image_jpeg, TRUE=1, QUALITY=100

a = 1.
s = EXP(-a * tau)
j = 0.5 * EXP(-tau) * (ALOG(2*tau) + 0.577215E0) + 0.5 * EXPINT(1, tau)
PLOT, tau, s, yr = [0, 1.2], xr = [0, 3], /ys, /xs, $
  xtit = '!4s!3', ytickint = 0.5, xtickint = 1., charsi = 1.2
OPLOT, tau, j, lines = 2
XYOUTS, 0.8, 0.5, 'S = EXP(-!4s!3)', charsi = 1.4
XYOUTS, 0.5, 0.35, '!4K!3[S]', charsi = 1.4

image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'kourganoff_exponential2.jpeg', image_jpeg, TRUE=1, QUALITY=100

; a = 5. 
; Note that in this case argument of the exponential function becomes negative.
; It requires an extended definition of En that is not covered by the IDL
; routine, see Kourganoff, Appendix I, eq. 36.9.

;-------------------------------------------------------------------------------
; Pulse source function:
;-------------------------------------------------------------------------------
tau1 = 1.
tau2 = 2.
int1 = WHERE(tau LT 1)
int2 = WHERE(tau GE 1 AND tau LT 2)
int3 = WHERE(tau GE 2)
s = tau * 0.
s[int2] = 1.
j = tau * 0.
j[int1] = 0.5 * (EXPINT(2, tau1-tau[int1]) - EXPINT(2, tau2-tau[int1]))
j[int2] = 1 - 0.5 * (EXPINT(2, tau[int2]-tau1) + EXPINT(2, tau2-tau[int2]))
j[int3] = 0.5 * (EXPINT(2, tau[int3]-tau2) - EXPINT(2, tau[int3]-tau1))
PLOT, tau, s, yr = [0, 1.5], xr = [0, 3], /ys, /xs, $
  xtit = '!4s!3', ytickint = 0.5, xtickint = 0.5, charsi = 1.2
OPLOT, tau, j, lines = 2
XYOUTS, 1.48, 1.02, 'S', charsi = 1.4
XYOUTS, 1.44, 0.71, '!4K!3[S]', charsi = 1.4

image_jpeg = TVRD(TRUE=1)
WRITE_JPEG, 'kourganoff_pulse.jpeg', image_jpeg, TRUE=1, QUALITY=100

END