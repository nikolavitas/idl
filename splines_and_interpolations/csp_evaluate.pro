;+
; NAME:
;
;       CSP_EVALUATE
;
; PURPOSE:
;
;       This function takes the coefficients of the periodic cubic splines
;       as an input and returns the values of the interpolation function
;       for a given array of x.
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
;       x0  =  Array, float. Array of the original x-values at which the unknown 
;             function is sampled.
;
;       x   = Array, float. Values of x at which the interpolation function will be
;             evaluated.
;
;       cp  = Structure. It contains the coefficients of the cubic splines for 
;             each segment defined by x0. These coefficients are assumed to be computed 
;             in advance using CSP_INTERPOLATION function.
;
; OUTPUTS:
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
; 
;      ; See the example in CSP_INTERPOLATION
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October, 2013.
;-
;================================================================================
; CSP_EVALUATE, IDL routine by Nikola Vitas is licensed under a Creative 
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

FUNCTION csp_evaluate, x0, x, cp

nx = N_ELEMENTS(x)
nx0 = N_ELEMENTS(x0)
period = x[nx-1]-x[0]
y0 = DBLARR(nx0)
FOR i = 0, nx0-1 DO BEGIN
  factor = FLOOR((x0[i]-x[0])/period)
  j = VALUE_LOCATE(x, x0[i]-factor*period) 
  xtemp = x0[i]-factor*period
  y0[i] = cp.a[j] + cp.b[j]*(xtemp-x[j]) + cp.c[j]*(xtemp-x[j])^2 + cp.d[j]*(xtemp-x[j])^3
ENDFOR

RETURN, y0
END