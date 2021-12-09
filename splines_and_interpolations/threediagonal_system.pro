;+
; NAME:
;       THREEDIAGONAL_SYSTEM
;
; PURPOSE:
;
;       This function solves threediagonal system of n linear equations. The
;       system is specified by its diagonal elements. It is solved using 
;       Gaussian elimination. 
;
;       The algorithm is described in details at
;       http://nikolavitas.blogspot.com.es/2013/08/three-diagonal-system-of-linear.html
;
;       Note: It can be coded without the explicit definition of the prim
;       variables (for the efficiency reasons). However, it seems to me
;       that the code is a bit more understandable when the prim variables
;       are explicitly introduced. The main purpose of this blog is to explain 
;       how things work and to provide an understandable and easy-to-adapt 
;       code (and not necessarily the most efficient solutiton - IDL already 
;       has many of them integrated in the intrinsic functions).
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
;       a  = Array, float or double. The lower diagonal. a[0] must = 0.
;       b  = Array, float or double. The main diagonal. 
;       c  = Array, float or double. The upper diagonal. c[n-1] must = 0.
;       d  = Array, float or double. The right-hand side.
; 
; OUTPUTS:
;
;       x  = Array, float or double. The result of the system.
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
; 
;       ; We form a system of five equations
;       aa = [0., 2.,  3., 4., 1.]
;       bb = [3., 4., 11., 7., 2.]
;       cc = [1., 1.,  1., 3., 0.]
;       
;       ;       Set the exact solution
;       x0 = FINDGEN(5)
; 
;       ; Compute the right-hand side to complete the system
;       dd = FLTARR(5)
;       dd[0] = bb[0]*x0[0] + cc[0]*x0[1] + aa[0]*x0[4]
;       FOR i = 1, 3 DO dd[i] = aa[i]*x0[i-1] + bb[i]*x0[i] + cc[i]*x0[i+1]
;       dd[4] = aa[4]*x0[3] + bb[4]*x0[4] + cc[4]*x0[0]
;       
;       ; Solve the system
;       x = threediagonal_system(aa, bb, cc, dd)
; 
;       ; Compare the solution of the system and the exact solution
;       PRINT, x
;       PRINT, x0
; 
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, August 2013.
;-
;================================================================================
; THREEDIAGONAL_SYSTEM, IDL routine by Nikola Vitas is licensed under a Creative 
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

FUNCTION threediagonal_system, a, b, c, d

n = N_ELEMENTS(d)
x = MAKE_ARRAY(size = size(d))

dprim = d
bprim = b

; Forth
FOR i = 1, n-1 DO BEGIN
  dprim[i] = d[i] - (a[i]*dprim[i-1])/bprim[i-1]
  bprim[i] = b[i] - (a[i]*c[i-1])/bprim[i-1]
ENDFOR

; Back
x[n-1] = dprim[n-1]/bprim[n-1]
FOR i = n-2, 0, -1 DO x[i] = (dprim[i] - c[i]*x[i+1])/bprim[i]

RETURN, x
END