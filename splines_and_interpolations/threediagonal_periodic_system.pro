;+
; NAME:
;       THREEDIAGONAL_PERIODIC_SYSTEM
;
; PURPOSE:
;
;       This function solves three-diagonal system of n periodic linear 
;       equations. The difference in respect to the ordinary three-diagonal
;       system is that the corners of the matrix of the system are now all
;       filled (so the coefficients are cyclic).
;
;       The algorithm of the solution is described in details at
;       http://nikolavitas.blogspot.com.es/2013/08/three-diagonal-system-of-linear.html
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
;       a  = Array, float or double. The lower diagonal. 
;       b  = Array, float or double. The main diagonal. 
;       c  = Array, float or double. The upper diagonal.
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
;       aa = [6., 2.,  3., 4., 1.]
;       bb = [3., 4., 11., 7., 2.]
;       cc = [1., 1.,  1., 3., 3.]
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
;       x = threediagonal_periodic_system(aa, bb, cc, dd)
; 
;       ; Compare the solution of the system and the exact solution
;       PRINT, x
;       PRINT, x0
; 
; DEPENDENCIES:
;
;       THREEDIAGONAL_SYSTEM
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, August 2013.
;-
;================================================================================
; THREEDIAGONAL_PERIODIC_SYSTEM, IDL routine by Nikola Vitas is licensed under a 
; Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
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

FUNCTION threediagonal_periodic_system, a, b, c, d

n = N_ELEMENTS(d)
x = MAKE_ARRAY(size = size(d))

; First threediagonal subproblem, x1 solution
a1 = [0, a[1:n-2]]
b1 = b[0:n-2]
c1 = [c[0:n-3], 0]
d1 = d[0:n-2]
x1 = threediagonal_system(a1, b1, c1, d1)

; Second threediagonal subproblem, x2 solution
d2 = d1
d2[0] = -a[0]
d2[1:n-3] = 0.
d2[n-2] = -c[n-1]
x2 = threediagonal_system(a1, b1, c1, d2)

; Evaluate the last x-value
x[n-1] = (d[n-1] - c[n-1]*x1[0] - a[n-1]*x1[n-2])/(b[n-1] + c[n-1]*x2[0] + a[n-1]*x2[n-2])

; Evaluate other x
x[0:n-2] = x1 + x2*x[n-1]

RETURN, x
END