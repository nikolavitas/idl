;+
; NAME:
;       CARLSON_QUADRATURE_A
;
; PURPOSE:
;
;       This procedure finds the directional cosine and their weights for
;       an angular quadrature described in 
;
;         Carlson B.G.,1963, in Alder B., Fernbach S.(eds.) Methods in 
;         Computational Physics. Vol. 1 p.1
;
;       The problem is solved for the so-called Set A that is frequently
;       used for discretization of the radiative transfer equation (see
;       Bruls et al, 1999). The equation for W1 coefficient is solved by
;       Newton-Raphson method and the set of linear equations for the 
;       Wm weights (weights per ray) is solved using built-in Cramer 
;       method. For given n the code finds out possible variations of 
;       directional cosines and use them for the weight evaluation. See
;       CARLSON_QUADRATURE_B.PRO for set B.
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
;       CARLSON_QUADRATURE_A, n, wl=wl, mu=mu, wm=wm, variations=variations
;
; INPUTS:
;
;       n = Scalar, integer. This is the number of mu-levels. The number of rays
;           is n(n+2)/8.
;
; OUTPUTS:
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
;       wl         = Array(n/2), float. Weights of directional cosines.
;       wm         = Array(n/2-1), float. Weights of rays.
;       mu         = Array(n/2), float. Directional cosines.
;       variations = Array(n(n+2)/8, 3), integer. Possible variations of directional 
;                    cosines
;
;  EXAMPLE:
;
;       CARLSON_QUADRATURE_A, 8, wl=wl, mu=mu, wm=wm, variations=variations
;
; DEPENDENCIES:
;
;       Requires VARIATIONS3 and W1NR_A.
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October 2012. Based on Bruls et al (1999) and 
;         Carlson (1963).
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

PRO carlson_quadrature_a, n, wl=wl, mu=mu, wm=wm, variations=variations

lmax = n/2

values = INDGEN(lmax)+1 ; This is the set of all possible values

; We define sum of indices for one direction, for given n (we still assume only set A)
IF FLOAT(lmax)/2 EQ FLOOR(lmax/2) THEN suma = (2*lmax-4)/2.+4 ELSE suma = (2*lmax-2)/2.+3

; Find the suitabel variations:
variations = VARIATIONS3(values)

; Extract a subset that has a given value as a sum.
totvar = TOTAL(variations, 2)      ; sum of each variation
varind = WHERE(totvar EQ suma)     ; those with correct sum
variations = variations[varind, *] ; we keep only the ones that we need
nvar = (SIZE(variations))[1]

; Now we need to find coefficients w1 and mu1 fpr the given value of n.
delta = 2.D/(n-1.D)

; Compute directional angles
mu = DBLARR(lmax)
mu[0] = SQRT(1.D/(3.D*(n-1.D)))
FOR l = 2, lmax DO $
  mu[l-1] = SQRT(mu[0]^2 + (l-1)*delta)

; Compute weight that correspond to deirectionals
omega = W1NR_A(n)
FOR l = 0, lmax-1 DO PRINT, omega[l], mu[l]

; Identify classes
varinclass = INTARR(nvar)
FOR iv = 0, nvar-1 DO BEGIN
  set = REFORM(variations[iv, *])
  set = set[SORT(set)]
  set = STRJOIN(STRCOMPRESS(STRING(set), /rem))
  IF iv EQ 0 THEN BEGIN
    classes = [set] 
  ENDIF ELSE BEGIN
    test = WHERE(set EQ classes)
    IF test EQ -1 THEN BEGIN
      classes = [classes, set]
      varinclass[iv] = N_ELEMENTS(classes)-1
    ENDIF ELSE BEGIN
      varinclass[iv] = test
    ENDELSE
  ENDELSE
ENDFOR
nclasses = N_ELEMENTS(classes) ; must be = lmax-1

; At the end we still have to find the weights for particular rays. 
; The matrix of the system we construct by desceding through an octant,
; see Fig.2 of Carlson (1963). 
mtrx = FLTARR(lmax-1, lmax-1)
wl = REVERSE(omega[1:lmax-1])
FOR im = 0, lmax-2 DO BEGIN 
  melem = varinclass(WHERE(REFORM(variations[*, 0]) EQ lmax-im))
  nm = N_ELEMENTS(melem)
;   PRINT, melem
  FOR ie = 0, nm-1 DO mtrx[melem[ie], im] = mtrx[melem[ie], im]+1
ENDFOR
wm = CRAMER(mtrx, wl, /double)

; Prepare the output
wm = REVERSE(wm)
wl = omega

END

