;+
; NAME:
;       CARLSON_QUADRATURE_B
;
; PURPOSE:
;
;       This procedure finds the directional cosine and their weights for
;       an angular quadrature described in 
;
;         Carlson B.G.,1963, in Alder B., Fernbach S.(eds.) Methods in 
;         Computational Physics. Vol. 1 p.1
;
;       The problem is solved for the so-called Set B that is an alternative
;       to the commonly used Set A. For more details, please see the paper
;       of Carlson, and an older paper by Carlson and Bell (1958). 
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
;       CARLSON_QUADRATURE_B, n, wl=wl, mu=mu, wm=wm, variations=variations
;
; INPUTS:
;
;       n = Scalar, integer. This is the number of mu-levels. The number of rays
;           is (n+2)(n+4)/8-3.
;
; OUTPUTS:
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
;       wl         = Array, float. Weights of directional cosines.
;       wm         = Array, float. Weights of rays.
;       mu         = Array, float. Directional cosines.
;       variations = Array, integer. Possible variations of directional 
;                    cosines
;
;  EXAMPLE:
;
;       CARLSON_QUADRATURE_B, 8, wl=wl, mu=mu, wm=wm, variations=variations
;
; DEPENDENCIES:
;
;       Requires VARIATIONS3, MU1NR_B and W1NR_B.
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, November 2012. Based on Carlson (1963).
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


PRO carlson_quadrature_b, n, wl=wl, mu=mu, wm=wm, variations=variations

nprim = n+2

lmax = n/2

; The set B for n is constructed using set A for n + 2 (by cutting the
; corners, Carlson, 1963, p.15).
lmaxprim = nprim/2

values = INDGEN(lmaxprim)+1 ; This is the set of all possible values

; We define sum of indices for one direction
IF FLOAT(lmaxprim)/2 EQ FLOOR(lmaxprim/2) THEN suma = (2*lmaxprim-4)/2.+4 $
                                          ELSE suma = (2*lmaxprim-2)/2.+3

; Find the suitabel variations:
variations = VARIATIONS3(values)

; Extract a subset that has a given value as a sum.
totvar = TOTAL(variations, 2)      ; sum of each variation
maxvar = MAX(variations, dim = 2)  ; max of each variation
varind = WHERE(totvar EQ suma AND maxvar LT lmaxprim) ; those with correct sum 
                                                      ; minus the corners 

variations = variations[varind, *] ; we keep only the ones that we need
nvar = (SIZE(variations))[1]

; Compute weight that correspond to directionals
omega = W1NR_B(n)

; Compute directional angles
mu = MU1NR_B(omega, n)

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

