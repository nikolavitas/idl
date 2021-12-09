;+
; NAME:
;
;   PARTITION_SAUVAL_AND_TATUM
;
; PURPOSE:
;
;   This procedure interpolates computes the partition function of a given
;   element and a given temperature using the pretabulated coefficients 
;   of Sauval and Tatum (1984, ApJ Suppl, 56, 193; 1984ApjS...56..193S)
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
;   Atomic data.
;
; CALLING SEQUENCE:
;
;   pf = PARTITION_SAUVAL_AND_TATUM(t, z)
;
; INPUTS:
;
;   t =      Scalar or array, float. Temperature (K)
;
;   z =      Scalar. Atomic number (1 = H, 2 = He, ...)
;
; OUTPUTS:
; 
;   pf =  Structure with 2 tags, u1 and u2 (for the partition functions
;         of the first two ionization stages.
;
; KEYWORDS:
;
;  range     = If it is set and the temperature is out of the range, the value 
;              of pf at the edge of the range adopted.
;
;   missing  = The value that is returned if the data is  missing in the table 
;            (note that the default value is -1).
;
; DEPENDENCIES:
;
;
; EXAMPLE:
;
;   pf = PARTITION_SAUVAL_AND_TATUM([6E3, 7E3, 1E4], 7)
;   PRINT, pf.u1
;   PRINT, pf.u2
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (July2013)
;-
;================================================================================
; PARTITION_WITTMANN, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in this routine are taken from Sauval and Tatum (1984, ApJ Suppl, 56, 
; 193). If you use the code and the data please acknowledge it by citing the 
; original paper.
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

FUNCTION partition_sauval_and_tatum, t, z, range = range, missing = missing

@pf_path
filename = path+'pf_sauval_and_tatum.dat'

IF NOT(KEYWORD_SET(missing)) THEN missing = -1

tarray = t
zeroarray = t*0.D0

trange = [1.D3, 9.D3] ; K
IF KEYWORD_SET(range) THEN tarray = (tarray > MIN(trange)) < MAX(trange)

IF NOT(FILE_TEST(path+'pf_sauval_and_tatum.sav')) THEN BEGIN
  coeffs = DBLARR(92, 2, 5)
  OPENR, u, filename, /get_lun
  num = FILE_LINES(filename)-1
  data = FLTARR(8, num)
  a = ''
  READF, u, a
  dataline = ''
  FOR i = 0, num-1 DO BEGIN
    READF, u, dataline
    sa = STRSPLIT(dataline, ' ', /extract, count = nc)
    atom = FIX(sa[0])
    ion =  FIX(sa[1])
    coeffs[atom-1, ion, 0:nc-4] = DOUBLE(sa[3:nc-1])
  ENDFOR
  FREE_LUN, u
  SAVE, coeffs, file = path+'pf_sauval_and_tatum.sav'
ENDIF ELSE RESTORE, path+'pf_sauval_and_tatum.sav'
  
; tt = REFORM(t)
; tt = [tt]
; nt = N_ELEMENTS(tt)
tsize = SIZE(t)
theta = 5040.D0/tarray

result = zeroarray

theta = 5040.D0/tarray

pf = {u1:zeroarray, u2:zeroarray}

FOR r = 0, 1 DO BEGIN
  result = zeroarray
  cfs = REFORM(coeffs[z-1, r, *])
  n = FIX(TOTAL(cfs NE 0))
  IF n EQ 0 THEN BEGIN
    result = result + missing
  ENDIF ELSE BEGIN
    FOR i = 0, n-1 DO BEGIN
      result = result + cfs[i]*(ALOG10(theta))^i
    ENDFOR
    result = 10.D0^result
  ENDELSE
  IF z EQ 1 AND r EQ 1 THEN result = zeroarray + 1.
  IF z EQ 2 AND r EQ 0 THEN result = zeroarray + 1.
  IF z EQ 3 AND r EQ 1 THEN result = zeroarray + 1.
  
  IF r EQ 0 THEN pf.u1 = result
  IF r EQ 1 THEN pf.u2 = result  

ENDFOR  
  


RETURN, pf

END