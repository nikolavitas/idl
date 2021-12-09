;+
; NAME:
;
;   PARTITION_COWLEY
;
; PURPOSE:
;
;   This procedure interpolates computes the partition function of a given
;   element and a given temperature using the pretabulated coefficients 
;   of Cowling. I adopted the data from a routine that I found in a file
;   codes.ada. I got the routine from A.Voegler back in 2007. I am not 
;   able to trace back this data. The only reference I found for Cowley's 
;   PF is in 2001A&A...374..265W and this is a reference to an ftp address
;   that does not exist any more in 2013.
;      Cowley, C. R., 1998,  
;          ftp://astro.lsa.umich.edu/pub/get/cowley/partition/bolpfn.f
;   If you use this routine for a publication, please try to trace the 
;   data yourself to find a proper reference.
;
;   The range of the PF validity is also not known to me. From the original
;   routine I concluded that the upper temperature limit is 31 000 K, so I 
;   added that number to the test range.
; 
;   Note that this is the only set (of those in this library) that has data 
;   for the fourth stage of ionization of some elements
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
;   pf = PARTITION_COWLEY(t, z)
;
; INPUTS:
;
;   t =      Scalar or array, float. Temperature (K)
;
;   z =      Scalar. Atomic number (1 = H, 2 = He, ...)
;
; OUTPUTS:
; 
;   pf =  Structure with 3 tags, u1, u2 and u3 (for the partition functions
;         of the first three ionization stages (neutral, + and ++).
;
; KEYWORDS:
;
;  range     = If it is set and the temperature is out of the range, the value 
;              of pf at the edge of the range adopted.;
;
; COMMENT:
;
; DEPENDENCIES:
;
;
; EXAMPLE:
;
;   pf = PARTITION_COWLEY([6E3, 7E3, 1E4], 7)
;   PRINT, pf.u1
;   PRINT, pf.u2
;   PRINT, pf.u3
;   PRINT, pf.u4
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (October2012)
;-
;================================================================================
; PARTITION_COWLEY, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; For the data in pf_cowley.dat, please see the remark in the header of this
; file.
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

FUNCTION partition_cowley, t, z, range = range

trange = [1.D3, 3.1D4] ; K
tarray = t
zeroarray = t*0.D0
IF KEYWORD_SET(range) THEN tarray = (tarray > MIN(trange)) < MAX(trange) 

theta = ALOG10(5040.D0/tarray)
;-------------------------------------------------------------------------------
; Read the data
;-------------------------------------------------------------------------------
@pf_path

IF NOT(FILE_TEST(path+'pf_cowley.sav')) THEN BEGIN
  line = ''
  nlines = 92*4
  pfcoeffs = DBLARR(92, 4, 6)
  coeffs = DBLARR(6)
  zz = 0
  rr = 0
  ion1 = 0
  OPENR, un, path+'pf_cowley.dat', /get_lun
  FOR il = 0, 3 DO READF, un, line
  FOR iz = 0, 91 DO BEGIN
    FOR ir = 0, 3 DO BEGIN
      READF, un, zz, rr, coeffs, ion1
      pfcoeffs[iz, ir, *] = coeffs
    ENDFOR  
  ENDFOR
  FREE_LUN, un
  SAVE, pfcoeffs, file = path+'pf_cowley.sav'
ENDIF ELSE BEGIN
  RESTORE, path+'pf_cowley.sav'
ENDELSE

;u = DBLARR(4, N_ELEMENTS(theta))
u1 = zeroarray
u2 = zeroarray
u3 = zeroarray
u4 = zeroarray
utemp = zeroarray

; FOR ir = 0, 3 DO BEGIN
; ;   PRINT, pfcoeffs[z-1, ir, *]
;   IF pfcoeffs[z-1, ir, 1] LE -23. THEN BEGIN
;     u[ir, *] = -1.
;   ENDIF ELSE BEGIN
;     u[ir, *] = pfcoeffs[z-1, ir, 0]
;     IF TOTAL(pfcoeffs[z-1, ir, 1:5]) NE 0 THEN BEGIN
;       u[ir, *] = u[ir, *] + EXP(pfcoeffs[z-1, ir, 1] + theta * (pfcoeffs[z-1, ir, 2] + $
;                   theta * (pfcoeffs[z-1, ir, 3] + theta * (pfcoeffs[z-1, ir, 4] + $
;                   theta * (pfcoeffs[z-1, ir, 5])))))
;     ENDIF
;   ENDELSE
; ENDFOR  

FOR ir = 0, 3 DO BEGIN
  IF pfcoeffs[z-1, ir, 1] LE -23. THEN BEGIN
    utemp = -1.
  ENDIF ELSE BEGIN
    utemp = pfcoeffs[z-1, ir, 0]
    IF TOTAL(pfcoeffs[z-1, ir, 1:5]) NE 0 THEN BEGIN
      utemp = utemp + EXP(pfcoeffs[z-1, ir, 1] + theta * (pfcoeffs[z-1, ir, 2] + $
                  theta * (pfcoeffs[z-1, ir, 3] + theta * (pfcoeffs[z-1, ir, 4] + $
                  theta * (pfcoeffs[z-1, ir, 5])))))
    ENDIF
  ENDELSE
  
  CASE ir+1 OF
    1 : u1 = utemp
    2 : u2 = utemp
    3 : u3 = utemp
    4 : u4 = utemp
  ENDCASE
ENDFOR  

result = {u1:u1, u2:u2, u3:u3, u4:u4}

RETURN, result
END  
