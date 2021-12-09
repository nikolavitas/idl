;+
; NAME:
;
;   PARTITION_IRWIN
;
; PURPOSE:
;
;   This procedure interpolates computes the partition function of a given
;   element and a given temperature using the pretabulated coefficients 
;   of Alan W. Irwin (see Irwin, 1981, ApJ Suppl, 45, 621). The coefficients
;   used here are updated values read from the file IRWIN_atoms_v07.3.dat
;   (see the file for the further information and reference).
;
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
;   pf = PARTITION_IRWIN(t, z)
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
;   pf = PARTITION_IRWIN([6E3, 7E3, 1E4], 7)
;   PRINT, pf.u1
;   PRINT, pf.u2
;   PRINT, pf.u3
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (August2013)
;-
;================================================================================
; PARTITION_IRWIN, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in IRWIN_atoms_v07.3.dat are not integral part of this routine. If you 
; use the data please contact Alan W. Irwin (irwin@uvastro.phys.uvic.ca) and cite
; the proper reference (see the header above).
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

FUNCTION partition_irwin, t, z, range = range

trange = [1.D3, 1.6D4] ; K
IF KEYWORD_SET(range) THEN tt = (t > MIN(trange)) < MAX(trange) ELSE tt = t

;-------------------------------------------------------------------------------
; Read the data
;-------------------------------------------------------------------------------
@pf_path

IF NOT(FILE_TEST(path+'pf_irwin.sav')) THEN BEGIN
  line = ''
  pfcoeffs = DBLARR(92, 3, 6)
  coeffs = DBLARR(6)
  OPENR, un, path+'IRWIN_atoms_v07.3.dat', /get_lun
  FOR il = 0, 27 DO READF, un, line
  FOR il = 0, 1 DO BEGIN
    READF, un, line
    READS, STRMID(line, 56, strlen(line)-56), coeffs
    pfcoeffs[0, il, *] = coeffs
  ENDFOR  

  FOR iz = 1, 91 DO BEGIN
    FOR il = 0, 2 DO BEGIN
      READF, un, line
      READS, STRMID(line, 56, STRLEN(line)-56), coeffs
      pfcoeffs[iz, il, *] = coeffs
    ENDFOR
  ENDFOR
  FREE_LUN, un
  SAVE, pfcoeffs, file = path+'pf_irwin.sav'
ENDIF ELSE BEGIN
  RESTORE, path+'pf_irwin.sav'
ENDELSE
  
;-------------------------------------------------------------------------------
; Compute the partition function
;-------------------------------------------------------------------------------

nt = N_ELEMENTS(tt)
zz = ALOG(tt)

IF z EQ 1 THEN BEGIN

  u1 = POLY(zz, pfcoeffs[0, 0, *])
;   pos = WHERE(tt GT 16000.)
;   IF (pos[0] NE -1) THEN u1[pos] = EXP(7.19420668D-1)
  u2 = DBLARR(nt) + 1.
  u3 = DBLARR(nt)
ENDIF ELSE BEGIN
  u1 = POLY(zz, pfcoeffs[z-1, 0, *])
  u2 = POLY(zz, pfcoeffs[z-1, 1, *])
  u3 = POLY(zz, pfcoeffs[z-1, 2, *])
  u2 = EXP(u2)
  u3 = EXP(u3)

ENDELSE
u1 = EXP(u1)

IF z EQ 1 THEN u3 = 0.*u3 - 1.
IF z EQ 2 THEN u3 = 0.*u3 + 1.

result = {u1:u1, u2:u2, u3:u3}

RETURN, result
          

END
