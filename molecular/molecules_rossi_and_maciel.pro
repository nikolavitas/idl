;+
; NAME:
;
;   MOLECULES_ROSSI_AND_MACIEL
;
; PURPOSE:
;
;   This routine evaluates the interpolation polynomial the chemical 
;   equilibrium constant for a given temperature and molecular specie. 
;   The coefficients of the polynomial are from Rossi and Maciel (1983).
;   For more information, check the header of diatomics.pro. This 
;   routine includes multi-nuclear molecules and corresponds to 
;   Chem_Equi version 5 (April 2019). In this version the data from
;   particular surveys is included in the catalogue file. 
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
;   f = MOLECULES_SAUVAL_AND_MACIEL(t, name, type = type, d0 = d0, mab = mab)
;
; INPUTS:
;
;   t =      Scalar or array, float. Temperature (K)
;
;   name =   Scalar string. Name = chemical formula of the molecule
;
; OUTPUTS:
; 
;   f  =     Array, float. The output function (the partition function,
;            the chemical equilibrium constant or the internal energy).
;
; KEYWORDS:
;
;   type = String. Specifies the output function ('pf', the partition function,
;          default; 'kp', the chemical equilibrium constant for pressure;
;          or 'eint', the internal energy).
;
;   d0   = Scalar, float. The dissociation constant in eV.
;
;   mab  = Scalar, float. The reduced mass, mab = ma*mb/(ma+mb) in a.m.u.
;   
; COMMENT:
;
; DEPENDENCIES:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (February 2014)
;
;-
;================================================================================
; MOLECULES_ROSSI_AND_MACIEL, IDL routine by Nikola Vitas is licensed under a  
; Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in sauvalandtatum1984_*.sav are not integral part of this routine. If  
; you use the data in a publication, please acknowledge it by citing the proper  
; reference (1984ApJS...56..193S).
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

FUNCTION molecules_rossi_and_maciel, t, name, type = type, d0 = d0, mab = mab

IF NOT KEYWORD_SET(type) THEN type = 'kp'

kkev    = 8.6173324D-5         ; eV K^-1                  (from NIST)
theta = ALOG10(EXP(1.0))/(kkev * t)

@molecules_path

IF type EQ 'kp' THEN BEGIN

  RESTORE, path + 'catalogue_of_molecules.sav'
  data = rm_eqc
  index = WHERE(name EQ data.name)

  kp = t*0.
  coeffs = data[index].coeffs
  nc = N_ELEMENTS(coeffs)

  FOR ii = 0, nc-1 DO $
      kp = kp + coeffs[ii]*theta^ii
  kp = 10.^kp

  result = kp

ENDIF

IF type EQ 'pf' OR type EQ 'eint' THEN BEGIN
  PRINT, 'Not available.'
  result = -1
ENDIF

RETURN, result
END
