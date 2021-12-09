;+
; NAME:
;
;   PARTITION_GRAY
;
; PURPOSE:
;
;   This procedure interpolates the partition function table from Gray (The 
;   observations and analysis of stellar photospheres, CUP, 2005) for a given
;   temperature T and for a element (Z). The data is available only for the 
;   first three ionization stages (not for all elements.). 
;
;   Gray: "The actual computation of partition function can be somewhat 
;   laboriuos and requires a detailed knowledge of the energy levels. For model
;   photosphere computation, interpolation within (this) Table is easy and 
;   convinient. In many publications. polynomials are given, and while these 
;   are also convinient, they can be misleading when used outside the temperature 
;   range for which they were intended. (...) The entries (in the Table) are 
;   accurate to $\approx 1\%$, and come from several sources: Aller (1963), Evans 
;   (1966), Bolton (1970), Irwin (1981), Cowley and Adelman (1983), Sauval and 
;   Tatum (1984), Milone and Merlo (1998), Halenka et al. (2001), and NIST atomic 
;   data base. Sauval and Tatum (1984) also give polynomials for partition 
;   functions for 300 diatomic molecules."
;
;   The interpolation is linear in the logarithmic scale.
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
;   pf = PARTITION_GRAY(t, z)
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
;         of the first three ionization stages.
;
; KEYWORDS:
;
;   log      = If set, the result is log u(T). If not, it is just u(T).
;
;   missing  = The value that is returned if the data is  missing in the table 
;            (note that the default value is 2). .
;
;  range     = If it is set and the temperature is out of the range, the value 
;              of pf at the edge of the range adopted.
; 
;
; COMMENT:
;
;   Use of the table is restricted to the elements and stages of ionisation 
;   listed in Gray (2005). Table file pf_gray.dat that has to be in 
;   the same folder as this procedure.
;
; DEPENDENCIES:
;
;
; EXAMPLE:
;
;   pf = PARTITION_FUNCTION([6E3, 7E3, 1E4], 7)
;   PRINT, pf.u1
;   PRINT, pf.u2
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (Jan2006)
;               # as a function, NV (Nov2007)
;               # output for non-existing atom NV (Dec2008)
;               # r removed from the input NV (Jul2013)
;
;-
;================================================================================
; PARTITION_GRAY, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in pf_gray.dat are taken from Table D.2 (p.514) of "The observations
; and analysis of stellar photospheres", David F. Gray, 3rd edition, 2005. If 
; you use this data please cite that book as the reference.
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

PRO prep_gray, t, z, log = log, missing = missing, range = range









@pf_path
filename = path + 'pf_gray.dummy'
IF NOT(KEYWORD_SET(missing)) THEN missing = 2

GET_LUN, u
OPENR, u, filename
num = FILE_LINES(filename)
data = FLTARR(13, 247)
READF, u, data
FREE_LUN, u

OPENW, 1, 'dummy.txt'
FOR j = 0, 246 DO BEGIN
ll = ''
FOR i = 0, 12 DO BEGIN
  ll = ll + ' ' + STRING(data[i, j], format = '(F6.3)') 
  IF i LT 12 THEN ll = ll + ',' ELSE ll = ll + ', &'
ENDFOR
PRINTF, 1, ll 
ENDFOR

CLOSE, 1
  

END
