;+
; NAME:
;       VARIATIONS3
;
; PURPOSE:
;
;       This function finds all the variations without repetition for the
;       elements of a given set. Order of variations is set to 3, though
;       the code can easily be adapted for any other order.
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
;       Radiative transfer, variability, statistics, arrays.
;
; CALLING SEQUENCE:
;
;       vars = VARIATIONS3(set)
;
; OPTIONAL INPUTS:
;
;       set = Array, integer. Set of integers. 
;
; OUTPUTS:
;
;       vars = Array, integer. Array that contains all the possible 
;              variations without repetition made out of the elements of
;              the input array.
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
;       vars = VARIATIONS3(INDGEN(5)+1)
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October 2012. 
;-
;================================================================================
; Carlson Quadrature by Nikola Vitas is licensed under a Creative Commons 
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

FUNCTION variations3, set

order = 3 ; number of dimensions, since we are in 3D space, we choose 3 directionals
          ; that define one direction.

setsize = SIZE(set)

IF setsize[0] NE 1 THEN MESSAGE, 'Set cannot have more than one dimension.'

n = N_ELEMENTS(set)

; Number of variations with repetition
nvar = n^order 

newsize = [2, nvar, order, setsize[2], nvar*order]

; Array variations contain all the variation with repetition
variations = MAKE_ARRAY(size = newsize)

FOR ix = 0, order-1 DO BEGIN
  nbl = n^(ix+1)             ; number of blocks
  nel = nvar/nbl             ; number of elements in each block
  FOR iy = 0, nbl-1 DO BEGIN ; now fill the blocks
    variations[iy*nel:(iy+1)*nel -1, ix] = set[iy MOD n]
  ENDFOR
ENDFOR

RETURN, variations
END