;+
; NAME:
;       LOAD_ELECTRON_AFFINITY
;
; PURPOSE:
;
;       This function loads the values of the electron affinity in [eV] for the 
;       92 atomic elements
;
;       The data is downloaded using the interface at:
;       https://en.wikipedia.org/wiki/Electron_affinity_(data_page)
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
;       Atomic data.
;
; CALLING SEQUENCE:
;
;       chimin = LOAD_ELECTRON_AFFINITY(z)
;
; INPUTS:
;
;       z        = String or array, integer. The atomic number(s). If not 
;                  specified, the abundances for the first 92 elements are
;                  returned. It must not be specified if the PPM keyword is set.
; 
; OUTPUTS:
;
;        eaffi  = String or 1D or 2D array, double. The electron affinity in eV.
;
; INPUT KEYWORDS:
; 
; OUTPUT KEYWORDS:
;
; EXAMPLE:
; 
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October 2012.
;       - Data update (NV), October 2019.
;-
;================================================================================
; LOAD_ELECTRON_AFFINITY, IDL routine by Nikola Vitas is licensed under a 
; Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; This routine contains the data from several sources listed above. If any of 
; that data is used in a publication, please cite the corresponding reference.
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

FUNCTION load_electron_affinity, z

COMPILE_OPT idl2, HIDDEN

IF KEYWORD_SET(z) THEN BEGIN
  IF MAX(z) GT 92 OR MIN(z) LT 1 THEN MESSAGE, 'Z out of range. STOP.'
ENDIF

eaffi = [   0.754, -0.500,  0.618, -0.500,  0.279,  1.262, -0.070,  1.461,  3.401, -1.200, $
            0.548, -0.400,  0.432,  1.389,  0.747,  2.077,  3.612, -1.000,  0.501,  0.024, $
            0.188,  0.075,  0.528,  0.676, -0.500,  0.153,  0.662,  1.157,  1.235, -0.600, $
            0.430,  1.232,  0.804,  2.020,  3.364, -1.000,  0.486,  0.052,  0.307,  0.433, $
            0.917,  0.747,  0.550,  1.045,  1.142,  0.562,  1.304, -0.700,  0.384,  1.112, $
            1.047,  1.971,  3.059, -0.800,  0.472,  0.144,  0.550,  0.570,  0.962,  1.916, $
            0.129,  0.162,  0.116,  0.137,  1.165,  0.352,  0.338,  0.312,  1.029, -0.020, $
            0.239,  0.178,  0.323,  0.816,  0.060,  1.078,  1.564,  2.125,  2.309, -0.500, $
            0.377,  0.357,  0.942,  1.400,  2.420, -0.700,  0.486,  0.100,  0.350,  1.170, $
            0.550,  0.530 ]
            
    
RETURN, DOUBLE(eaffi[z-1])

END
