;+
; NAME:
;
;   PF
;
; PURPOSE:
;
;   This is a wrapper procedure for computing the atomic partition function
;   from precomputed and tabulated parameters of the best fit. The routine 
;   takes temperature T and atomic number Z as input parameters. There are
;   four available datasets and four functions that use them:
;
;   Function:                     Reference:             Data:
;   PARTITION_GRAY                2005oasp.book.....G    pf_gray.dat
;   PARTITION_IRWIN               1981ApJS...45..621I    pf_irwin.sav
;   PARTITION_SAUVAL_AND_TATUM    1984ApjS...56..193S    pf_sauval_and_tatum.dat
;   PARTITION_WITTMANN            1974SoPh...35...11W    hardcoded
;   PARTITION_COWLEY         see: 2001A&A...374..265W    pf_cowley.dat
;
;   If you use any of these data, please acknowledge it by citing the 
;   corresponding paper.
;
;   Note that each of these datasources strictly applies only to a limited
;   range of temperatures: Gray (2520 - 25200 K), Irwin (1000 - 16000 K), 
;   Sauval and Tatum (1000 - 9000 K). Out of these ranges the values returned 
;   by this function may not be reliable and it is up to the user to decide
;   on the extrapolation procedure in that case. If keyword RANGE is set in
;   a call to a partition function, then the temperatures out of the range
;   are clipped to it.
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
;   pf = PF(t, z, data = data)
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
;         of the first three ionization stages(only two tags, u1 and u2 if the
;         data of Sauval and Tatum is used).
;
; KEYWORDS:
;
;   data = String, scalar. Reference to the data that has to be used. 
;          Possibilities are: irwin, wittmann, gray, sauval_and_tatum, cowley.
;          Only the first letter matters!
;
; COMMENT:
;
; DEPENDENCIES:
;
;
; EXAMPLE:
;
;   t = DINDGEN(5)*100+5000.
;   z = 3 ; Litium
;   pfw = PF(t, z, data = 'w')
;   pfi = PF(t, z, data = 'i')
;   pfg = PF(t, z, data = 'g')
;   pfs = PF(t, z, data = 's')
;   pfc = PF(t, z, data = 'c')
;   
;   PRINT, 'Gray:'
;   PRINT, pfg.u1
;   PRINT, pfg.u2
;   PRINT, pfg.u3
;   
;   PRINT, 'Wittmann:'
;   PRINT, pfw.u1
;   PRINT, pfw.u2
;   PRINT, pfw.u3
;   
;   PRINT, 'Irwin:'
;   PRINT, pfi.u1
;   PRINT, pfi.u2
;   PRINT, pfi.u3
;   
;   PRINT, 'Sauval and Tatum:'
;   PRINT, pfs.u1
;   PRINT, pfs.u2
;
;   PRINT, 'Cowley'
;   PRINT, pfc.u1
;   PRINT, pfc.u2
;   PRINT, pfc.u3
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (October 2013)
;     # All subroutines adapted to work with 1D/2D/3D temperature arrays,
;       (Dec, 2015)
;
;-
;================================================================================
; PF, IDL routine by Nikola Vitas is licensed under a Creative Commons 
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

FUNCTION pf, t, z, data = data

letter = STRLOWCASE(STRMID(data, 0, 1))

CASE letter OF 
  'i' : pf = PARTITION_IRWIN(t, z, /range)
  'w' : pf = PARTITION_WITTMANN(t, z)
  'g' : pf = PARTITION_GRAY(t, z)
  's' : pf = PARTITION_SAUVAL_AND_TATUM(t, z)
  'c' : pf = PARTITION_COWLEY(t, z)    
ENDCASE  
  
RETURN, pf
END
