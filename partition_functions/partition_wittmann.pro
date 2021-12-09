;+
; NAME:
;
;   PARTITION_WITTMANN
;
; PURPOSE:
;
;   This procedure interpolates computes the partition function of a given
;   element and a given temperature using the pretabulated coefficients 
;   of Wittmann (1974, SoPh, 35, 11). 
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
;   pf = PARTITION_WITTMANN(t, z)
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
;
; COMMENT:
;
; DEPENDENCIES:
;
;
; EXAMPLE:
;
;   pf = PARTITION_WITTMANN([6E3, 7E3, 1E4], 7)
;   PRINT, pf.u1
;   PRINT, pf.u2
;   PRINT, pf.u3
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (July2013)
;               # based on a routine from the MISMA code by Jorge Sanchez
;                 Almeida (IAC).
;
;-
;================================================================================
; PARTITION_WITTMANN, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in this routine are taken from Wittmann (1974, SoPh, 35, 11). If you
; use the code and the data please acknowledge it by citing the original paper.
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

FUNCTION partition_wittmann, t, z, hirwin = hirwin

x = ALOG(5040./t)
y = 1.D-3*t

zerroarray = t*0.

CASE z OF

  ; H
  1 : BEGIN
        u1 = zerroarray+2 
        pos = WHERE((t GT 1.3D4) AND (t LT 1.62D4))
        IF (pos[0] NE -1) THEN u1[pos] = 1.51+3.8D-5*t[pos]  
        pos = WHERE(t GT 1.62D4)
        IF (pos[0] NE -1) THEN u1[pos] = 11.41d0+t[pos]*(-1.1428d-3+t[pos]*3.52d-8)
        u2 = zerroarray+1
        u3 = zerroarray-1
      END

  ; He
  2 : BEGIN
        u1 = zerroarray + 1.
        u2 = zerroarray + 2.
        u3 = zerroarray + 1.
      END

  ; Li
  3 : BEGIN
        u1 = 2.081D0 - y*(6.8926D-2 - y*1.4081D-2) 
        pos = WHERE(t GT 6.D3) 
        IF (pos[0] NE -1) THEN u1[pos] = 3.4864D0 + t[pos]*(-7.3292D-4 + t[pos]*8.5586D-8)
        u2 = zerroarray + 1.
        u3 = zerroarray + 2.
      END

  ; Be
  4 : BEGIN
        u1 = (zerroarray + 1.) > (0.631D0+7.032D-5*t)  ; IDL choose the largest at each 
        u2 = zerroarray + 2.
        u3 = zerroarray + 1.
  END

  ; B
  5 : BEGIN
        u1 = 5.9351D0 + 1.0438D-2*y 
        u2 = zerroarray + 1.
        u3 = zerroarray + 2.
      END

  ; C
  6 : BEGIN
        u1 = 8.6985D0 + y*(2.0485D-2 + y*(1.7629D-2 - 3.9091D-4*y)) 
        pos = WHERE(t GT 1.2D4)
        IF (pos[0] NE -1) THEN u1[pos] = 13.97D0 + t[pos]*(-1.3907D-3 + t[pos]*9.0844D-8)
        u2 = 5.838D0 + 1.6833D-5*t
        pos = WHERE(t GT 2.4D4)
        IF (pos[0] NE -1) THEN u2[pos] = 10.989D0 + t[pos]*(-6.9347D-4 + t[pos]*2.0861D-8)
        u3 = zerroarray + 1.
        pos = WHERE(t GT 1.95D4)
        IF (pos[0] NE -1) THEN u3[pos] = -0.555D0 + 8D-5*t[pos] ELSE u3 = zerroarray + 1.
      END

  ; N
  7 : BEGIN
        u1 = 3.9914D0 + y*(1.7491D-2 - y*(1.0148D-2 - y*1.7138D-3)) 
        pos = WHERE(t GT 8800.D0 AND t LE 1.8D4)
           IF (pos[0] NE -1) THEN u1[pos] = 2.171D0 + 2.54D-4*t[pos]
        pos = WHERE(t GT 1.8D4)
        IF (pos[0] NE -1) THEN u1[pos] = 11.396D0 + t[pos]*(-1.7139D-3 + t[pos]*8.633D-8)
        u2 = 8.060D0 + 1.420D-4*t
        pos = WHERE(t GT 3.3D4) 
        IF (pos[0] NE -1) THEN u2[pos] = 26.793D0 + t[pos]*(-1.8931D-3 + t[pos]*4.4612D-8)
        u3 = 5.9835D0 + t*(-2.6651D-5 + t*1.8228D-9)
        pos = WHERE(t LT 7310.5D0)
        IF (pos[0] NE -1) THEN u3[pos] = 5.89D0
      END

  ; O
  8 : BEGIN
        u1 = 8.29D0 + 1.10D-4*t 
        pos = WHERE(t GT 1.9D4)
        IF (pos[0] NE -1) THEN u1[pos] = 66.81D0 + t[pos]*(-6.019D-3 + t[pos]*1.657D-7)
        u2 = (zerroarray + 4.) > (3.51D0 + 8.D-5*t)
        pos = WHERE(t GT 3.64D4)
        IF (pos[0] NE -1) THEN u2[pos] = 68.7D0 + t[pos]*(-4.216D-3 + t[pos]*6.885D-8)
        u3 = 7.865D0 + 1.1348D-4*t
      END

  ; F
  9 : BEGIN
        u1 = 4.5832D0 + y*(0.77683D0 + y*(-0.20884D0 + y*(2.6771D-2 - 1.3035D-3*y))) 
        pos = WHERE(t GT 8750.D0 AND t LE 2.D4)
        IF (pos[0] NE -1) THEN u1[pos] = 5.9D0 + 0.*pos
        pos = WHERE(t GT 2D4)
        IF (pos[0] NE -1) THEN u1[pos] = 15.16D0 + t[pos]*(-9.229D-4 + t[pos]*2.312D-8)
        u2 = 8.15D0 + 8.9D-5*t
        u3 = 2.315D0 + 1.38D-4*t
      END

  ; Ne
  10 : BEGIN
         u1 = zerroarray + 1.
         pos = WHERE(t GT 2.69D4)
         IF (pos[0] NE -1) THEN u1[pos] = 26.3D0 + t[pos]*(-2.113D-3 + t[pos]*4.359D-8)
         u2 = 5.4D0 + 4.D-5*t
         u3 = 7.973D0 + 7.956D-5*t
       END

  ; Na
  11 : BEGIN
         u1 =  (zerroarray + 2.) > (1.72D0 + 9.3D-5*t)
         pos = WHERE(t GT 5400.D0 AND t LE 8.5D3)
         IF (pos[0] NE -1) THEN u1[pos] = -0.83D0 + 5.66D-4*t[pos]
         pos = WHERE(t GT 8.5D3)
         IF (pos[0] NE -1) THEN u1[pos] = 4.5568D0 + t[pos]*(-1.2415D-3 + t[pos]*1.3861D-7)
         u2 = zerroarray + 1.
         u3 = 5.69D0 + 5.69D-6*t
       END

  ; Mg
  12 : BEGIN
         u1 = 1. + EXP(-4.027262D0 - x*(6.173172D0 + x*(2.889176D0 + x*(2.393895D0 + 0.784131D0*x))))
         pos = WHERE(t GT 8D3)
         IF (pos[0] NE -1) THEN u1[pos] = 2.757D0 + t[pos]*(-7.8909D-4 + t[pos]*7.4531D-8)
         u2 = 2.D0 + EXP(-7.721172D0 - x*(7.600678D0 + x*(1.966097D0 + 0.212417D0*x)))
         pos = WHERE(t GT 2.D4)
         IF (pos[0] NE -1) THEN u2[pos] = 7.1041D0 + t[pos]*(-1.0817D-3 + t[pos]*4.7841D-8)
         u3 = zerroarray + 1.
       END

  ; Al
  13 : BEGIN
         u1 = 5.2955D0 + y*(0.27833D0 - y*(4.7529D-2 - y*3.0199D-3)) 
         u2 = (zerroarray + 1.) > (0.725D0 + 3.245D-5*t)
         pos = WHERE(t GT 2.24D4)
         IF (pos[0] NE -1) THEN u2[pos] = 61.06D0 + t[pos]*(-5.987D-3 + t[pos]*1.485D-7)
         u3 = (zerroarray + 2.) > (1.976D0 + 3.43D-6*t)
         pos = WHERE(t GT 1.814D4)
         IF (pos[0] NE -1) THEN u3[pos] = 3.522D0 + t[pos]*(-1.59D-4 + t[pos]*4.382D-9)
       END

  ; Si
  14 : BEGIN
         u1 = 6.7868D0 + y*(0.86319D0 + y*(-0.11622D0 + y*(0.013109D0 - 6.2013D-4*y))) 
         pos = WHERE(t GT 1.04D4)
         IF (pos[0] NE -1) THEN u1[pos] = 86.01D0 + t[pos]*(-1.465D-2 + t[pos]*7.282D-7)
         u2 = 5.470D0 + 4.D-5*t
         pos = WHERE(t GT 1.8D4)
         IF (pos[0] NE -1) THEN u2[pos] = 26.44D0 + t[pos]*(-2.22D-3 + t[pos]*6.188D-8)
         u3 = (zerroarray + 1.) > (0.911D0 + 1.1D-5*t)
         pos = WHERE(t GT 3.33D4)
         IF (pos[0] NE -1) THEN u3[pos] = 19.14D0 + t[pos]*(-1.408D-3 + t[pos]*2.617D-8)
       END

  ; P
  15 : BEGIN
         u1 = 4.2251D0 + y*(-0.22476D0 + y*(0.057306D0 - y*1.0381D-3)) 
         pos = WHERE(t GT 6.D3)
         IF (pos[0] NE -1) THEN u1[pos] = 1.56D0 + 5.2D-4*t[pos]
         u2 = 4.4151D0 + y*(2.2494D0 + y*(-0.55371D0 + y*(0.071913D0 - y*3.5156D-3)))
         pos = WHERE(t GT 7250.D0)
         IF (pos[0] NE -1) THEN u2[pos] = 4.62D0 + 5.38D-4*t[pos]
         u3 = 5.595D0 + 3.4D-5*t
       END

  ; S
  16 : BEGIN
         u1 = 7.5D0 + 2.15D-4*t 
         pos = WHERE(t GT 1.16D4)
         IF (pos[0] NE -1) THEN u1[pos] = 38.76D0 + t[pos]*(-4.906D-3 + t[pos]*2.125D-7)
         u2 = 2.845D0 + 2.43D-4*t
         pos = WHERE(t GT 1.05D4)
         IF (pos[0] NE -1) THEN u2[pos] = 6.406D0 + t[pos]*(-1.68D-4 + t[pos]*1.323D-8)
         u3 = 7.38D0 + 1.88D-4*t
       END

  ; Cl
  17 : BEGIN
         u1 = 5.2D0 + 6.D-5*t 
         pos = WHERE(t GT 1.84D4)
         IF (pos[0] NE -1) THEN u1[pos] = -81.6D0 + 4.8D-3*t[pos]
         u2 = 7.0D0 + 2.43D-4*t
         u3 = 2.2D0 + 2.62D-4*t
       END

  ; Ar
  18 : BEGIN
         u1 = zerroarray + 1.
         u2 = 5.20D0 + 3.8D-5*t
         u3 = 7.474D0 + 1.554D-4*t
       END
 
  ; K
  19 : BEGIN
         u1 = 1.9909D0 + y*(0.023169D0 - y*(0.017432D0 - y*4.0938D-3)) 
         pos = WHERE(t GT 5800.D0)
         IF (pos[0] NE -1) THEN u1[pos] = -9.93D0 + 2.124D-3*t[pos]
         u2 = zerroarray + 1.
         u3 = 5.304D0 + 1.93D-5*t
       END

  ; Ca
  20 : BEGIN
         u1 = 1.D0 + EXP(-1.731273D0 - x*(5.004556D0 + x*(1.645456D0 + x*(1.326861D0+0.508553D0*x))))
         u2 = 2.D0 + EXP(-1.582112D0 - x*(3.996089D0 + x*(1.890737D0 + 0.539672D0*x)))
         u3 = zerroarray + 1.
       END

  ; Sc
  21 : BEGIN
         u1 = 4.D0 + EXP(2.071563D0 + x*(-1.2392D0 + x*(1.173504D0 + .0517796D0*x))) 
         u2 = 3.D0 + EXP(2.988362D0 + x*(-0.596238D0 + .054658D0*x))
         u3 = zerroarray + 10.
       END

  ; Ti
  22 : BEGIN
         u1 = 5. + EXP(3.200453D0 + x*(-1.227798D0 + x*(0.799613D0 + 0.278963*x))) 
         pos = WHERE(t LT 5.5D3)
         IF (pos[0] NE -1) THEN u1[pos] = 16.37D0 + t[pos]*(-2.838D-4 + t[pos]*5.819D-7)
         u2 = 4.D0 + EXP(3.94529D0 + x*(-0.551431D0 + 0.115693D0*x))
         u3 = 16.4D0 + 8.5D-4*t
       END

  ; V
  23 : BEGIN
         u1 = 4.D0 + EXP(3.769611D0 + x*(-0.906352D0 + x*(0.724694D0 + 0.1622D0*x))) 
         u2 = 1.D0 + EXP(3.755917D0 + x*(-0.757371D0 + 0.21043D0*x))
         u3 = -18.D0 + 1.03D-2*t
         pos = WHERE(t LT 2.25D3)
         IF (pos[0] NE -1) THEN u3[pos] = 2.4D-3*t[pos]
       END

  ; Cr
  24 : BEGIN
         u1 = 7.D0 + EXP(1.225042D0 + x*(-2.923459D0 + x*(0.154709D0 + 0.09527D0*x))) 
         u2 = 6.D0 + EXP(0.128752D0 - x*(4.143973D0 + x*(1.096548D0 + 0.230073D0*x)))
         u3 = 10.4D0 + 2.1D-3*t
       END

  ; Mn
  25 : BEGIN
         u1 = 6.D0 + EXP(-0.86963D0 - x*(5.531252D0 + x*(2.13632D0 + x*(1.061055D0 + 0.265557D0*x)))) 
         u2 = 7.D0 + EXP(-0.282961D0 - x*(3.77279D0 + x*(0.814675D0 + 0.159822D0*x)))
         u3 = zerroarray + 10.
       END

  ; Fe
  26 : BEGIN
         u1 = 9.D0 + EXP(2.930047D0 + x*(-0.979745D0 + x*(0.76027D0 + 0.118218D0*x))) 
         pos = WHERE(t LT 4D3)
         IF (pos[0] NE -1) THEN u1[pos] = 15.85D0 + t[pos]*(1.306D-3 + t[pos]*2.04D-7)
         pos = WHERE(t GT 9D3)
         IF (pos[0] NE -1) THEN u1[pos] = 39.149D0 + t[pos]*(-9.5922D-3 + t[pos]*1.2477D-6)
         u2 = 10.D0 + EXP(3.501597 + x*(-0.612094 + 0.280982*x))
         pos = WHERE(t GT 1.8D4)
         IF (pos[0] NE -1) THEN u2[pos] = 68.356D0 + t[pos]*(-6.1104D-3 + t[pos]*5.1567D-7)
         u3 = 17.336D0 + t*(5.5048D-4 + t*5.7514D-8)
       END

  ; Co
  27 : BEGIN
         u1 = 8.65D0 + 4.9D-3*t 
         u2 = 11.2D0 + 3.58D-3*t
         u3 = 15.0D0 + 1.42D-3*t
       END

  ; Ni
  28 : BEGIN
         u1 = 9.D0 + EXP(3.084552D0 + x*(-0.401323D0 + x*(0.077498D0 - 0.278468D0*x))) 
         u2 = 6.D0 + EXP(1.593047D0 - x*(1.528966D0 + 0.115654D0*x))
         u3 = 13.3D0 + 6.9D-4*t
       END

  ; Cu
  29 : BEGIN
         u1 =  (zerroarray + 2.) > (1.50D0 + 1.51D-4*t) 
         pos = WHERE(t GT 6250.D0)
         IF (pos[0] NE -1) THEN u1[pos] = -0.3D0 + 4.58D-4*t[pos]
         u2 =  (zerroarray + 1.) > (0.22D0 + 1.49D-4*t)
         u3 = 8.025D0 + 9.4D-5*t
       END

  ; Zn
  30 : BEGIN
         u1 =  (zerroarray + 1.) > (0.632D0 + 5.11D-5*t) 
         u2 = zerroarray + 2.
         u3 = zerroarray + 1.
       END

  ; Ga
  31 : BEGIN
         u1 = 1.7931D0 + y*(1.9338D0 + y*(-0.4643D0 + y*(0.054876D0 - y*2.5054D-3))) 
         pos = WHERE(t GT 6.D3)
         IF (pos[0] NE -1) THEN u1[pos] = 4.18D0 + 2.03D-4*t[pos]
         u2 = zerroarray + 1.
         u3 = zerroarray + 2.
       END

  ; Ge
  32 : BEGIN
         u1 = 6.12D0 + 4.08D-4*t 
         u2 = 3.445D0 + 1.78D-4*t
         u3 = zerroarray + 1.1
       END

  ; As
  33 : BEGIN
         u1 = 2.65D0 + 3.65D-4*t 
         u2 = -0.25384D0 + y*(2.284D0 + y*(-0.33383D0 + y*(0.030408D0 - y*1.1609D-3)))
         pos = WHERE(t GT 1.2D4)
         IF (pos[0] NE -1) THEN u2[pos] = 8. + 0.*pos
         u3 = zerroarray + 8.
       END

  ; Se
  34 : BEGIN
         u1 = 6.34D0 + 1.71D-4*t 
         u2 = 4.1786D0 + y*(-0.15392D0 + 3.2053D-2*y)
         u3 = zerroarray + 8.
       END

  ; Br
  35 : BEGIN
         u1 = 4.12D0 + 1.12D-4*t 
         u2 = 5.22D0 + 3.08D-4*t
         u3 = 2.3D0 + 2.86D-4*t
       END

  ; Kr
  36 : BEGIN
         u1 = zerroarray + 1.
         u2 = 4.11D0 + 7.4D-5*t
         u3 = 5.35D0 + 2.23D-4*t
       END

  ; Rb
  37 : BEGIN
         u1 = (zerroarray + 2.) > (1.38D0 + 1.94D-4*t) 
         pos = WHERE(t GT 6250.D0)
         IF (pos[0] NE -1) THEN u1[pos] = -14.9D0 + 2.79D-3*t[pos]
         u2 = zerroarray + 1.
         u3 = 4.207D0 + 4.85D-5*t
       END

  ; Sr
  38 : BEGIN
         u1 = 0.87127D0 + y*(0.20148D0 + y*(-0.10746D0 + y*(0.021424D0 - y*1.0231D-3))) 
         pos = WHERE(t GT 6500.D0)
         IF (pos[0] NE -1) THEN u1[pos] = -6.12D0 + 1.224D-3*t[pos]
         u2 = (zerroarray + 2.) > (0.84D0 + 2.6D-4*t)
         u3 = zerroarray + 1. 
       END

  ; Y
  39 : BEGIN
         u1 = 0.2D0 + 2.58D-3*t 
         u2 = 7.15D0 + 1.855D-3*t
         u3 = 9.71D0 + 9.9D-5*t
       END

  ; Zc
  40 : BEGIN
         u1 = 76.31D0 + t*(-1.866D-2 + t*2.199D-6)  
         pos = WHERE(t LT 6236.D0)
         IF (pos[0] NE -1) THEN u1[pos] = 6.8D0 + t[pos]*(2.806D-3 + t[pos]*5.386D-7)
         u2 = 4.D0 + EXP(3.721329D0 - 0.906502*x)
         u3 = 12.3D0 + 1.385D-3*t
       END

  ; Nb
  41 : BEGIN
         u1 = (zerroarray + 1.) > (-19.D0 + 1.43D-2*t)
         u2 = -4.D0 + 1.015D-2*t
         u3 = zerroarray + 25.
       END

  ;Mo
  42 : BEGIN
         u1 = (zerroarray + 7.) > (2.1D0 + 1.5D-3*t) 
         pos = WHERE(t GT 7.D3)
         IF (pos[0] NE -1) THEN u1[pos] = -38.1D0 + 7.28D-3*t[pos]
         u2 = 1.25D0 + 1.17D-3*t
         pos = WHERE(t GT 6900.)
         IF (pos[0] NE -1) THEN u2[pos] = -28.5D0 + 5.48D-3*t[pos]
         u3 = 24.04D0 + 1.464D-4*t
       END

  ; Tc
  43 : BEGIN
         u1 = 4.439D0 + y*(0.30648D0 + y*(1.6525D0 + y*(-0.4078D0 + y*(0.048401 - y*2.1538D-3)))) 
         pos = WHERE(t GT 6.D3)
         IF (pos[0] NE -1) THEN u1[pos] = 24.D0 + 0.*pos
         u2 = 8.1096D0 + y*(-2.963D0 + y*(2.369D0 + y*(-0.502D0 + y*(0.049656D0 - y*1.9087D-3))))
         pos = WHERE(t GT 6.D3)
         IF (pos[0] NE -1) THEN u2[pos] = 17.D0 + 0.*pos
         u3 = zerroarray + 220.
       END

  ; Ru
  44 : BEGIN
         u1 = -3.D0 + 7.17D-3*t 
         u2 = 3.D0 + 4.26D-3*t
         u3 = zerroarray + 22.
       END

  ; Rh
  45 : BEGIN
         u1 = 6.9164D0 + y*(3.8468D0 + y*(0.043125D0 - y*(8.7907D-3 - y*5.9589D-4))) 
         u2 = 7.2902D0 + y*(1.7476D0 + y*(-0.038257D0 + y*(2.014D-3 + y*2.1218D-4)))
         u3 = zerroarray + 30.
       END

  ; Pd
  46 : BEGIN
         u1 = (zerroarray + 1.) > (-1.75D0 + 9.86D-4*t) 
         u2 = 5.60D0 + 3.62D-4*t
         u3 = zerroarray + 20.
       END

  ; Ag
  47 : BEGIN
         u1 = (zerroarray + 2.) > (1.537D0 + 7.88D-5*t) 
         u2 = (zerroarray + 1.) > (0.73D0 + 3.4D-5*t)
         u3 = 6.773D0 + 1.248D-4*t
       END

  ; Cd
  48 : BEGIN
         u1 = (zerroarray + 1.) > (0.43D0 + 7.6D-5*t) 
         u2 = zerroarray + 2.
         u3 = zerroarray + 1.
       END

  ; In
  49 : BEGIN
         u1 = 2.16D0 + 3.92D-4*t 
         u2 = zerroarray + 1.
         u3 = zerroarray + 2.
       END

  ; Sn
  50 : BEGIN
         u1 = 2.14D0 + 6.16D-4*t 
         u2 = 2.06D0 + 2.27D-4*t
         u3 = zerroarray + 1.05
       END

  ; Sb
  51 : BEGIN
         u1 = 2.34D0 + 4.86D-4*t 
         u2 = 0.69D0 + 5.36D-4*t
         u3 = zerroarray + 3.5
       END

  ; Te
  52 : BEGIN
         u1 = 3.948D0 + 4.56D-4*t 
         u2 = 4.2555D0 + y*(-0.25894D0 + y*(0.06939D0 - y*2.4271D-3))
         pos = WHERE(t GT 1.2D4)
         IF (pos[0] NE -1) THEN u2[pos] = 7.D0 + 0.*pos
         u3 = zerroarray + 5.
       END

  ; I
  53 : BEGIN
         u1 = (zerroarray + 4.0) > (3.8D0 + 9.5D-5*t) 
         u2 = 4.12D0 + 3.D-4*t
         u3 = zerroarray + 7.0
       END

  ; Xe
  54 : BEGIN
         u1 = zerroarray + 1. 
         u2 = 3.75D0 + 6.876D-5*t
         u3 = 4.121D0 + 2.323D-4*t
       END

  ; Cs
  55 : BEGIN
         u1 = (zerroarray + 2.) > (1.56D0 + 1.67D-4*t) 
         pos = WHERE(t GT 4850.D0)
         IF (pos[0] NE -1) THEN u1[pos] = -2.680D0 + 1.04D-3*t[pos]
         u2 = zerroarray + 1.
         u3 = 3.769D0 + 4.971D-5*t
       END

  ; Ba
  56 : BEGIN
         u1 = (zerroarray + 1.) > (-1.8D0 + 9.85D-4*t) 
         pos = WHERE(t GT 6850.D0)
         IF (pos[0] NE -1) THEN u1[pos] = -16.2D0 + 3.08D-3*t[pos]
         u2 = 1.11D0 + 5.94D-4*t
         u3 = zerroarray + 1.
       END

  ; La
  57 : BEGIN
         u1 = 15.42D0 + 9.5D-4*t 
         pos = WHERE(t GT 5060.D0)
         IF (pos[0] NE -1) THEN u1[pos] = 1.D0 + 3.8D-3*t[pos]
         u2 = 13.2D0 + 3.56D-3*t
         u3 = zerroarray + 12. 
       END

  ; Ce
  58 : BEGIN
         u1 = 9.D0 + EXP(5.202903D0 + x*(-1.98399D0 + x*(0.119673D0 + 0.179675D0*x))) 
         u2 = 8.D0 + EXP(5.634882D0 - x*(1.459196D0 + x*(0.310515D0 + 0.052221D0*x)))
         u3 = 9.D0 + EXP(3.629123D0 - x*(1.340945D0 + x*(0.372409D0 + x*(0.03186D0 - 0.014676D0*x))))
       END

  ; Pr 
  59 : BEGIN
         u2 = 9.D0 + EXP(4.32396D0 - x*(1.191467D0 + x*(0.149498D0 + 0.028999*x))) 
         u1 = u2 
         u3 = 10.D0 + EXP(3.206855D0 + x*(-1.614554D0 + x*(0.489574D0 + 0.277916*x)))
       END

  ; Nb
  60 : BEGIN
         u1 = 9.D0 + EXP(4.456882D0 + x*(-2.779176D0 + x*(0.082258D0 + x*(0.50666D0 + 0.127326D0*x)))) 
         u2 = 8.D0 + EXP(4.689643D0 + x*(-2.039946D0 + x*(0.17193D0 + x*(0.26392D0 + 0.038225D0*x))))
         u3 = u2 
       END

  ; Pm
  61 : BEGIN
         u1 = zerroarray + 20. 
         u2 = zerroarray + 25.
         u3 = zerroarray + 100.
       END

  ; Sm
  62 : BEGIN
         u1 = 1.D0 + EXP(3.549595D0 + x*(-1.851549D0 + x*(0.9964D0 + 0.566263D0*x))) 
         u2 = 2.D0 + EXP(4.052404D0 + x*(-1.418222D0 + x*(0.358695D0 + 0.161944D0*x)))
         u3 = 1.D0 + EXP(3.222807D0 - x*(0.699473D0 + x*(-0.056205 + x*(0.533833D0 + 0.251011D0*x))))
       END

  ; Eu
  63 : BEGIN
         u1 = 8.D0 + EXP(1.024374D0 - x*(4.533653D0 + x*(1.540805D0 + x*(0.827789D0 + 0.286737D0*x))))
         u2 = 9.D0 + EXP(1.92776D0 + x*(-1.50646D0 + x*(0.379584D0 + 0.05684D0*x)))
         u3 = zerroarray + 8.
       END

  ; Gd
  64 : BEGIN
         u1 = 5.D0 + EXP(4.009587D0 + x*(-1.583513D0 + x*(0.800411D0 + 0.388845D0*x))) 
         u2 = 6.D0 + EXP(4.362107D0 - x*(1.208124D0 + x*(-0.074813D0 + x*(0.076453D0 + 0.055475D0*x))))
         u3 = 5.D0 + EXP(3.412951D0 - x*(0.50271D0 + x*(0.042489D0 - 4.017D-3*x)))
       END

  ; Tb
  65 : BEGIN
         u1 = 16.D0 + EXP(4.791661D0 + x*(-1.249355D0 + x*(0.570094D0 + 0.240203D0*x))) 
         u2 = 15.D0 + EXP(4.472549D0 - x*(0.295965D0 + x*(5.88D-3 + 0.131631D0*x)))
         u3 = u2 
       END

  ; Dy
  66 : BEGIN
         u1 = 17.D0 + EXP(3.029646D0 - x*(3.121036D0 + x*(0.086671D0 - 0.216214D0*x))) 
         u2 = 18.D0 + EXP(3.465323D0 - x*(1.27062D0 + x*(-0.382265D0 + x*(0.431447D0 + 0.303575D0*x))))
         u3 = u2 
       END

  ; Ho
  67 : BEGIN
         u3 = 16.D0 + EXP(1.610084D0 - x*(2.373926D0 + x*(0.133139D0 - 0.071196D0*x))) 
         u1 = u3
         u2 = u3 
       END

  ; Er
  68 : BEGIN
         u1 = 13.D0 + EXP(2.895648D0 - x*(2.968603D0 + x*(0.561515D0 + x*(0.215267D0 + 0.095813D0*x)))) 
         u2 = 14.D0 + EXP(3.202542D0 - x*(0.852209D0 + x*(-0.226622D0 + x*(0.343738D0 + 0.186042D0*x))))
         u3 = u2 
       END

  ; Tm
  69 : BEGIN
         u1 = 8.D0 + EXP(1.021172D0 - x*(4.94757D0 + x*(1.081603D0 + 0.034811D0*x))) 
         u2 = 9.D0 + EXP(2.173152D0 + x*(-1.295327D0 + x*(1.940395D0 + 0.813303D0*x)))
         u3 = 8.D0 + EXP(-0.567398D0 + x*(-3.383369D0 + x*(0.799911D0 + 0.554397D0*x)))
       END

  ; Yb
  70 : BEGIN
         u1 = 1.D0 + EXP(-2.350549D0 - x*(6.688837D0 + x*(1.93869D0 + 0.269237D0*x))) 
         u2 = 2.D0 + EXP(-3.047465D0 - x*(7.390444D0 + x*(2.355267D0 + 0.44757D0*x)))
         u3 = 1.D0 + EXP(-6.192056D0 - x*(10.560552D0 + x*(4.579385D0 + 0.940171D0*x)))
       END

  ; Lu
  71 : BEGIN
         u1 = 4.D0 + EXP(1.537094D0 + x*(-1.140264D0 + x*(0.608536D0 + 0.193362D0*x))) 
         u2 = MAX([1.D0, 0.66D0 + 1.52D-4*t])
         pos = WHERE(t GT 5250.D0)
         IF (pos[0] NE -1) THEN u2[pos] = -1.09D0 + 4.86D-4*t[pos]
         u3 = zerroarray + 5. 
       END

  ; Hf
  72 : BEGIN
         u1 = 4.1758D0 + y*(0.407D0 + y*(0.57862D0 - y*(0.072887D0 - y*3.6848D-3))) 
         u2 = -2.979D0 + 3.095D-3*t
         u3 = zerroarray + 30.
       END

  ; Ta
  73 : BEGIN
         u1 = 3.0679D0 + y*(0.81776D0 + y*(0.34936D0 + y*(7.4861D-3 + y*3.0739D-4))) 
         u2 = 1.6834D0 + y*(2.0103D0 + y*(0.56443D0 - y*(0.031036D0 - y*8.9565D-4)))
         u3 = zerroarray + 15.
       END

  ; W
  74 : BEGIN
         u1 = 0.3951D0 + y*(-0.25057D0 + y*(1.4433D0 + y*(-0.34373D0 + y*(0.041924D0 - y*1.84D-3)))) 
         pos = WHERE(t GT 1.2D4)
         IF (pos[0] NE -1) THEN u1[pos] = 23.D0 + 0.*pos
         u2 = 1.055D0 + y*(1.0396D0 + y*(0.3303D0 - y*(8.4971D-3 - y*5.5794D-4)))
         u3 = zerroarray + 20.
       END

  ; Re
  75 : BEGIN
         u1 = 5.5671D0 + y*(0.72721D0 + y*(-0.42096D0 + y*(0.09075D0 - y*3.9331D-3))) 
         pos = WHERE(t GT 1.2D4)
         IF (pos[0] NE -1) THEN u1[pos] = 29.D0 + 0.*pos
         u2 = 6.5699D0 + y*(0.59999D0 + y*(-0.28532D0 + y*(0.050724D0 - y*1.8544D-3)))
         pos = WHERE(t GT 1.2D4)
         IF (pos[0] NE -1) THEN u2[pos] = 22.D0 + 0*pos
         u3 = zerroarray + 20.
       END

  ; Os
  76 : BEGIN
         u1 = 8.6643D0 + y*(-0.32516D0 + y*(0.68181D0 - y*(0.044252D0 - y*1.9975D-3))) 
         u2 = 9.7086D0 + y*(-0.3814D0 + y*(0.65292D0 - y*(0.064984D0 - y*2.8792D-3)))
         u3 = zerroarray + 10.
       END

  ; Ir
  77 : BEGIN
         u1 = 11.07D0 + y*(-2.412D0 + y*(1.9388D0 + y*(-0.34389D0 + y*(0.033511D0-1.3376D-3*y)))) 
         pos = WHERE(t GT 1.2D4)
         IF (pos[0] NE -1) THEN u1[pos] = 30.D0 + 0.*pos
         u2 = zerroarray + 15.
         u3 = zerroarray + 20.
       END

  ; Pt
  78 : BEGIN
         u1 = 16.4D0 + 1.27D-3*t 
         u2 = 6.5712D0 + y*(-1.0363D0 + y*(0.57234D0 - y*(0.061219D0 - 2.6878D-3*y)))
         u3 = zerroarray + 15.
       END

  ; Au
  79 : BEGIN
         u1 = 1.24D0 + 2.79D-4*t 
         u2 = 1.0546D0 + y*(-0.040809D0 + y*(2.8439D-3 + y*1.6586D-3))
         u3 = zerroarray + 7.
       END

  ; Hg
  80 : BEGIN
         u1 = zerroarray + 1.
         u2 = zerroarray + 2.
         u3 = (zerroarray + 1.) > (0.669D0 + 3.976D-5*t)
       END

  ; Tl
  81 : BEGIN
         u1 = (zerroarray + 2.) > (0.63D0 + 3.35D-4*t) 
         u2 = zerroarray + 1.
         u3 = zerroarray + 2.
       END

  ; Pb
  82 : BEGIN
         u1 = (zerroarray + 1.) > (0.42D0 + 2.35D-4*t) 
         pos = WHERE(t GT 6125.D0)
         IF (pos[0] NE -1) THEN u1[pos] = -1.2D0 + 5.D-4*t[pos]
         u2 = (zerroarray + 2.) > (1.72D0 + 7.9D-5*t)
         u3 = zerroarray + 1.
       END

  ; Bi
  83 : BEGIN
         u1 = 2.78D0 + 2.87D-4*t 
         u2 = (zerroarray + 1.) > (0.37D0 + 1.41D-4*t)
         u3 = zerroarray + 2.5
       END

  ; Po
  84 : BEGIN
         u1 = zerroarray + 5.
         u2 = zerroarray + 5.
         u3 = zerroarray + 4.
       END

  ; At
  85 : BEGIN
         u1 = zerroarray + 4.
         u2 = zerroarray + 6.
         u3 = zerroarray + 6.
       END

  ; Rn
  86 : BEGIN
         u1 = zerroarray + 1.
         u2 = zerroarray + 4.
         u3 = zerroarray + 6.
       END

  ; Fr
  87 : BEGIN
         u1 = zerroarray + 2.
         u2 = zerroarray + 1.
         u3 = zerroarray + 4.5
       END

  ; Ra
  88 : BEGIN
         u1 = zerroarray + 1.
         u2 = zerroarray + 2.
         u3 = zerroarray + 1.
       END

  ; Ac
  89 : BEGIN
         u1 = zerroarray + 6.
         u2 = zerroarray + 3.
         u3 = zerroarray + 7.
       END

  ; Th
  90 : BEGIN
         u1 = zerroarray + 8.
         u2 = zerroarray + 8.
         u3 = zerroarray + 8.
       END

  ; Pa
  91 : BEGIN
         u1 = zerroarray + 50.
         u2 = zerroarray + 50.
         u3 = zerroarray + 50.
       END

  ; U
  92 : BEGIN
         u1 = zerroarray + 25.
         u2 = zerroarray + 25.
         u3 = zerroarray + 25.
       END

  ELSE : PRINT, 'PARTITION: no available data, Z > 92'
ENDCASE

result = {u1:u1, u2:u2, u3:u3}
RETURN, result

END
