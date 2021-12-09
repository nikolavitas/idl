;+
; NAME:
;       LOAD_ABUNDANCES
;
; PURPOSE:
;
;       This function loads the photospheric abundances from various sources
;       for atomic numbers between Z = 1 and Z = 92. By default, the returned 
;       relative abundances are in the logarithmic scale where the abundance 
;       of hydrogen is set to 12. If the keyword PPM is set, the returned 
;       abundances are given in parts per milion relative to the hydrogen
;       abundance. If the keyword XYZ is specified, instead of the abundances
;       the function computes and returns the mass fraction X, Y and Z (X for
;       H, Y for He and Z for all other metals. The keyword SOURCE is 
;       mandatory. User has to choose one of the available measurements. If
;       the keyword is not specified, the function asks for it. The available
;       sources are:
;       
;       Name                        Code            ADS Reference
;       ---------------------------------------------------------------
;       Withbroe (1971)             W71             1971spas.conf..127W
;       Grevesse (1984)             G84             1984PhST....8...49G  
;       Anders & Grevesse (1989)    AG89            1989GeCoA..53..197A
;       Grevesse & Sauval (1998)    GS98            1998SSRv...85..161G
;       Asplund et al (2005)        A05             2005ASPC..336...25A
;       Asplund et al (2009)        A09             2009ARA&A..47..481A
;       Lodders et al (2009)        L09             2009LanB...4B...44L
;       
;       If you use any of these data please acknowledge it appropriately.
;
;       An additiotal set of abundances by Thevenin (1989A&AS...77..137T,  
;       1990A&AS...82..179T) is hidden in the code. I copied the data from
;       another code. However, after checking the original papers it came
;       clear that the set of values present in this function corresponds
;       mainly to the paper of Holweger and Mueller. Thevenin measured 
;       the abundance of ~20 elements, so even the correct set does not
;       cover the complete periodic system.
; 
;       Missing values for some elements are replaced by corresponding
;       values of the meteoritic abundances. If you need to know where it
;       is done, go to the original papers and compare the numbers.
;
;       Abundances from Withbroe (1971) are included solely for the 
;       historical reasons - this set was adopted for building the famous
;       VAL3C model of the solar atmosphere. However, note that Whitbroe
;       has values for only 20 elements in the list.
;
;       Warning: Any set of the solar abundances should not be used without
;       knowing the assumptions used to get it. Please read the original 
;       papers for that.
;
;       The atomic weights are loaded using the function LOAD_ATOMIC_WEIGHTS.
;       See the header of that function for the information on the data that
;       is used.
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
;       abund = LOAD_ABUNDANCES([z = z, source = source, ppm = ppm, xyz = xyz])
;
; INPUTS:
;
;       z      = String or array, integer. The atomic number(s). If not 
;                specified, the abundances for the first 92 elements are
;                returned. It must not be specified if the PPM keyword is set.
; 
; OUTPUTS:
;
;       abund = Double array (default), logarithmic abundances relative to H.
;             = Double array (if ppm is set), abundances relative to H in 
;               parts per milion.
;             = Structure with three floating tags (x, y, z) containing the
;               mass fractions if the XYZ keyword is specified.                     
;
; INPUT KEYWORDS:
; 
;       source = String. Code of the abundance reference (see the table above).
;
;       ppm    = Integer. If true, the output is changed to linear relative
;                abundances in ppm.
; 
;       xyz    = Integer. If true, the output is changed to a structure 
;                contraining the mass fractions.
; 
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
;       ; Logarithimic abundances for C. N and O relative to H from 
;       ; Asplund et al, 2009
;       IDL> a = LOAD_ABUNDANCES([6, 7, 8], source = 'A09')
;       IDL> PRINT, a
;       8.4300000       7.8300000       8.6900000
; 
;       ; Linear abundances for C. N and O relative to H [in ppm].
;       IDL> a = LOAD_ABUNDANCES([6, 7, 8], source = 'A09', /ppm)
;       IDL> PRINT, a
;       269.15348       67.608298       489.77882
; 
;       ; Mass fractions X, Y, Z computed using the abundances of
;       ; Asplund et al, 2009. 
;       IDL> a = LOAD_ABUNDANCES(source = 'A09', /xyz)               
;       IDL> PRINT, a
;       {      0.73738727      0.24924153     0.013371199}
;       
;
; DEPENDENCIES:
;
;       LOAD_ATOMIC_WEIGHTS
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October 2007.
;       - Withbroe added, August 2013.
;       - Thevenin removed, Asplund 09 and Lodders added (NV), other small 
;         fixes, October 2013.
;-
;================================================================================
; LOAD_ABUNDANCES, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
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

FUNCTION load_abundances, z, source = source, ppm = ppm, xyz = xyz
                          
options = ['W71', 'AG89', 'G84', 'GS98', 'A05', 'A09', 'L09']                          
IF ~KEYWORD_SET(source) THEN BEGIN
  source = ''
  PRINT, ' Specify the set of abundances to be loaded (type the code).'
  PRINT, ' Name                        Code            ADS Reference'
  PRINT, '---------------------------------------------------------------'
  PRINT, ' Withbroe (1971)             W71             1971spas.conf..127W'
  PRINT, ' Grevesse (1984)             G84             1984PhST....8...49G'  
  PRINT, ' Anders & Grevesse (1989)    AG89            1989GeCoA..53..197A'
  PRINT, ' Grevesse & Sauval (1998)    GS98            1998SSRv...85..161G'
  PRINT, ' Asplundet et al (2005)      A05             2005ASPC..336...25A'
  PRINT, ' Asplundet et al (2009)      A09             2009ARA&A..47..481A'
  PRINT, ' Lodders et al (2009)        L09             2009LanB...4B...44L'
  READ, source, prompt = 'Source: '
ENDIF

source = STRUPCASE(source)

IF TOTAL(source EQ options) EQ 0 THEN MESSAGE, 'Source is not available. STOP'

CASE source OF

  'W71'   :  abundance = [ 12.00D0, 11.00D0,  0.00D0,  0.00D0,  0.00D0,  8.73D0,  8.06D0,  8.83D0, $
                            0.00D0,  7.55D0,  6.45D0,  7.65D0,  6.45D0,  7.65D0,  5.45D0,  7.21D0, $
                            0.00D0,  6.75D0,  5.75D0,  6.40D0,  0.00D0,  0.00D0,  0.00D0,  5.85D0, $
                            5.55D0,  7.40D0,  5.35D0,  6.53D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.00D0,  0.00D0,  0.00D0 ]
                            
  'G84'   :  abundance = [ 12.00D0, 11.00D0,  1.00D0,  1.15D0,  2.60D0,  8.55D0,  7.99D0,  8.77D0, $
                            4.56D0,  8.00D0,  6.18D0,  7.48D0,  6.40D0,  7.55D0,  5.45D0,  7.21D0, $
                            5.50D0,  6.58D0,  5.12D0,  6.36D0,  3.10D0,  5.02D0,  4.00D0,  5.67D0, $
                            5.45D0,  7.50D0,  4.92D0,  6.25D0,  4.21D0,  4.60D0,  2.88D0,  3.63D0, $ 
                            2.39D0,  3.35D0,  2.63D0,  3.21D0,  2.60D0,  2.90D0,  2.24D0,  2.56D0, $ 
                            2.10D0,  1.92D0,  0.00D0,  1.84D0,  1.12D0,  1.69D0,  0.94D0,  1.86D0, $ 
                            1.66D0,  2.00D0,  1.00D0,  2.25D0,  1.51D0,  2.19D0,  1.12D0,  2.13D0, $ 
                            1.22D0,  1.55D0,  0.71D0,  1.34D0,  0.00D0,  0.80D0,  0.51D0,  1.12D0, $ 
                            0.20D0,  1.10D0,  0.26D0,  0.93D0,  0.00D0,  1.08D0,  0.76D0,  0.88D0, $
                           -0.09D0,  1.11D0,  0.26D0,  1.45D0,  1.35D0,  1.80D0,  1.13D0,  1.27D0, $ 
                            0.90D0,  1.90D0,  0.71D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.02D0,  0.00D0, -0.47D0 ]

  'AG89'  :  abundance = [ 12.00D0, 10.99D0,  1.16D0,  1.15D0,  2.88D0,  8.56D0,  8.05D0,  8.93D0, $
                            4.56D0,  8.09D0,  6.33D0,  7.58D0,  6.47D0,  7.55D0,  5.45D0,  7.21D0, $
                            5.50D0,  6.56D0,  5.12D0,  6.36D0,  3.10D0,  4.99D0,  4.00D0,  5.67D0, $
                            5.39D0,  7.67D0,  4.92D0,  6.25D0,  4.21D0,  4.60D0,  2.88D0,  3.41D0, $ 
                            2.37D0,  3.35D0,  2.63D0,  3.23D0,  2.60D0,  2.90D0,  2.24D0,  2.60D0, $
                            1.42D0,  1.92D0,  0.00D0,  1.84D0,  1.12D0,  1.69D0,  1.24D0,  1.86D0, $
                            0.82D0,  2.00D0,  1.00D0,  2.24D0,  1.51D0,  2.23D0,  1.12D0,  2.13D0, $
                            1.22D0,  1.55D0,  0.71D0,  1.50D0,  0.00D0,  1.00D0,  0.51D0,  1.12D0, $
                            0.33D0,  1.10D0,  0.50D0,  0.93D0,  0.13D0,  1.08D0,  0.12D0,  0.88D0, $
                            0.13D0,  0.68D0,  0.27D0,  1.45D0,  1.35D0,  1.80D0,  0.83D0,  1.09D0, $
                            0.82D0,  1.85D0,  0.71D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.12D0,  0.00D0, -0.49D0 ]
                            
  'GS98'  :  abundance = [ 12.00D0, 10.93D0,  1.10D0,  1.40D0,  2.55D0,  8.52D0,  7.92D0,  8.83D0, $
                            4.56D0,  8.08D0,  6.33D0,  7.58D0,  6.47D0,  7.55D0,  5.45D0,  7.33D0, $
                            5.50D0,  6.40D0,  5.12D0,  6.36D0,  3.17D0,  5.02D0,  4.00D0,  5.67D0, $
                            5.39D0,  7.50D0,  4.92D0,  6.25D0,  4.21D0,  4.60D0,  2.88D0,  3.41D0, $
                            2.37D0,  3.41D0,  2.63D0,  3.31D0,  2.60D0,  2.97D0,  2.24D0,  2.60D0, $
                            1.42D0,  1.92D0,  1.84D0,  1.12D0,  1.69D0,  0.94D0,  1.77D0,  1.66D0, $
                            2.00D0,  0.00D0,  1.00D0,  2.24D0,  1.51D0,  2.17D0,  1.13D0,  2.13D0, $
                            1.17D0,  1.58D0,  0.71D0,  1.50D0,  1.01D0,  0.51D0,  1.12D0, -0.10D0, $
                            1.14D0,  0.26D0,  0.93D0,  0.00D0,  1.08D0,  0.06D0,  0.88D0, -0.13D0, $
                            1.11D0,  0.28D0,  1.45D0,  1.35D0,  1.80D0,  1.01D0,  1.13D0,  0.90D0, $
                            1.95D0,  0.71D0,  0.09D0, -0.50D0 ]
                            
  'A05'   :  abundance = [ 12.00D0, 10.93D0,  1.05D0,  1.38D0,  2.70D0,  8.39D0,  7.78D0,  8.66D0, $
                            4.56D0,  7.84D0,  6.17D0,  7.53D0,  6.37D0,  7.51D0,  5.36D0,  7.14D0, $
                            5.50D0,  6.18D0,  5.08D0,  6.31D0,  3.05D0,  4.90D0,  4.00D0,  5.64D0, $
                            5.39D0,  7.45D0,  4.92D0,  6.23D0,  4.21D0,  4.60D0,  2.88D0,  3.58D0, $
                            2.30D0,  3.33D0,  2.56D0,  3.28D0,  2.60D0,  2.92D0,  2.21D0,  2.59D0, $
                            1.42D0,  1.92D0,  0.00D0,  1.84D0,  0.12D0,  1.69D0,  0.94D0,  1.77D0, $
                            1.60D0,  2.00D0,  1.00D0,  2.19D0,  1.51D0,  2.27D0,  1.07D0,  2.17D0, $
                            1.13D0,  1.58D0,  0.71D0,  1.45D0,  0.00D0,  1.01D0,  0.52D0,  1.12D0, $
                            0.28D0,  1.14D0,  0.51D0,  0.93D0,  0.08D0,  1.08D0,  0.06D0,  0.88D0, $
                           -0.17D0,  1.11D0,  0.23D0,  1.45D0,  1.38D0,  1.64D0,  1.01D0,  1.13D0, $
                            0.90D0,  2.00D0,  0.65D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $ 
                            0.00D0,  0.06D0,  0.00D0, -0.54D0 ]

  'A09'   :  abundance = [ 12.00D0, 10.93D0,  1.05D0,  1.38D0,  2.70D0,  8.43D0,  7.83D0,  8.69D0, $
                            4.56D0,  7.93D0,  6.24D0,  7.60D0,  6.45D0,  7.51D0,  5.41D0,  7.12D0, $
                            5.50D0,  6.40D0,  5.03D0,  6.34D0,  3.15D0,  4.95D0,  3.93D0,  5.64D0, $
                            5.43D0,  7.50D0,  4.99D0,  6.22D0,  4.19D0,  4.56D0,  3.04D0,  3.65D0, $
                            2.30D0,  3.34D0,  2.54D0,  3.25D0,  2.52D0,  2.87D0,  2.21D0,  2.58D0, $
                            1.46D0,  1.88D0,  0.00D0,  1.75D0,  0.91D0,  1.57D0,  0.94D0,  1.71D0, $
                            0.80D0,  2.04D0,  1.01D0,  2.18D0,  1.55D0,  2.24D0,  1.08D0,  2.18D0, $
                            1.10D0,  1.58D0,  0.72D0,  1.42D0,  0.00D0,  0.96D0,  0.52D0,  1.07D0, $
                            0.30D0,  1.10D0,  0.48D0,  0.92D0,  0.10D0,  0.84D0,  0.10D0,  0.85D0, $
                           -0.12D0,  0.85D0,  0.26D0,  1.40D0,  1.38D0,  1.62D0,  0.92D0,  1.17D0, $
                            0.90D0,  1.75D0,  0.65D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.02D0,  0.00D0, -0.54D0 ] 
                           
  'L09'   :  abundance = [ 12.00D0, 10.93D0,  3.28D0,  1.32D0,  2.81D0,  8.39D0,  7.86D0,  8.73D0, $
                            4.44D0,  8.05D0,  6.29D0,  7.54D0,  6.46D0,  7.53D0,  5.45D0,  7.16D0, $
                            5.25D0,  6.50D0,  5.11D0,  6.31D0,  3.07D0,  4.93D0,  3.99D0,  5.65D0, $
                            5.50D0,  7.46D0,  4.90D0,  6.22D0,  4.27D0,  4.65D0,  3.10D0,  3.59D0, $
                            2.32D0,  3.36D0,  2.56D0,  3.28D0,  2.38D0,  2.90D0,  2.20D0,  2.57D0, $
                            1.42D0,  1.94D0,  0.00D0,  1.78D0,  1.10D0,  1.67D0,  1.22D0,  1.73D0, $
                            0.78D0,  2.09D0,  1.03D0,  2.20D0,  1.57D0,  2.27D0,  1.10D0,  2.18D0, $
                            1.19D0,  1.60D0,  0.77D0,  1.47D0,  0.00D0,  0.96D0,  0.53D0,  1.09D0, $
                            0.34D0,  1.14D0,  0.49D0,  0.95D0,  0.14D0,  0.94D0,  0.11D0,  0.73D0, $
                           -0.14D0,  0.67D0,  0.28D0,  1.37D0,  1.36D0,  1.64D0,  0.82D0,  1.19D0, $
                            0.79D0,  2.06D0,  0.67D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
                            0.00D0,  0.08D0,  0.00D0, -0.52D0 ]

;   'T90'   :  abundance = [ 12.00D0, 11.00D0,  1.00D0,  1.15D0,  2.60D0,  8.69D0,  7.99D0,  8.91D0, $ 
;                             4.56D0,  8.00D0,  6.28D0,  7.53D0,  6.43D0,  7.50D0,  5.45D0,  7.21D0, $
;                             5.50D0,  6.58D0,  5.05D0,  6.36D0,  2.99D0,  4.88D0,  3.91D0,  5.61D0, $
;                             5.47D0,  7.46D0,  4.85D0,  6.18D0,  4.24D0,  4.60D0,  2.88D0,  3.57D0, $
;                             2.39D0,  3.35D0,  2.63D0,  3.21D0,  2.60D0,  2.93D0,  2.18D0,  2.46D0, $
;                             1.46D0,  2.10D0,  0.00D0,  1.78D0,  1.10D0,  1.69D0,  0.94D0,  1.86D0, $
;                             1.66D0,  2.00D0,  1.00D0,  2.25D0,  1.51D0,  2.19D0,  1.12D0,  2.18D0, $
;                             1.07D0,  1.58D0,  0.76D0,  1.40D0,  0.00D0,  0.88D0,  0.48D0,  1.13D0, $
;                             0.20D0,  1.07D0,  0.26D0,  0.93D0,  0.00D0,  1.08D0,  0.76D0,  0.88D0, $
;                            -0.09D0,  0.98D0,  0.26D0,  1.45D0,  1.36D0,  1.80D0,  1.13D0,  1.27D0, $
;                             0.90D0,  1.90D0,  0.71D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0, $
;                             0.00D0,  0.02D0,  0.00D0, -0.47D0 ]

ENDCASE

natoms = N_ELEMENTS(abundance)

abundance = DOUBLE(abundance)

IF KEYWORD_SET(ppm) THEN abundance = (10.D0^abundance)/(10.D0^12.)*1.d6

IF KEYWORD_SET(z) THEN abundance = abundance[z-1]  

IF KEYWORD_SET(xyz) THEN BEGIN
  aw = LOAD_ATOMIC_WEIGHTS()
  relabund = 10.D^(abundance)/TOTAL(10.D^(abundance))
  xx = relabund[0]*aw[0]/(TOTAL(relabund*aw))
  yy = relabund[1]*aw[1]/(TOTAL(relabund*aw))
  zz = TOTAL(relabund[2:natoms-1]*aw[2:natoms-1])/(TOTAL(relabund*aw))
  abundance = {x:xx, y:yy, z:zz}
ENDIF
  
RETURN, abundance
END
