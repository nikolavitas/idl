;+
; NAME:
;       LOAD_ATOMIC_WEIGHTS
;
; PURPOSE:
;
;       This function loads the atomic weights of the first 92 elements
;       in the periodic system.
;       
;       The data is downloaded from www.nist.gov/pml/data/comp.cfm
;       The original source is Wieser & Berglund (2009) in Pure Appl. 
;       Chem., Vol.81, No.11, p.2131-2156 (published as an IUPAC Tecnical 
;       Report:
;       http://pac.iupac.org/publications/pac/pdf/2009/pdf/8111x2131.pdf
;      
;       For more information on the dataset, check the reference. If
;       you use the data for a publication, please acknowledge it 
;       appropriately.
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
;       a = LOAD_ATOMIC_WEIGHTS([z])
;
; INPUTS:
;
;       z      = String or array, integer. The atomic number(s). If not 
;                specified, the atomic wieghts for the first 92 elements are
;                returned. 
; 
; OUTPUTS:
;
;       a     = Double array or string. Atomic weight(s).
;
; INPUT KEYWORDS:
; 
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
;       IDL> a = LOAD_ATOMIC_WEIGHTS([1, 2, 6, 7, 8])
;       IDL> PRINT, a
;       1.00790      4.00260      12.0107      14.0067      15.9994
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October 2007.
;       - small fixes and data update (NV), October 2013.
;-
;================================================================================
; LOAD_ATOMIC_WEIGHTS, IDL routine by Nikola Vitas is licensed under a Creative 
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

FUNCTION load_atomic_weights, z

aweight = [   1.0079,   4.0026,   6.9410,   9.0122,  10.8110,  12.0107,  14.0067,  15.9994, $
             18.9984,  20.1797,  22.9898,  24.3050,  26.9815,  28.0855,  30.9738,  32.0650, $
             35.4530,  39.9480,  39.0983,  40.0780,  44.9559,  47.8670,  50.9415,  51.9961, $
             54.9380,  55.8450,  58.9332,  58.6934,  63.5460,  65.3800,  69.7230,  72.6400, $
             74.9216,  78.9600,  79.9040,  83.7980,  85.4678,  87.6200,  88.9059,  91.2240, $
             92.9064,  95.9600,  98.0000, 101.0700, 102.9060, 106.4200, 107.8680, 112.4110, $
            114.8180, 118.7100, 121.7600, 127.6000, 126.9040, 131.2930, 132.9050, 137.3270, $
            138.9050, 140.1160, 140.9080, 144.2420, 145.0000, 150.3600, 151.9640, 157.2500, $
            158.9250, 162.5000, 164.9300, 167.2590, 168.9340, 173.0540, 174.9670, 178.4900, $
            180.9480, 183.8400, 186.2070, 190.2300, 192.2170, 195.0840, 196.9670, 200.5900, $
            204.3830, 207.2000, 208.9800, 209.0000, 210.0000, 222.0000, 223.0000, 226.0000, $
            227.0000, 232.0380, 231.0360, 238.0290 ]

IF KEYWORD_SET(z) THEN aweight = aweight[z-1]          

RETURN, aweight          
END