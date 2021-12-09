;+
; NAME:
;
;   MOLECULES
;
; PURPOSE:
;
;   This is a wrapper procedure for molecular data. There is an internal catalogue
;   with 322 DIATOMIC and 27 MULTI-NUCLEAR molecules. The catalogue (an IDL 
;   structure) is saved in catalogue_of_molecules.sav. The structure contains the 
;   following tags for every molecule:
;
;     NAME            STRING    Chemical formula, e.g. 'H2', 'CO', 'H2+'
;     NC              INT       Number of constituents, always 2 for diatomics
;     OMEGA           INT       Array of stoichiometric coefficients
;     CHARGE          INT       Charge, e.g. -1, 0, +1
;     RH_CODE         STRING    Code as in the RH code (Uitenbroek)
;     EION            FLOAT     Ionization energy of constituent with the lowest Eion
;     D0              FLOAT     Dissociation energy in eV
;     SOURCE          STRING    List of available data sets for this molecule
;     TYPE            STRING    Type of data, pf for partition function, kp for
;                               chemical equilibrium constant, eint for internal
;                               energy
;     REFERENCE       STRING    List of ADS references (for each dataset) 
;
;   Input to the routine is the temperature (in Kelvin) and the chemical
;   formula of the molecule (e.g. 'H2', 'CO'). The output is specified
;   by the keywords DATA and TYPE. Depending on TYPE, the output can be
;   the partition function ('pf'), the chemical equilibrium constant for pressure
;   ('kp') or the internal energy due to the rotational-vibrational bound
;   states ('eint'). The chemical equilibrium constant is given for the
;   pressures (to get its form for the number densities divide it by kT. 
;   The internal energy in erg per particle. (All the units in the output
;   and used internally in the code are in CGS). 
; 
;   The code does not compute these functions but employs various fits
;   available in the literature and evaluates them for the given temperature
;   and the specified molecule. The keyword DATA specifies the source from
;   which the data is taken and the particular fit. The following sources
;   are available:
;
;   Vardya 1961 (Kp, only H2 and H2+), 'vardya'
;   Vardya 1965 (Eint, only H2 and H2+), 'vardya'
;   Tsuji 1973 (U, 95 diatomics), 'tsuji1973'
;   Sauval and Tatum 1984 (U, Kp, D0, 294 diatomics), 'sauvalandtatum1984'
;   Irwin 1981 (U, only H2 and CO), 'irwin1981'
;   Bohn and Wolf 1984 (U, only H2 and CO), 'bohn1984)
;   Rossi and Maciel 1983 (53 molecules, some multi-nuc), 'rossi_and_maciel1983'
;
;   None of the datasets contain all three functions. The missing functions 
;   when there is enough data the code computes using the following equations:
;  
;   Kp = (2 pi mab kT/h^2)^1.5 exp(-D0/kT) Ua Ub/Uab
;   
;   Eint = kT d ln U / d ln T
; 
;   where Ua and Ub are the partition functions of the atoms A and B,
;   Uab is the partition function of the diatomic, mab is the reduced
;   mass (mab = ma*mb/(ma+mb)*m_h) and D0 is the dissociation energy.
;
;   The default data is sauvalandtatum1984 for all the molecules that
;   are included in their list and tsuji1973 for the rest.
;
;   The dissociation constants are available only for the molecules
;   of Sauval and Tatum. Therefore the chemical equilibrium constant
;   cannot be computed for the molecules that are exclusiv for Tsuji
;   if one does not supply D0 manually. 
;
;   Note: The set of Irwin 1981 for the diatomics replicates an earlier
;   set of Tatum 1966 that is already included and updated in Sauval and
;   Tatum 1984. The set of Kurucz can be easily added to the catalogue.
;   There are also updated sets of Sauval and Irwin (priv.comm), but 
;   I do not include them here as I miss the full reference and the
;   data description.
;
;   D0 values are taken from Barklem 2016 for all molecules present in
;   his list.
;
;   OUTDATED:
;   Note that each of these datasources strictly applies only to a limited
;   range of temperatures: Gray (2520 - 25200 K), Irwin (1000 - 16000 K), 
;   Sauval and Tatum (1000 - 9000 K). Out of these ranges the values returned 
;   by this function may not be reliable and it is up to the user to decide
;   on the extrapolation procedure in that case. If keyword RANGE is set in
;   a call to a partition function, then the temperatures out of the range
;   are clipped to it.
;
;   If you use any of the data, please acknowledge it by citing the 
;   corresponding paper.
;
;   Installation
;
;   OUTDATED
;   Before running the code for the first time, copy all the functions
;   into a directory within your IDL path. Open diatomics_path.pro and
;   set the path to that exact directory. Make sure that you already
;   have my partition_function routines downloaded and installed, see
;   http://nikolavitas.blogspot.com.es/2013/12/atomic-partition-functions-in-idl.html 

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
;   f = MOLECULES(t, name, data = data, type = type, d0 = d0, mab = mab)
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
;   data = String, scalar. Reference to the data set that has to be used. 
;          The options are: irwin1987, sauvalandtatum1984, vardya, tsuji1973
;          and bohn1984. Only the first letter matters!
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
;   Libraries for the atomic data and the atomic partition functions.
;
; EXAMPLE:
;
;   t = 2000. + DINDGEN(6)*2000.
;   name = 'H2'
;
;   Chemical equilibrium constant
;   kps = MOLECULES(t, name, data = 's', type = 'kp')
;   kpv = MOLECULES(t, name, data = 'v', type = 'kp')
;   kpt = MOLECULES(t, name, data = 't', type = 'kp')
;   kpi = MOLECULES(t, name, data = 'i', type = 'kp')
;   kpb = MOLECULES(t, name, data = 'b', type = 'kp')
;   
;   Partition function
;   pfs = MOLECULES(t, name, data = 's', type = 'pf')
;   pfi = MOLECULES(t, name, data = 'i', type = 'pf')
;   pfb = MOLECULES(t, name, data = 'b', type = 'pf')
;   
;   Energy
;   es = MOLECULES(t, name, data = 's', type = 'eint')
;   ei = MOLECULES(t, name, data = 'i', type = 'eint')
;   eb = MOLECULES(t, name, data = 'b', type = 'eint')
;   ev = MOLECULES(t, name, data = 'v', type = 'eint')
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (February 2014)
;
;-
;================================================================================
; MOLECULES, IDL routine by Nikola Vitas is licensed under a Creative Commons 
; Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The routine called by MOLECULES use data from different sources. The references
; for the are in their headers. If you use any of them in a publication, please
; acknowledge by citing the corresponding paper.
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

FUNCTION molecules, t, name, data = data, type = type, d0 = d0, mab = mab

@molecules_path
filename = path + 'catalogue_of_molecules.sav'
RESTORE, filename

ind = WHERE(name EQ cat.name)

IF ind EQ -1 THEN MESSAGE, 'Data for ' + name + ' is not available.'

d0 = cat[ind].d0
mab = cat[ind].mab

letter = STRLOWCASE(STRMID(data, 0, 1))

CASE letter OF 
  'i' : f = DIATOMICS_IRWIN(t, name, type = type, d0 = d0, mab = mab)
  's' : f = DIATOMICS_SAUVAL_AND_TATUM(t, name, type = type, d0 = d0, mab = mab)
  'v' : f = DIATOMICS_VARDYA(t, name, type = type, d0 = d0, mab = mab)
  'b' : f = DIATOMICS_BOHN(t, name, type = type, d0 = d0, mab = mab)
  't' : f = DIATOMICS_TSUJI(t, name, type = type, d0 = d0, mab = mab)    
  'r' : f = MOLECULES_ROSSI_AND_MACIEL(t, name, type = type, d0 = d0, mab = mab)    
ENDCASE  
  
RETURN, f
END
