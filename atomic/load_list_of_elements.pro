;+
; NAME:
;       LOAD_LIST_OF_ELEMENTS
;
; PURPOSE:
;
;       This function returns the elemental symbols for the first 92 elements
;       or for a given subset of atomic numbers. If keyword element is set
;       to an undefined variable, the function loads it with the names of the 
;       elements.
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
;       symbols = LOAD_LIST_OF_ELEMENTS([z, element = element])
;
; INPUTS:
;
;       z     = Scalar or array, integer. Atomic number(s).
; 
; OUTPUTS:
;
;       symbols = String array, Symbols of the elements.
;
; INPUT KEYWORDS:
;
;      printout = True/false. If set, the list of elemens is printed.
;
; OUTPUT KEYWORDS:
;
;       element = String array. Names of the elements
;
; EXAMPLE:
;
;       IDL> sym = LOAD_LIST_OF_ELEMENTS([6, 7, 8], element = names)
;       IDL> PRINT, sym
;       C N O
;       IDL> PRINT, names
;       Carbon Nitrogen Oxygen
;
; DEPENDENCIES:
;
; MODIFICATION HISTORY:
;
;       Written by Nikola Vitas, October 2007.
;       - small fixes (NV), October 2013.
;-
;================================================================================
; LOAD_LIST_OF_ELEMENTS, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
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

FUNCTION load_list_of_elements, z, element = element, printout = printout

COMPILE_OPT idl2, HIDDEN

element = ['Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon', $
           'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium', 'Magnesium', $
           'Aluminum', 'Silicon', 'Phosphorus', 'Sulfur', 'Chlorine', 'Argon', $
           'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium', $
           'Manganese', 'Iron', 'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', $
           'Germanium', 'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium', $
           'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum', 'Technetium', $
           'Ruthenium', 'Rhodium', 'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin', $
           'Antimony', 'Tellurium', 'Iodine', 'Xenon', 'Cesium', 'Barium', 'Lanthanum', $
           'Cerium', 'Praseodymium', 'Neodymium', 'Promethium', 'Samarium', 'Europium', $
           'Gadolinium', 'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium', $
           'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium', $
           'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury', 'Thallium', 'Lead', 'Bismuth', $
           'Polonium', 'Astatine', 'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium', $
           'Protactinium', 'Uranium']
IF KEYWORD_SET(z) THEN element = element[z-1]


symbol = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', $
          'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',  $
          'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', $
          'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', $
          'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', $
          'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', $
          'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U']

IF KEYWORD_SET(z) THEN symbol = symbol[z-1]

IF KEYWORD_SET(printout) THEN BEGIN
  PRINT, ' Z  Symbol            Name'
  PRINT, '--------------------------'
  FOR iz = 0, 91 DO $
    PRINT, iz+1, ' ', symbol[iz], ' ', element[iz], format = '(I3, A1, A6, A1, A15)'
ENDIF
    
RETURN, symbol

END