; This function should serve as a wrapper for other atomic data routines.
; It may include other atomic data as:
; - istopic composition
; - excitation levels
; - configuration
; - partition functions

FUNCTION load_atoms, z, abund_source = abund_source

IF ~KEYWORD_SET(z) THEN z = INDGEN(92)+1

natoms = N_ELEMENTS(z)

eion = LOAD_IONIZATION_ENERGIES(z, [0, 1])
ion1 = REFORM(eion[*, 0])
ion2 = REFORM(eion[*, 1])

label = LOAD_LIST_OF_ELEMENTS(z)

abund = LOAD_ABUNDANCES(z, source = abund_source)

w = LOAD_ATOMIC_WEIGHTS(z)

atom = {z:0, label:'', w:0.D0, ion1:0.D0, ion2:0.D0, abund:0.D0}
atoms = REPLICATE(atom, natoms)
atoms.z = z
atoms.label = label
atoms.w = w
atoms.ion1 = ion1
atoms.ion2 = ion2
atoms.abund = abund


RETURN, atoms
END
