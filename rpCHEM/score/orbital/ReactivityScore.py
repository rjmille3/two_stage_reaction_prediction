"""Qualitative Reactivity Estimators, not Meant to directly Calculate TS Energy.
Just give some numerical measure that correlates to the likelihood that
different orbitals are going to react.

Key functions for external use:

    ReactiveOrbitalFactory:
        Given an input molecule, will output all of the orbitals on the molecule
        that could be worth considering for reactivity analysis.
s
"""
#from sets import Set;
#from openeye.oechem import OEGetHybridization, OEHybridization_sp, OEHybridization_sp2, OEHybridization_sp3;
from openeye.oechem import OECreateIsoSmiString, OEAddExplicitHydrogens, OEGetPathLength, OEPerceiveSymmetry;
from openeye.oechem import OEAssignAromaticFlags, OEClearAromaticFlags;




class ReactiveOrbitalFactory:
    """
    Given an input molecule, will output all of the orbitals on the molecule
    that could be worth considering for reactivity analysis based on the
    results of the provided orbital reactivity calculator.
    """
    orbitalReactivityCalc = None;
    
    def __init__(self, orbitalReactivityCalc, restrictSymmetry=True):
        self.orbitalReactivityCalc = orbitalReactivityCalc;
        self.restrictSymmetry = restrictSymmetry;



    @staticmethod
    def prepareChainOrbital(orbital, parentOrb, occupied):
        """Prepare an extended orbital chain combination of the 
        given orbital and parent orbital candidates.
        
        Determine whether the given orbital is a potential candidate
        for use as an extended orbital chain combined with the parent orbital,
        with different consideration depending on whether we're considering 
        occupied or unoccupied orbitals.
        
        If the orbital is a valid candidate, return the extended orbital chain,
        otherwise return None to indicate an invalid candidate.
        """
        chainOrbital = None;
        
        hasExtendedPotential = (not occupied and orbital.neighbor is not None and orbital.neighbor.GetAtomicNum() == 1); # Adjacent proton elimination candidates (E1, E2, enolate)
        #hasExtendedPotential = (not occupied and ???); # Retro-aldol, Grob fragmentations, etc. Sigma[C-C] bonds can be extended candidates, if there is a respective source filled orbital directly attached.
        #                                                   Incorporate into orbital pair proposal for special cases, like 1,2 (carbocation) rearrangements where C-C / C-H sigma bond could be source, but only for directly adjacent sinks
        hasExtendedPotential = hasExtendedPotential or orbital.type == "pi";    # Pi system conjugation is most common type

        if hasExtendedPotential:
            # Ensure that none of the proposed orbital atoms overlap
            parentOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet( parentOrb );
            candidateOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet( orbital );

            if len(parentOrbIdxSet.intersection(candidateOrbIdxSet)) > 0:
                # Orbitals have overlapping atoms, does not make sense for extended chain?
                #   Beware of confounding cases with pericyclic reactions and multi-bond forming cycles
                #   (e.g., bromonium ion formation).
                pass;
            else:
                # Orbital is a good candidate and no overlapping atoms with parent, 
                #   should be able to prepare a potential extended orbital candidate.
                chainOrbital = orbital.flippedBondOrbital();
                chainOrbital.extOrbital = parentOrb;

        return chainOrbital;


    @staticmethod
    def orbitalAtomIdxSet(orbital, atomIdxSet=None):
        """Update the atomIdxSet with all of the atoms covered by the given orbital, incl. extensions."""
        if atomIdxSet is None:
            atomIdxSet = set()
        atomIdxSet.add(orbital.atom.GetIdx())
        if orbital.neighbor is not None:
            atomIdxSet.add(orbital.neighbor.GetIdx())
        if getattr(orbital, "extOrbital", None) is not None:
            # Recursively check extended orbitals
            ReactiveOrbitalFactory.orbitalAtomIdxSet(orbital.extOrbital, atomIdxSet)
        return atomIdxSet

