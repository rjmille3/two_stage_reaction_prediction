#!/usr/bin/env python
import sys, os;

from openeye.oechem import OEGraphMol, OEParseSmiles;
from openeye.oechem import OEAddMols;
from openeye.oechem import OEAddExplicitHydrogens;
from openeye.oechem import OESubSearch;


from rpCHEM.Common.Util import createStandardSmiString;
from rpCHEM.Common.OrbitalModel import Orbital, orbitalIter;
from rpCHEM.Common.OrbitalModel import ORBITAL_LABEL_INCREMENT;
from rpCHEM.Common.OrbitalModel import atomHybridization;
from rpCHEM.CombiCDB.OrbitalInteraction import moveOrbitalElectrons, undoMoveOrbitalElectrons;
from rpCHEM.CombiCDB.OrbitalInteraction import compatibleOrbitalPair;
#from rpCHEM.score.orbital.ReactivityScore import ReactiveOrbitalFactory;
#from rpCHEM.score.orbital.ReactivityScore import OrbitalReactivityProbability, AcceptAllOrbitalReactivityProbability;
from rpCHEM.CombiCDB.MechanismModel import clearMechanismLabels;

from rpCHEM.CombiCDB.Util import log;



class BaseOrbitalPairFilter:
    """Simple base class for filters that analyze a proposed orbital
    pair and return True/False whether the proposal should pass the filter or not.
    """
    def __call__( self, filledOrb, unfilledOrb ):
        raise NotImplementedError();

class OverlapAtomsFilter(BaseOrbitalPairFilter):
    """Disallow any orbital pairs where the atoms in one orbital are in another orbital"""
    def __call__(self, filledOrb, unfilledOrb):
        
        return compatibleOrbitalPair(filledOrb, unfilledOrb);

class ExtremeChargeFilter(BaseOrbitalPairFilter):
    """Disallow any orbital pairs that would result in a
    product with extreme formal charges on any single atom.
    """
    
    def __call__( self, filledOrb, unfilledOrb ):
        moveOrbitalElectrons( filledOrb, unfilledOrb )
        
        violationFound = False
        for atom in filledOrb.mol.GetAtoms():
            if not self.validAtomCharge( atom ):
                # Already found a violation, don't need to look further
                violationFound = True
                break
        
        undoMoveOrbitalElectrons( filledOrb, unfilledOrb )
        
        return not violationFound

    def validAtomCharge( atom ):
        valid = True;
        valid = valid and atom.GetFormalCharge() >= ExtremeChargeFilter.minAllowedCharge( atom.GetAtomicNum() );
        valid = valid and atom.GetFormalCharge() <= ExtremeChargeFilter.maxAllowedCharge( atom.GetAtomicNum() );
        return valid;
    validAtomCharge = staticmethod(validAtomCharge);
    
    def minAllowedCharge( atomicNum ):
        return -1;
    minAllowedCharge = staticmethod(minAllowedCharge);

    def maxAllowedCharge( atomicNum ):
        # Weird, but third row elements like S and P could be written in +2 form like sulfuric acid
        if atomicNum in (15,16,12,):
            return +2;
        return +1;
    maxAllowedCharge = staticmethod(maxAllowedCharge);


class ExtremeBondFilter(BaseOrbitalPairFilter):
    """Disallow orbital pairs that create bonds with order >=4"""
    MAX_BOND_ORDER = 3
    
    def __call__( self, filledOrb, unfilledOrb ):
        moveOrbitalElectrons( filledOrb, unfilledOrb )
        
        violationFound = False
        for bond in filledOrb.mol.GetBonds():
            if not self.validBondOrder( bond ):
                # Already found a violation, don't need to look further
                violationFound = True
                break
        
        undoMoveOrbitalElectrons( filledOrb, unfilledOrb )
        
        return not violationFound;
    
    
    def validBondOrder(self, bond):
        """Test that a bond is within an acceptable range"""
        return bond.GetOrder() <= self.MAX_BOND_ORDER;
    
    

class SubstitutionHybridizationFilter(BaseOrbitalPairFilter):
    """Disallow any substitution reactions (sigma* targets)
    at sp2 or sp hybridized centers.  Should only be
    reasonable at sp3 hybridized centers.
    
    This IS reasonable in a halogenation of an alkyne.  One step here is the 
    substitution to break the halonium ion cycle at an sp2 center.
    
    Allow a special case where there is a bridged halonium.
    
    """
    def __call__( self, filledOrb, unfilledOrb ):
        if unfilledOrb.type == "sigma" and atomHybridization(unfilledOrb.atom) in (1,2):
            # Looks bad with substitution (sigma*) target hybridization, 
            #   but could be okay if this is a local orbital interaction between directly neighboring atoms
            #   or could be ok if the substitution is to break up a bridged halonium ion.
            isLocal = filledOrb.atom.GetBond( unfilledOrb.atom ) is not None;
            is3Ring = any([(unfilledOrb.atom.GetBond(atom) is not None) for atom in unfilledOrb.neighbor.GetAtoms() \
                        if atom is not unfilledOrb.atom])
            #log.info('fOrb: %s, uOrb: %s, isLocal:%s, is3Ring:%s' % \
            #    (pformat(filledOrb.toLabeledSmiAndInfoStr(10)), pformat(unfilledOrb.toLabeledSmiAndInfoStr(20)),
            #    str(isLocal), str(is3Ring)))
            if not (isLocal or is3Ring):
                # Atoms are not direct neighbors, or in a 3-ring
                # looks like a bad substitution then
                return False;
        return True;



class ReactionProposerParameters:
    """Simple struct to encapsulate parameters that could be set at the UI level."""
    
    useLocalOrbitalPairs = False;
    acceptAllOrbitals = False;
    
    def __init__(self, useLocalOrbitalPairs=False, acceptAllOrbitals=False):
        """Basic Constructor"""
        self.useLocalOrbitalPairs = useLocalOrbitalPairs;
        self.acceptAllOrbitals = acceptAllOrbitals;
