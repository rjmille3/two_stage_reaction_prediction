#!/usr/bin/env python
# encoding: utf-8
"""
Class to XX, callable from Command Line.

SimpleAtomObject.py

Created by Matt Kayala
"""

import sys
import os
from optparse import OptionParser;

from rpCHEM.ML.orbDB.simplepred.Util import log, canonicalKekule

from rpCHEM.Common.Util import splitCompositeSmilesToList, molBySmiles
from rpCHEM.Common.Util import joinSmilesListToCompositeSmiles
from openeye.oechem import OEGraphMol, OEParseSmiles, OECreateIsoSmiString
from openeye.oechem import OEAssignAromaticFlags, OEClearAromaticFlags, OEKekulize
from openeye.oechem import OEPerceiveChiral, OEPerceiveSymmetry, OESubSearch
from openeye.oechem import OEFindRingAtomsAndBonds
from rpCHEM.Common.Util import standardizeSmiles, createAtomMapSmiString
from rpCHEM.Common.Util import createStdAtomMapSmiString, clearAtomMapsSmiStr
from rpCHEM.Common.CanonicalAtomMapSmiles import (canonicalizeAtomMapSmiString,
                                                createCanonicalAtomMapSmiString
                                                )
from rpCHEM.Common.MolExt import clearAtomMaps, removeNonsenseStereo
from rpCHEM.Common.MolExt import setSingleExplicitHydrogens, getAtmMappedSmilesFromCompositeSmiles


class SimpleAtomObject(object):
    """Simple object to capture a representation of an atom.
    
    Basically, given an atom object, store a copy of the atom object, full smiles, and the smiles of 
    only connected component containing self.
    """
    def __init__(self, oeAtom, oeMol, symmetricAtoms):
        """Set the basic params"""
        self.atom = oeAtom
        self.mol = oeMol
        
        self.symmetricAtoms = symmetricAtoms;
        #self.symmetryClass = symmetricAtoms;
        
        self.fullSmiles = None
        self.connectedSmiles = None
        self.connectedNonMappedSmiles = None
        
        
        clearAtomMaps(self.mol)
        self.atom.SetMapIdx(1)
        self.fullSmiles = createCanonicalAtomMapSmiString(self.mol)
        self.atom.SetMapIdx(0)
        self.connectedSmiles = getAtmMappedSmilesFromCompositeSmiles(self.fullSmiles, clearMaps=False)
        self.connectedSmiles = canonicalizeAtomMapSmiString(self.connectedSmiles)
        self.connectedNonMappedSmiles = clearAtomMapsSmiStr(self.connectedSmiles);
        self.connectedKekuleSmiles = canonicalKekule(self.connectedSmiles)
        
        # Eventual ML stuff
        self.rawFeatDict = None
        self.dbIdFeatDict = None
        self.normFeatDict = None
        self.fillPredValue = None
        self.unfillPredValue = None
        self.fillPredDecision = None
        self.unfillPredDecision = None
        
    
    def __str__(self):
        """Just set this to return the connectedSmiles"""
        return self.connectedSmiles
    
    def __repr__(self):
        return "'%s'" % self.__str__();
    
    @staticmethod
    def atomObjFromReactantSmi(reactantSmiles):
        """Convenience method to make a list of SimpleAtomObjects"""
        """Given a single reactant (connected component), return the list of atoms."""
        # First, ensure a single copy!
        #reactantSmiles = canonicalKekule(reactantSmiles);
        
        reactantSmiles = clearAtomMapsSmiStr(reactantSmiles)
        mol = molBySmiles(reactantSmiles)
        OEAssignAromaticFlags(mol)
        reactantSmiles = createCanonicalAtomMapSmiString(mol)
        smilesSet = set(splitCompositeSmilesToList(reactantSmiles))
        newReactantSmiles = joinSmilesListToCompositeSmiles(list(smilesSet))
        newReactantSmiles = canonicalizeAtomMapSmiString(newReactantSmiles)
        
        mol = molBySmiles(newReactantSmiles);
        setSingleExplicitHydrogens(mol)
        OEAssignAromaticFlags(mol)
        clearAtomMaps(mol)
        OEPerceiveChiral(mol)
        removeNonsenseStereo(mol)
        OEPerceiveSymmetry(mol)
        OEFindRingAtomsAndBonds(mol)
        
        resBySymmetryClass = {}
        seenSymClassSet = set([])
        for atm in mol.GetAtoms():
            theSymClass = atm.GetSymmetryClass()
            if theSymClass not in resBySymmetryClass:
                resBySymmetryClass[theSymClass] = SimpleAtomObject(atm, mol, [atm])
            else:
                resBySymmetryClass[theSymClass].symmetricAtoms.append(atm);
        return resBySymmetryClass.values();
    
    
    @staticmethod
    def combineAtomObjListToSmiles(atomObjList):
        """Given some 'filtered' atomObjs in a list, return a composite smiles with all of these atoms labeled."""
        if len(atomObjList) == 0:
            return None
        mol = atomObjList[0].mol
        clearAtomMaps(mol)
        for atomObj in atomObjList:
            [atm.SetMapIdx(1) for atm in atomObj.symmetricAtoms]
            #atomObj.atom.SetMapIdx(1)
        smi = createStdAtomMapSmiString(mol)
        clearAtomMaps(mol)
        return smi;
