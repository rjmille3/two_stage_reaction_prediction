#!/usr/bin/env python
# encoding: utf-8
"""
Class to encapsulate a pair of orbital objects with associated SimpleAtomObjects.

Used to track things like the feature dicts and the associated ranking

SimpleOrbPairObject.py

Created by Matt Kayala
"""

import sys
import os
from optparse import OptionParser;

from rpCHEM.CombiCDB.OrbitalInteraction import reactionSmilesFromOrbitalPair, ElementaryStepData
from rpCHEM.CombiCDB.OrbitalInteraction import moveOrbitalElectrons, undoMoveOrbitalElectrons
from rpCHEM.CombiCDB.OrbitalInteraction import ElementaryStepData
from rpCHEM.CombiCDB.proposers.PolarProposer import PolarProposer
from rpCHEM.CombiCDB.proposers.Util import isBondDissociation
from rpCHEM.CombiCDB.MechanismModel import ElectronArrow
from rpCHEM.CombiCDB.Const import REACTION_DELIM
from rpCHEM.Common.OrbitalModel import Orbital
from rpCHEM.Common.Util import createStdAtomMapSmiString
from rpCHEM.Common.MolExt import clearAtomMaps


from Util import log;

class SimpleOrbPairObject:
    """Class to represent a basic Orb Pair"""
    def __init__(self, srcOrb, sinkOrb, srcAtom, sinkAtom):
        """Constructor"""
        self.srcOrb = srcOrb
        self.sinkOrb = sinkOrb
        self.srcAtom = srcAtom
        self.sinkAtom = sinkAtom
        
        self.labelOrbPair(self.srcOrb, self.sinkOrb)
        
        self.connectedNonMappedSmiles = self.srcAtom.connectedNonMappedSmiles;
        
        self.arrowObjectList = [];
        self.reactionSmiles = reactionSmilesFromOrbitalPair(self.srcOrb, self.sinkOrb, False, self.arrowObjectList)
        self.reactantSmiles, self.productSmiles = self.reactionSmiles.split(REACTION_DELIM)
        self.arrowCodes = ElectronArrow.prepareArrowCodes(self.arrowObjectList)
        
        # For eventual ML
        self.rawFeatDict = None
        self.normFeatDict = None
        self.predValue = None
        self.rawPredValues = None
    
    
    def toElementaryStepData(self):
        """Make an elementary step object version of self"""
        eStep = ElementaryStepData();
        eStep['compositeMol'] = self.srcOrb.mol
        eStep['filledOrb'] = self.srcOrb
        eStep['unfilledOrb'] = self.sinkOrb
        return eStep;
    
    
    def __str__(self):
        """Show something"""
        return str((self.reactionSmiles, self.arrowCodes))
    
    def __repr__(self):
        return "'%s'" % self.__str__();
    
    @staticmethod
    def orbPairObjectsFromAtoms(srcAtomList, sinkAtomList):
        """Given some atoms, run the proposer and make a list of SimpleOrbPairObjects"""
        proposer = PolarProposer(srcAtomList, sinkAtomList);
        orbPairObjectList = []
        for srcOrb, sinkOrb, srcAtom, sinkAtom in proposer:
            orbPairObjectList.append(SimpleOrbPairObject(srcOrb, sinkOrb, srcAtom, sinkAtom))
        return orbPairObjectList;
    
    @staticmethod
    def labelOrbPair(src, sink):
        """Convenience to turn an orb pair into a string tuple"""
        clearAtomMaps(src.mol)
        if not isBondDissociation(src, sink):
            src.labelOrbitalAtoms(10)   
        sink.labelOrbitalAtoms(20)

    
