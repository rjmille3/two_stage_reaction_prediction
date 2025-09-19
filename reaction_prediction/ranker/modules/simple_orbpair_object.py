#!/usr/bin/env python
# encoding: utf-8
"""
Module to encapsulate an orbital pair for reaction prediction.

Contains convenience methods to propose reasonable reactions.

File: simple_orbpair_object.py

"""

from rpCHEM.CombiCDB.OrbitalInteraction import reactionSmilesFromOrbitalPair
from rpCHEM.CombiCDB.proposers.Util import isBondDissociation
from rpCHEM.CombiCDB.MechanismModel import ElectronArrow
from rpCHEM.CombiCDB.Const import REACTION_DELIM
from rpCHEM.Common.MolExt import clearAtomMaps
from rpCHEM.Common.Util import exact_mass

from rpCHEM.CombiCDB.proposers.PolarProposer import PolarProposer

from openeye.oechem import *  # for mass calculation added 02.16.16
from reaction_prediction.ranker.modules.Filters import Filters

from reaction_prediction.atom.utils import *

class SimpleOrbPairObject:

    # Filtering
    pass_filters = Filters()

    """Class to represent a basic Orb Pair"""
    def __init__(self, srcOrb, sinkOrb, srcAtom, sinkAtom, radical=False, pericyclic=False):
        """Constructor"""
        self.srcOrb = srcOrb
        self.sinkOrb = sinkOrb
        self.srcAtom = srcAtom
        self.sinkAtom = sinkAtom

        self.radical = radical
        self.pericyclic = pericyclic
        
        self.labelOrbPair(self.srcOrb, self.sinkOrb)
        
        self.connectedNonMappedSmiles = self.srcAtom.connectedNonMappedSmiles
        
        self.arrowObjectList = []
        self.reactionSmiles = reactionSmilesFromOrbitalPair(self.srcOrb, self.sinkOrb, False, self.arrowObjectList)
        self.reactantSmiles, self.productSmiles = self.reactionSmiles.split(REACTION_DELIM)
        self.product_masses = self.get_masses(self.productSmiles)
        self.arrowCodes = ElectronArrow.prepareArrowCodes(self.arrowObjectList)
        
    
    def __str__(self):
        """provide reaction and arrow pushing codes"""
        return str((self.reactionSmiles, self.arrowCodes))
    
    def __repr__(self):
        return "'%s'" % self.__str__()
    
    
    @classmethod
    def orbPairObjectsFromAtoms(cls, srcAtomList, sinkAtomList, photo=False, radical=False, pericyclic=False, checkIntra=True):
        """Given some atoms, run the proposer and make a list of SimpleOrbPairObjects"""
        proposer = None

        proposer = PolarProposer(srcAtomList, sinkAtomList, checkIntra=checkIntra)
        
        orbPairObjectList = []
        

        for srcOrb, sinkOrb, srcAtom, sinkAtom in proposer:

            new_op = cls(srcOrb, sinkOrb, srcAtom, sinkAtom, radical, pericyclic)
            if cls.pass_filters(new_op):
                orbPairObjectList.append(new_op)

        return orbPairObjectList
    
    @classmethod
    def orbPairObjectsFromAtoms_bounded(cls, srcAtomList, sinkAtomList, max_orbs=50, photo=False, radical=False, pericyclic=False, checkIntra=True):
        
        """Given some atoms, run the proposer and make a list of SimpleOrbPairObjects"""

        proposer = PolarProposer(srcAtomList, sinkAtomList, checkIntra=checkIntra)      

        orbPairObjectList = []
        i = 0

        for srcOrb, sinkOrb, srcAtom, sinkAtom in proposer:

            new_op = cls(srcOrb, sinkOrb, srcAtom, sinkAtom, radical, pericyclic)
            if cls.pass_filters(new_op):
                orbPairObjectList.append(new_op)

            i += 1
            if i > max_orbs-1:
                break

        return orbPairObjectList
    
    @staticmethod
    def labelOrbPair(src, sink):
        """Convenience to turn an orb pair into a string tuple"""
        clearAtomMaps(src.mol)
        if not isBondDissociation(src, sink):
            src.labelOrbitalAtoms(10)   
        sink.labelOrbitalAtoms(20)

    def get_masses(self, smi_str):
        mass_list = []
        smi_list = smi_str.split('.')
        mol = OEGraphMol()
        for smi in smi_list:
            OEParseSmiles(mol, smi)
            mass = exact_mass(mol)
            mass_list.append( round(mass, 1) )
            mol.Clear()

        return mass_list


