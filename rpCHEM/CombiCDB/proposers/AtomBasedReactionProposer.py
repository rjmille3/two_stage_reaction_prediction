#!/usr/bin/env python
# encoding: utf-8
"""
Class to make orbital pair reaction proposals based on a restricted number of atoms.

AtomBasedReactionProposer.py

"""

from pprint import pformat;
from openeye.oechem import OEDetermineComponents

from rpCHEM.CombiCDB.proposers.Util import log;

from rpCHEM.Common.OrbitalModel import Orbital
from rpCHEM.Common.MolExt import atomLocalKekuleIterator, kekulizeMol
from rpCHEM.Common.CanonicalAtomMapSmiles import createCanonicalAtomMapSmiString

## Working with SimpleAtomObjects
from rpCHEM.ML.orbDB.simplepred.SimpleAtomObject import SimpleAtomObject

from rpCHEM.CombiCDB.proposers.ResonanceStructUtil import ResonanceStructUtil
from rpCHEM.CombiCDB.proposers.Util import compatibleOrbitalPair, labelOrbPair

class AtomBasedReactionProposer(object):
    """Propose orbital interactions based on restricted atoms.
    
    Assumes we are given two lists of SimpleAtomObjects, one for potential source, one for sinks
    """
    def __init__(self, srcAtomList, sinkAtomList, orbitalPairFilters=None, checkResonance=True, checkIntra=True, checkRadical=False):
        """Constructor"""
        self.srcAtomList = srcAtomList
        self.sinkAtomList = sinkAtomList
        self.orbitalPairFilters = orbitalPairFilters
        self.mol = None
        self.seenOrbPairSet = set([])
        self.componentParts = None
        self.checkRadical = checkRadical
        self.checkIntra = checkIntra
        self.interOnly = False

        
        if len(self.srcAtomList) > 0:
            #self.mol = self.srcAtomList[0].mol
            self.mol = list(self.srcAtomList)[0].mol
            self.numParts, self.partsMap = OEDetermineComponents(self.mol)
            self.mol = kekulizeMol(self.mol)

            # big errors in original code here. used to be called "intraOnly" which meant the opposite.
            # set interOnly to false, meaning we want intra considered no matter how many reactants.
            #self.interOnly = self.numParts > 1 and self.checkIntra
            self.interOnly = False
        
        #self.checkResonance = checkResonance
        #self.resonanceStructUtil = None
        #if self.checkResonance and self.mol is not None:
        #    self.resonanceStructUtil = ResonanceStructUtil(self.mol, checkRadical=checkRadical)
        
        self.srcOrbFactory = None
        self.sinkOrbFactory = None
        self.setOrbitalFactories()
        
        #self.inter = inter

    def setOrbitalFactories(self):
        """Abstract - Must implement"""
        raise NotImplementedError()
        
    def passOrbitalPairFilters(self, src, sink):
        """Check if any of the orbital pair filters do NOT pass"""
        for filterFunc in self.orbitalPairFilters:
            if not filterFunc(src, sink):
                return False

        #if self.checkResonance:
        #    if not self.resonanceStructUtil.checkOrbPair(src, sink):
        #        return False
        
        strRep = labelOrbPair(src, sink)
        if strRep in self.seenOrbPairSet:
            return False
        
        self.seenOrbPairSet.add(strRep)
        return True
    
    def proposeOrbitalPairs(self):
        """Main work of proposing.  Loop over the factories and make inter/intra orb pairs.
        
        NOTE: To get full coverage of the kekule structures when we have potential intra-molecular 
        reactions within the same aromatic system, we need to NOT pre-compute the sinkOrbList.  
        """
        
        ## Pre-make the sinkOrbs
        #sinkOrbList = [orb for orb in self.sinkOrbFactory]
        
        for srcOrb, srcAtom in self.srcOrbFactory:
            for sinkOrb, sinkAtom in self.sinkOrbFactory: #sinkOrbList:

                #print("sinkOrb ", sinkOrb)
                #print("sinkAtom ", sinkAtom)
                #log.info('srcOrb: %s' % pformat(srcOrb.toLabeledSmiAndInfoStr(10)))
                #log.info('sinkOrb: %s' % pformat(sinkOrb.toLabeledSmiAndInfoStr(20)))

                ## First do the inter check.  May want to skip if not from different components
                if self.interOnly and (self.partsMap[srcOrb.atom.GetIdx()] == self.partsMap[sinkOrb.atom.GetIdx()]):
                    #log.info('Skipping because NOT INTER')
                    continue
                

                #print("compatible orb pair ", compatibleOrbitalPair(srcOrb, sinkOrb))

                # First intra:
                #if not self.inter:
                if compatibleOrbitalPair(srcOrb, sinkOrb):

                    #print("compatible srcOrb was ", srcOrb)
                    #print("compatible sinkOrb was ", sinkOrb)
                    #log.info('About to return an intra mol orb pair: %s, %s' %
                    # (pformat(srcOrb.toLabeledSmiAndInfoStr(10)),
                    # pformat(sinkOrb.toLabeledSmiAndInfoStr(20))))
                    #log.info('Src Orb.mol: %s' % pformat(srcOrb.mol))
                    #log.info('sinkOrb.mol: %s' % pformat(sinkOrb.mol))
                    ((clonedSrcOrb, clonedSinkOrb), copyMol) = \
                        Orbital.cloneOrbitalsWithCommonMol([srcOrb, sinkOrb])
                    #log.info('Returning an intra mol orb pair: %s, %s' % \
                    # (pformat(clonedSrcOrb.toLabeledSmiAndInfoStr(10)),
                    # pformat(clonedSinkOrb.toLabeledSmiAndInfoStr(20))))
                    yield (clonedSrcOrb, clonedSinkOrb, srcAtom, sinkAtom)
           
                    #NOTE: this appears to control whether duplicate copies of species are used for src/sink
                    # Then inter:
                    ## Only do this if from same mol
               # else:
                    #if self.partsMap[srcOrb.atom.GetIdx()] == self.partsMap[sinkOrb.atom.GetIdx()]:
                    #( (compositeSrcOrb, compositeSinkOrb), compositeMol ) =  \
                    #                Orbital.compositeOrbitalsFromDistinctMols( [srcOrb, sinkOrb] );
                    #yield (compositeSrcOrb, compositeSinkOrb, srcAtom, sinkAtom)
        
    def __iter__(self):
        """Iterator to yield the orbitals.  Propose and check if we pass orbital filters"""
        for srcOrb, sinkOrb, srcAtom, sinkAtom in self.proposeOrbitalPairs():
            #print("in iter loop ", srcOrb)
            yield (srcOrb, sinkOrb, srcAtom, sinkAtom)
            #removed filtering due to step "C[O+:11](C)[Na:10]>>C[O:11]C.[Na+:10] 10,11=11". when providing correct src and sink orbital, it does not pass filters
            #if self.passOrbitalPairFilters(srcOrb, sinkOrb):
            #    yield (srcOrb, sinkOrb, srcAtom, sinkAtom)
        

class AtomBasedOrbitalFactory(object):
    """Abstract Virtual class defining the expected interface.

    Handles how to run the atomIterator.
    User must implement iter 
    """
    def __init__(self, atomObjList, allSymmetry=True, allKekule=True):
        """Basic Constructor"""
        self.atomObjList = atomObjList
        self.allSymmetry = allSymmetry
        self.allKekule = allKekule

    def atomIterator(self):
        
        """Convenience to handle how we iterate over the atoms in the """
        for atomObj in self.atomObjList:
            atomList = [atomObj.atom]
            if self.allSymmetry:
                atomList = atomObj.symmetricAtoms
            
            for oeAtom in atomList:
                # If we need to iterate over the kekule structs, do so:
                mol = atomObj.mol
                #log.info('Starting with mol with smiles: %s' % (createCanonicalAtomMapSmiString(mol)))
                if self.allKekule:
                    for iMol, nMol in enumerate(atomLocalKekuleIterator(mol, oeAtom)):
                        oeAtom.SetMapIdx(0)
                        yield oeAtom, atomObj
                else:
                    yield oeAtom, atomObj
    
    def __iter__(self):
        """Needs to be implemented.

        Expected to return the orbital and the current atomObj
        """
        raise NotImplementedError()
