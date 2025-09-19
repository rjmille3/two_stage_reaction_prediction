#!/usr/bin/env python
# encoding: utf-8
"""
Util.py

Created by Matt Kayala
"""
import Const
import sys, os
import logging
import pickle

log = logging.getLogger(Const.APPLICATION_NAME)
log.setLevel(Const.LOGGER_LEVEL)

handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter(Const.LOGGER_FORMAT)

handler.setFormatter(formatter)
log.addHandler(handler)




## Some actual code utils
from rpCHEM.Common.Util import splitCompositeSmilesToList, molBySmiles
from openeye.oechem import OEGraphMol, OEParseSmiles, OECreateIsoSmiString
from openeye.oechem import OEAssignAromaticFlags, OEClearAromaticFlags, OEKekulize
from openeye.oechem import OECanonicalOrderAtoms, OECanonicalOrderBonds
from rpCHEM.Common.Util import standardizeSmiles, createAtomMapSmiString
from rpCHEM.Common.Util import createStdAtomMapSmiString
from rpCHEM.Common.CanonicalAtomMapSmiles import (canonicalizeAtomMapSmiString,
                                                createCanonicalAtomMapSmiString
                                                )

from rpCHEM.Common.MolExt import clearAtomMaps, removeNonsenseStereo

def canonicalKekule(smi):
    """Simple method to canonicalize the kekule form"""
    mol = molBySmiles(smi);
    removeNonsenseStereo(mol)
    OEAssignAromaticFlags(mol)
    for bond in mol.GetBonds():
        if bond.IsAromatic():
            bond.SetIntType(5)
    OECanonicalOrderAtoms(mol)
    OECanonicalOrderBonds(mol)
    OEClearAromaticFlags(mol)
    OEKekulize(mol)
    return createCanonicalAtomMapSmiString(mol)



def loadPickledObject(filename):
    """Convenience to load up a single simple object from a pickled file"""
    ifs = open(filename)
    obj = pickle.load(ifs)
    ifs.close();
    return obj;
