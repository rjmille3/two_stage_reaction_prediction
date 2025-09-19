#!/usr/bin/env python
# encoding: utf-8
"""
A convenience module with some code to help parse marvin xml data into
smiles and arrow codes.

Include some code to make the 

Also, will include some code to map to orbital strings
"""
from urllib.parse import urlencode
from urllib.request import urlopen
from cgi import parse_qsl

from rpCHEM.Common.Util import molBySmiles
from rpCHEM.CombiCDB.ArrowCodesToOrbitals import orbitalPairFromArrowStr, ArrowConversionError
from rpCHEM.CombiCDB.ArrowCodesToOrbitals import BondDoesNotExistException

def arrowToOrbital(smiles, arrowCodes):
    """Convert arrow representations into orbitals format"""
    mol = molBySmiles(smiles)
    orbInfoTuple = orbitalPairFromArrowStr(mol, arrowCodes)

    srcOrb, sinkOrb, nElectrons = orbInfoTuple

    srcInfoStr = srcOrb.formatInfoStr()
    sinkInfoStr = sinkOrb.formatInfoStr()

    return srcInfoStr, sinkInfoStr, srcOrb, sinkOrb, nElectrons



