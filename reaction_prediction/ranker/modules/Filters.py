#!/usr/bin/env python
# encoding: utf-8
"""
Some simple fitlers to remove some nonsense reactions which get put up on top.

File: PathwayFilters.py
"""
import sys, os

from rpCHEM.Common.Util import molBySmiles
from openeye.oechem import *

class Filters(object):
    """A filter to throw out nonsense reactions
    """
    
    def __init__(self, filterList=None):
        """Setup default filterList
        """
        self.filterList = filterList
        if self.filterList is None:
            self.setupDefaultFilters()

    def setupDefaultFilters(self):
        """Set some reasonable defaults"""
        self.filterList = [
            BadProductFilter(),
            LoneHydrogenFilter(),
        ]
        
    def __call__(self, op):
        """Run through each of the filters, if any return false, return false"""
        for f in self.filterList:
            val = f(op)
            if not val:
                return val;
        return True

class BaseOPFilter(object):
    """Simple structure for how one of the filters should look"""
    def __init__(self):
        """Constructor"""

    def __call__(self, op):
        """Subclass must implement"""
        raise NotImplementedError()

class LoneHydrogenFilter(BaseOPFilter):
    """Return false if there is a lone hydrogen in the products"""
    
    def __call__(self, op):
        """Like it says above """
        mol = molBySmiles(op.productSmiles)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum()==1 and atom.GetDegree() == 0:
                return False

        return True

class BadProductFilter(BaseOPFilter):
    """Return false if there is a single unwanted product by formula"""

    def __init__(self):
        self.bad_product_formulae = set(['CH3+', 'CH3-', 'H+', 'H-'])
        self.failed_orbpair = []

    def __call__(self, op):
        mol_parts = op.productSmiles.split('.')
        mol = OEGraphMol()
        for part in mol_parts:
            OEParseSmiles(mol, part)
            if OEMolecularFormula(mol) in self.bad_product_formulae:
                return False
            mol.Clear()

        return True

