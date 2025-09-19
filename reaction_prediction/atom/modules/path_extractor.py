import numpy as np
import rdkit

from rdkit import Chem

import copy

import sys
import random
import os

from reaction_prediction.atom.utils import *

class PathExtractor():

    def __init__(self, length):
        self.length = length

    def extract_numerical_path(self, smi):
        """
        Args:
            smi (str): the smiles string
                NOTE: one and only one atom has to be labeled as
                1 within the SMILES string.
        Returns:
            list: list of all paths starting from the atom (idx)
                up to a certain length.
        """
        mol = mol_with_hydrogens(smi)
        label_check = False
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                label_check = True
                init_idx = atom.GetIdx()
        assert label_check

        path = []
        P = []
        for l in range(self.length):
            path = self.extract_immediate_numerical_paths(path, mol, init_idx)
            P.extend(path)

        return P

    def extract_immediate_numerical_paths(self, q, m, init_idx=None):
        """
        NOTE: We intentionally do NOT imlpement this as a BFS.

        Args:
            q (list): is the current queue of the atoms, extracted.
            is always starts with the atom of interest and ends with
            the last atom in the extracted path.

            m (rdkit molecule): the molecule.

            init_idx: the rdkit index of the starting atom.
                NOTE: if q is empty, this has to be a valid idx.
        
        Returns:
            (list): the updated list of extrated paths.
        """
        if len(q) == 0:     # starting point to extract features
            #assert init_idx
            idx = init_idx
            atom = m.GetAtomWithIdx(idx)
            
            new_paths = []
            
            for a in atom.GetNeighbors():
                a_idx = a.GetIdx()
                bond = m.GetBondBetweenAtoms(idx, a_idx)
                bond_str = bond_to_string(bond)
                new_paths.append([idx, bond_str, a_idx])
            
            return new_paths

        else:   # expand the already existing extracted path
            N = []  # all the new path are going to be stored here
            l = len(sorted(q, key= lambda x: len(x))[-1])
            longest_paths = [i for i in q if len(i)==l]
            d = (l-1)/2.    # the actual length of the path
            for path in longest_paths:
                last_idx = path[-1]
                atom = m.GetAtomWithIdx(last_idx)
                for a in atom.GetNeighbors():
                    new_path = copy.deepcopy(path)
                    a_idx = a.GetIdx()
                    if a_idx == path[-3]:
                        continue
                    else:
                        bond = m.GetBondBetweenAtoms(a_idx, last_idx)
                        bond_str = bond_to_string(bond)
                        new_path.append(bond_str)
                        new_path.append(a_idx)
                        N.append(new_path)
            return N


    def extract_path_names(self, smi):
        """
        Args:
            smi (str): the smiles string
                NOTE: one and only one atom has to be labeled as
                1 within the SMILES string.
        Returns:
            list: list of strings. Each is the name of the path feature
            example: =(-1rO)-(00nC) which means the atom is double bonded to 
            an oxygen in ring with charge -1 is signle bonded to a neutral
            carbon which is not in ring.
        """
        mol = mol_with_hydrogens(smi)
        label_check = False
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                label_check = True
                init_idx = atom.GetIdx()
        assert label_check

        path = []
        P = []
        for l in range(self.length):
            path = self.extract_immediate_numerical_paths(path, mol, init_idx)
            path_names = prepare_path_names(mol, path)
            P.extend(path_names)

        return P

