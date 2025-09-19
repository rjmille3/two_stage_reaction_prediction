import numpy as np
import rdkit

from rdkit import Chem

from reaction_prediction.atom.utils import *
from reaction_prediction.atom.modules.path_extractor import PathExtractor

import sys
import os
import argparse
import json
import csv

class FeatSel():
    
    def __init__(self, topn, extractor):
        self.topn = topn
        self.extractor = extractor
        assert self.extractor.length    # the extractor must have been already set up
        self.all_feats = set()
        self.meta_feat_dict = dict()

    def process_data(self, filename):
        with open(filename) as f:
            reader = csv.reader(f)  
            #extractor = setup_extractor(length)
            for i, row in enumerate(reader):
                reaction = row[0].split(" ")[0]
                reactants = reaction.split(">>")[0]
                mols = reactants.split(".")
                for smi in mols:
                    try:
                        L = label_each_atom(smi)    # some molecules are not valid and rdkit cannot handling sanitizing them
                    except:
                        continue
                    for atom_smi in L:
                        path_names = self.extractor.extract_path_names(atom_smi)
                        self.all_feats.update(path_names)
                if i%100==0:
                    print("%d reactions are processed..."%i)
            print("%d features are extracted from the entire dataset"%len(self.all_feats))
        meta_dict = {item: val for val, item in enumerate(self.all_feats)}
        self.meta_feat_dict = meta_dict
        
        return meta_dict
