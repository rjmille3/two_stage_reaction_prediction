import numpy as np
import rdkit

from rdkit import Chem
from rdkit.Chem import AllChem

from reaction_prediction.atom.utils import *
from reaction_prediction.atom.modules.path_extractor import PathExtractor as PE
from reaction_prediction.atom.modules.feature_extraction import FeatureExtraction as FE


class ReactionFeatureExtraction():
    def __init__(self, morgan_radi=2, morgan_bits=2048, atom_length=3, allid_file=''):
        self.atom_length = atom_length
        self.morgan_radi = morgan_radi
        self.morgan_bits = morgan_bits
        self.setup_atom_feat_extractor(allid_file, atom_length)

    def numpy_ecfp_fingerprint(self, smi):
        f_dict = {}
        mol = Chem.MolFromSmiles(smi)
        AllChem.GetMorganFingerprintAsBitVect(mol, radius=self.morgan_radi, nBits=self.morgan_bits, bitInfo=f_dict)
        fp = np.zeros(self.morgan_bits, dtype=np.float32)
        for item in f_dict.items():
            fp[item[0]] = float(len(item[1]))

        return fp

    def diff_morgan(self, react, prods):
        r = np.zeros(self.morgan_bits)
        p = np.zeros(self.morgan_bits)

        for molr in react:
            r += self.numpy_ecfp_fingerprint(molr)
        for molp in prods:
            p += self.numpy_ecfp_fingerprint(molp)

        return r-p

    def setup_atom_feat_extractor(self, allid_file, atom_length):
        path_extractor = PE(atom_length)
        feature_extractor = FE(allid_file, path_extractor)
        self.atom_feat_extractor = feature_extractor

    def extract_reactive_pair_features(self, source, sink):
        source_fv, _ = self.atom_feat_extractor.atom_feat_vec(source)
        sink_fv, _ = self.atom_feat_extractor.atom_feat_vec(sink)

        return np.concatenate((source_fv, sink_fv))

    def extract_rxn_rep(self, reaction, source, sink):
        reactants, prods = reaction.split(">>")
        react_mols = reactants.split(".")
        prod_mols = prods.split(".")

        morgan_part = self.diff_morgan(react_mols, prod_mols)
        atom_part = self.extract_reactive_pair_features(source, sink)

        return np.concatenate((morgan_part, atom_part))
        
