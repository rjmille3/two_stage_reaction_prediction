import numpy as np
import rdkit
from rdkit import Chem

from reaction_prediction.atom.utils import *

import random
import json

class FeatureExtraction():
    
    def __init__(self, allid_file, extractor):
        with open(allid_file) as f:
            self.allid = json.load(f)
        self.extractor = extractor
        self.num_non_reactive_atom_samples = 8

    def atom_feat_vec(self, smi):
        """
        smi is a atom labeled SMILES string of a molecule
        allid is json with all the path features and 
        their corresponding ids
        """
        path_fv = np.array([0.0 for i in range(len(self.allid))])
        
        fd = dict()
       
        mol = mol_with_hydrogens(smi)
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum()==1:
                a = atom
                break

        paths = self.extractor.extract_path_names(smi)
        path_dict = path_list_to_count_dict(paths)
    
        physio_chem_feats = list(atom_features(a))
        
        for item in path_dict:
            if item in self.allid:
                path_fv[self.allid[item]] = path_dict[item]
                   

        fv = np.concatenate([physio_chem_feats, path_fv])
        
        return fv, path_dict

    
    def reaction_to_feat_vecs_source(self, row):
        """
        row is supposed to be a list as follows:
        [reaction, arrows, source_atom, sink_atom]
        """
        reaction, arrows, source_atom, sink_atom = row
        
        fvs = []

        reactants = reaction.split(">>")[0]
        for mol_smi in reactants.split("."):
            try:
                atom_smis = label_each_atom(mol_smi)
            except:
                continue
            
            for smi in atom_smis:
                mol = mol_with_hydrogens(smi)
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum()==1:
                        a = atom
                        break
                
                if atom.GetNumRadicalElectrons()>0:
                    continue
                atom_fv, _ = self.atom_feat_vec(smi)
                fvs.append(atom_fv) 
        
        try:
            source_fv, _ = self.atom_feat_vec(source_atom) 
        except:
            return [], []

        final_fvs = []
        for item in fvs:
            if np.array_equal(item, source_fv):
                continue
            else:
                final_fvs.append(item)
        targs = [0.0 for i in range(len(final_fvs))]

        final_fvs.append(source_fv)
        for i in range(1):
            targs.append(1.0)

        return final_fvs, targs
    
    def reaction_to_feat_vecs_sink(self, row):
        """
        row is supposed to be a list as follows:
        [reaction, arrows, source_atom, sink_atom]
        """
        reaction, arrows, source_atom, sink_atom = row
        
        fvs = []

        reactants = reaction.split(">>")[0]
        for mol_smi in reactants.split("."):
            try:
                atom_smis = label_each_atom(mol_smi)
            except:
                continue
            
            for smi in atom_smis:
                mol = mol_with_hydrogens(smi)
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum()==1:
                        a = atom
                        break
                
                if atom.GetNumRadicalElectrons()>0:
                    continue
                atom_fv, _ = self.atom_feat_vec(smi)
                fvs.append(atom_fv)
        
        try:
            sink_fv, _ = self.atom_feat_vec(sink_atom)
        except:
            return [], []

        final_fvs = []
        for item in fvs:
            if np.array_equal(item, sink_fv):
                continue
            else:
                final_fvs.append(item)
        targs = [0.0 for i in range(len(final_fvs))]

        final_fvs.append(sink_fv)
        for i in range(1):
            targs.append(1.0)

        return final_fvs, targs

    def yield_samples(self, grp, num_samples):
        if len(grp) <= num_samples:
            return grp
        else:
            return random.sample(grp, num_samples)
    


