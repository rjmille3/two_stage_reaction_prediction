import rdkit
import numpy as np
import h5py
from rdkit import Chem
import sys
import csv
import argparse

from reaction_prediction.atom.utils import *

from reaction_prediction.atom.modules.feature_extraction import FeatureExtraction
from reaction_prediction.atom.modules.path_extractor import PathExtractor

def main(filename, out_hdf5, allid):
    extractor = PathExtractor(length=3)
    #extractor = PathExtractor(length=4)
    allid_file = allid
    fe = FeatureExtraction(allid_file, extractor)

    with open(filename, 'r') as f:
        all_fvs = []
        all_targs = [] 
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            fv, targ = fe.reaction_to_feat_vecs_sink(row)     
            all_fvs += fv
            all_targs += targ
            if i%100==0:
                print("%d number of reactions are processed..."%i, flush=True)
                print("the feature vector of %d atoms are created"%len(all_fvs), flush=True) 
        
        all_fvs, all_targs = shuffle_two_lists(all_fvs, all_targs)
        
        all_fvs = np.array(all_fvs)
        all_targs = np.array(all_targs)
        
        hf = h5py.File(out_hdf5, 'w', libver='latest')
        atom_features = hf.create_dataset('features', data=all_fvs,  dtype='float32')
        print("atom features are inserted into the hdf5 dataset")

        targets = hf.create_dataset('targets', data=all_targs,  dtype='float32')
        print("targets are inserted into the hdf5 dataset")
        
        hf.close()

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract features from reaction data and save to JSON."
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to the input file containing reaction data",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path to save the extracted features HDF5",
    )
    parser.add_argument(
        "--allid", "-a",
        required=True,
        help="Allid file containing graph topological features",
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.input, args.output, args.allid)
