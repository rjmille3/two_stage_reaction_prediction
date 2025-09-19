import numpy as np
import rdkit

from rdkit import Chem

from reaction_prediction.atom.utils import *

from reaction_prediction.atom.modules.path_extractor import PathExtractor
from reaction_prediction.atom.modules.feat_sel import FeatSel

import sys
import os
import argparse
import json
import csv


def all_features_to_json(input_file, output_file, length):
    extractor = PathExtractor(length)
    feat_selector = FeatSel(1500, extractor)
    meta_dict = feat_selector.process_data(input_file)
    with open(output_file, 'w') as f:
        json.dump(meta_dict, f)
    f.close()

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
        help="Path to save the extracted features JSON",
    )
    parser.add_argument(
        "--length", "-l",
        type=int,
        required=True,
        help="Path length to use for feature extraction",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    all_features_to_json(args.input, args.output, args.length)