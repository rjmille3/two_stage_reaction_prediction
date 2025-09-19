import sys
import csv
import traceback
import argparse

from openeye.oechem import *
from reaction_chem.modules.Reaction import *


def change_label_to_one(s):
    mol = OEGraphMol()
    OESmilesToMol(mol, s)
    for atom in mol.GetAtoms():
        if atom.GetMapIdx() == 10:
            atom.SetMapIdx(1)
    smiles = OEMolToSmiles(mol)
    return smiles

def main(input_filename, output_filename):
    f = open(input_filename, 'r')
    reader = csv.reader(f)
    
    g = open(output_filename, 'w')
    writer = csv.writer(g)
    
    for row in reader:
        #smirks, arrows = row[0].split(" ")
        try:
            smirks, arrows = row[0].split(" ")
            arrows = arrows.rstrip(',')

            rxn = Reaction(smirks, arrows)
        except Exception as e:
            print("row is unable to be processed", e)
            sys.exit()
            continue
        src = change_label_to_one(rxn.srcAtom)
        sink = change_label_to_one(rxn.sinkAtom)
        
        new_row = [smirks, arrows, src, sink]
        writer.writerow(new_row)

    f.close()
    g.close()


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
        help="Path to save the reformatted reactions",
    )
    return parser.parse_args()



if __name__ == "__main__":
    args = parse_args()
    main(args.input, args.output)


