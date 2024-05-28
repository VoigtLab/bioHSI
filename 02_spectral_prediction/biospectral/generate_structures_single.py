from rdkit import Chem
import pickle
import argparse
import os
import pandas as pd
import sys

from biospectral import utils

"""
Takes as input a list of molecules and outputs a 
3D structure for each molecule in the list.
"""


def get_args():
    options = argparse.ArgumentParser()
    options.add_argument("--smiles", dest="smiles", type=str)
    options.add_argument("--name", dest="name", type=str)
    options.add_argument("--out-dir", dest="out_dir", type=str)
    options.add_argument("--overwrite", dest="overwrite", action="store_true", default=False)
    args = options.parse_args()
    return args


def generate_structure(name: str, smile: str, out_dir: str, overwrite=False):
    try:
        file_name = os.path.join(out_dir, utils.clean_name(str(name)))
        if (not overwrite) and not (os.path.exists(file_name) or os.path.exists(file_name+'.coord') \
                                    or os.path.exists('.'.join(file_name.split('.')[:-1]) + '.com')):
            utils.save_cartesian(
                utils.embed_mol(Chem.MolFromSmiles(smile)),
                file_name,
                overwrite=overwrite
            )
            print(name, "done")
        else:
            print ("File already exists for {}. Not overwriting".format(name))
    except:
        print("Could not save coordinate file for", name)
        print(sys.exc_info()[0])

def main():
    args = get_args()
    smile = args.smiles
    name = args.name
    out_dir = args.out_dir
    overwrite = args.overwrite

    try:
        if not os.path.isdir(out_dir):
             os.mkdir(out_dir)
    except:
         pass
    
    generate_structure(name, smile, out_dir, overwrite)
