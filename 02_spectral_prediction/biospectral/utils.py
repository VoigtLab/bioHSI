import pandas as pd
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from openbabel import pybel
import pickle
import argparse
import os
import re

"""
List of functions necessary for processing SMILES
and generating molecule 3D structures.
"""


def count_conjugated(mol: Chem.rdchem.Mol) -> int:
    """
    Returns the number of conjugated pi-bonds in a molecule
    """

    bonds = mol.GetBonds()
    return np.sum([bond.GetIsConjugated() for bond in bonds])


def np_and(s1, s2):
    """
    Returns the elementwise product of two boolean vectors
    """
    x = [s1, s2]
    return np.all(x, axis=0)


def is_metal(at: Chem.rdchem.Atom) -> bool:
    """
    Returns True if at is a metal atom, else False
    """

    n = at.GetAtomicNum()
    return (
        (n >= 3 and n <= 3)
        or (n >= 11 and n <= 12)
        or (n >= 19 and n <= 29)
        or (n >= 40 and n <= 47)
        or (n >= 72 and n <= 79)
    )


def contains_metal(mol: Chem.rdchem.Mol) -> bool:
    """
    Returns True if mol contains a metal atom, else False
    """
    return np.any([is_metal(at) for at in mol.GetAtoms()])


def embed_mol(mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """
    Generate a MMF94 optimized 3D embedding for mol
    """
    mol2 = Chem.AddHs(mol)

    try:
        embed = AllChem.EmbedMolecule(mol2)
        if embed == -1:
            raise RuntimeError
        AllChem.MMFFOptimizeMolecule(mol2)
    except (RuntimeError, ValueError):
        bel_mol = pybel.readstring("mol", Chem.MolToMolBlock(mol2, includeStereo=True))
        bel_mol.localopt(forcefield="uff", steps=10000)
        bel_mol.localopt(forcefield="mmff94", steps=10000)
        mol2 = Chem.MolFromMolBlock(bel_mol.write("mol"), removeHs=False)

    return mol2


def save_cartesian(mol: Chem.rdchem.Mol, file_name: str, overwrite: bool = False):
    """
    Save the 3D cartesian coordinates of a molecules
    """

    if not (file_name[-4:] == ".coord"):
        file_name = file_name + ".coord"

    xyz = []
    charge = str(Chem.GetFormalCharge(mol))
    multiplicity = str(
        np.sum([atm.GetNumRadicalElectrons() for atm in mol.GetAtoms()]) + 1
    )
    xyz.append(charge + " " + multiplicity + "\n")

    if (not overwrite) and not (os.path.exists(file_name) or os.path.exists('.'.join(file_name.split('.')[:-1]) + '.com')):
        for c in mol.GetConformers():
            pos = c.GetPositions()
            atms = mol.GetAtoms()
            for p, a in zip(pos, atms):
                xyz.append(
                    str(a.GetSymbol()) + "\t" + "\t".join([str(x) for x in p]) + "\n"
                )
    else:
        with open(file_name, "r") as f:
            existing = f.readlines()
        if re.match("^[A-Z]", existing[0]) != None:
            xyz = xyz + existing
        else:
            xyz = existing

    with open(file_name, "w") as f:
        f.writelines(xyz)
        f.close()


def standardize_smiles(smiles: str) -> str:
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    except:
        return smiles


def get_name(smiles: str, dct: dict) -> str:
    try:
        return dct[smiles]
    except:
        return smiles


def clean_name(string: str) -> str:
    new = re.sub("[\s\.,]", "_", string)
    new = re.sub("[\[\]\(\)']", "", new)
    new = re.sub('&rarr;', '->', new)
    new = re.sub('<[a-zA-z]*>|<\/[a-zA-z]*>|;|&|^[Aa]n |^[0-9]* ','', new)
    new = re.sub('\+', 'plus', new) 
    new = re.sub('^-', 'minus', new)
    new = re.sub(',', '-', new)
    return new


def get_mols(smiles):
    mols = [Chem.MolFromSmiles(smiles) for smile in smiles]
