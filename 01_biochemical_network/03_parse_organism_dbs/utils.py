from rdkit import Chem
import re

def neutralize_atoms(smiles):
    #pulled from http://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return Chem.MolToSmiles(mol)
    
    
def standardize_smiles(smiles):
    if ' ' in smiles or 'R' in smiles:
        return smiles
    try:
        smiles = neutralize_atoms(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
        return smiles
    except:
        if Chem.MolFromSmiles(smiles)==None:
            m = Chem.MolFromSmiles(smiles, sanitize=False)
            fix_c = AllChem.ReactionFromSmarts('[#6-:1]>>[C;+0:1]')
            if m:
                f = fix_c.RunReactants([m])
                if len (f)>0:
                    f = f[0][0]
                    print ('fixed smiles :', smiles,Chem.MolToSmiles(f))
                    return neutralize_atoms(Chem.MolToSmiles(f))
                else:
                    return smiles
            else:
                return smiles
        else:
            return smiles