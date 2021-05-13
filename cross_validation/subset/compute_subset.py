"""
Using Rdkit to select a diverse subset of sampled molecules of 
a pocket in the validation set. 
"""
import argparse
import os
from os import listdir
from os.path import isfile, join
import yaml
from rdkit import Chem
# from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

# suppress rdkit error
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-mol_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample/",
                        help="directory of yaml files where each yaml is a list of sampled SMILES")

    parser.add_argument("-out_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered/",
                        help="directory of yaml files after being selected by clustering algorithm")

    return parser.parse_args()


def subset(smiles_list, picker, subset_size):
    """Select a diverse subset from a list of SMILES
    
    Arguments:
        smiles_list - a list of SMILES strings
        picker - a picker object of rdkit.SimDivFilters.rdSimDivPickers
    
    """
    mol_list = []
    valid_smiles = []
    for smiles in smiles_list:
        if len(smiles) > 0:
            mol = Chem.MolFromSmiles(smiles) 
            if mol:
                valid_smiles.append(smiles)
                mol_list.append(mol)
    
    # fingerprints = [GetMorganFingerprint(x, 3) for x in mol_list]
    fingerprints = [MACCSkeys.GenMACCSKeys(x) for x in mol_list]
    num_fingerprints = len(fingerprints)
    assert(num_fingerprints == len(mol_list))

    def distij(i, j, fps=fingerprints):
        """The distance between fingerprint i and j"""
        # The default similarity metric used by 
        # rdkit.DataStructs.FingerprintSimilarity() 
        # is the Tanimoto similarity.
        return 1-DataStructs.FingerprintSimilarity(fps[i], fps[j])

    pick_indices = picker.LazyPick(
        distij, num_fingerprints, subset_size, seed=23)
    
    pick_smiles = [valid_smiles[i] for i in pick_indices]

    return pick_smiles


if __name__ == "__main__":
    args = get_args()
    mol_dir = args.mol_dir
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print('input dir:', mol_dir)
    print('output dir', out_dir)

    # for each pocket file, read the SMILES of the sampled molecules
    pocket_files = [f for f in listdir(
        mol_dir) if isfile(join(mol_dir, f))]

    for pocket_file in pocket_files:
        pocket_name = pocket_file.split('_')[0]
        sampled_mol_path = join(mol_dir, pocket_file)
        with open(sampled_mol_path, 'r') as f:
            sampled_smiles = yaml.full_load(f)
        sampled_smiles = list(sampled_smiles.keys())
        picker = MaxMinPicker()
        subset_smiles = subset(sampled_smiles, picker, 100)

        subset_path = join(out_dir, pocket_name + '_subset.yaml')
        with open(subset_path, 'w') as f:
            yaml.dump(subset_smiles, f)

