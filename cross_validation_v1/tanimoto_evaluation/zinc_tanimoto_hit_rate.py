"""
Compute the "hit rate" of the genrated molecules according to their Tanimoto Coefficients 
with the label ligands of the pockets
"""
import argparse
import yaml
from os import listdir
from os.path import isfile, join
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import FingerprintSimilarity
from multiprocessing import Pool
import random
import pickle
from collections import defaultdict

# suppress rdkit error
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

class TSEvaluation:
    def __init__(self, mol_dir, label_smiles_path, control_path, seq_identity_path=None, num_workers=20):
        self.mol_dir = mol_dir
        self.num_workers = num_workers
        # tanimoto thresholds for hit
        self.ths = [0.7, 0.8, 0.9, 0.95, 1.0]
        self.cnt = {th:0 for th in self.ths}

        # list of pockects 
        self.pocket_files = [f for f in listdir(mol_dir) if isfile(join(mol_dir, f))]
        print(f"total number of pockets to evaluate: {len(self.pocket_files)}")

        if seq_identity_path:
            # load sequence identity
            with open(seq_identity_path, "r") as f:
                lines = f.readlines()
            seq_identity = defaultdict(lambda: float("inf"))
            for line in lines:
                line = line.rstrip()
                line = line.split()
                seq_identity[line[1]] = float(line[0])

            # filter out pockets that have similar pockets in train set
            temp = []
            for pocket_file in self.pocket_files:
                pocket = pocket_file.split('_')[0]
                protein = pocket[0:-2]
                seq_identity_score = seq_identity[protein]
                if seq_identity_score >= 0.5: continue
                temp.append(pocket_file)
            self.pocket_files = temp

        # true label SMILES
        with open(label_smiles_path, 'r') as f:
            self.label_smiles_dict = yaml.full_load(f)

        # read zinc dataset
        with open(control_path, 'rb') as f:
            bins = pickle.load(f)
        self.control_list = []
        for bin in bins:
            self.control_list.extend(bin)

    def compute_hit_rates(self):
        with Pool(self.num_workers) as p:
            # best tanimoto similarities for all pockets
            best_tcs = p.map(self.compute_single_pocket_ts, self.pocket_files)
        
        # compute hit rates
        error_pockets = 0
        for tc in best_tcs:
            if tc is None:
                error_pockets += 1
            else:
                for th in self.ths:
                    if tc >= th:
                        self.cnt[th] += 1
        print(f"Number of pockets with errors {error_pockets}")
        for th in self.ths:
            hit_rate = self.cnt[th] / (len(self.pocket_files) - error_pockets)
            print(f"Threshold: {th}, hit rate: {hit_rate}")


    def compute_single_pocket_ts(self, pocket_file):
        # generated molecules 
        with open(join(self.mol_dir, pocket_file), 'r') as f:
            mols = yaml.full_load(f)

        # sample same number of molecules from the vanilla rnn output
        num_mols = len(mols)
        mols = random.sample(self.control_list, num_mols)

        pocket = pocket_file.split('_')[0]
        try:
            label_mol = self.label_smiles_dict[pocket]
            label_mol = Chem.MolFromSmiles(label_mol)
            label_maccs = MACCSkeys.GenMACCSKeys(label_mol)
        except:
            print(f"something went wrong, pocket: {pocket}")
            return None

        # compute all the tanimoto similarities
        pocket_ts = []
        for mol in mols:
            try:
                ts = self.compute_similarity(label_maccs, mol)
                pocket_ts.append(ts)
            except:
                continue
        
        # get the best tanimoto similarity
        if not pocket_ts:
            print(f"no valid tc, pocket: {pocket}")
            return None
        print(pocket, max(pocket_ts))
        return max(pocket_ts)

    @staticmethod
    def compute_similarity(target_maccs, mol):
        """Compute the similarity between two smiles."""
        mol = Chem.MolFromSmiles(mol)
        mol_maccs = MACCSkeys.GenMACCSKeys(mol)
        return FingerprintSimilarity(target_maccs, mol_maccs)

def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument(
        "-mol_dir",
        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_81920/",
        required=False,
        help="directory of molecules to evaluate"
    )

    parser.add_argument(
        "-label_smiles_path",
        required=False,
        default="../../data/pocket-smiles.yaml",
        help="the label smiles of pockets"
    )

    parser.add_argument(
        "-control_path",
        required=False,
        default='../../../../zinc-sampler-project/zinc-dataset-bins.pickle',
        help="input molecules for molecular weight distribution"
    )

    parser.add_argument(
        "-seq_identity_path",
        required=False,
        default="../case_study/max_seq_identities.txt",
        help="the label smiles of pockets"
    )

    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    mol_dir = args.mol_dir
    label_smiles_path = args.label_smiles_path
    control_path = args.control_path
    seq_identity_path = args.seq_identity_path
    evaluation = TSEvaluation(mol_dir, label_smiles_path, control_path, seq_identity_path, 20)
    evaluation.compute_hit_rates()