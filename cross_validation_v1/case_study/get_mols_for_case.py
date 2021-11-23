"""Get molecules for the input case pocket"""
import argparse
from os.path import join
import yaml
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import FingerprintSimilarity
from tqdm import tqdm
from collections import defaultdict


def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument(
        "-fold",
        required=True,
        help="which fold is used for validation"
    )

    return parser.parse_args()

def compute_similarity(target_maccs, mol):
    """Compute the similarity between two smiles."""
    mol = Chem.MolFromSmiles(mol)
    mol_maccs = MACCSkeys.GenMACCSKeys(mol)
    return FingerprintSimilarity(target_maccs, mol_maccs)

if __name__ == "__main__":
    # load sequence identity
    with open("./max_seq_identities.txt", "r") as f:
        lines = f.readlines()
    seq_identity = defaultdict(lambda: float("inf"))
    for line in lines:
        line = line.rstrip()
        line = line.split()
        seq_identity[line[1]] = float(line[0])

    # load validation pockets
    args = get_args()
    fold = args.fold
    with open(f"../../cross_validation/folds/pockets_fold{fold}.yaml", "r") as f:
        val_pockets = yaml.full_load(f)

    # directory of all sampled molecules
    mol_dir = "../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_81920/"
    out_dir = "./tanimoto_cases/"
    good = 0
    bad = 0
    for pocket in tqdm(val_pockets):
        # max sequence identity of this pocket to training set
        protein = pocket[0:-2]
        seq_identity_score = seq_identity[protein]
        if seq_identity_score > 0.5: continue

        # label molecule
        label_mol = val_pockets[pocket]
        try:
            label_mol = Chem.MolFromSmiles(label_mol)
            label_maccs = MACCSkeys.GenMACCSKeys(label_mol)
        except:
            print("fsomething went wrong with the label molecule, pocket: {pocket}")
            continue
        
        # load input mols
        mol_path = join(mol_dir, f"{pocket}_sampled_temp1.0.yaml")
        with open(mol_path, "r") as f:
            mols = yaml.full_load(f)

        th_8_9_list = []
        th_9_10_list = []
        th_10_list = []
        # for each molecule
        for mol in mols:
            # compute tanimoto similarity
            ts = compute_similarity(label_maccs, mol)

            if 0.8 <= ts < 0.9:
                th_8_9_list.append(mol)
            elif 0.9 <= ts < 1.0:
                th_9_10_list.append(mol)
            elif ts == 1.0:
                th_10_list.append(mol)

        if th_10_list:
            good += 1
            with open(join(out_dir, f"{pocket}_mols.txt"), "w") as f:
                f.write("0.8 <= ts < 0.9\n")
                for mol in th_8_9_list:
                    f.write(f"{mol}\n")
                f.write("0.9 <= ts < 1.0\n")
                for mol in th_9_10_list:
                    f.write(f"{mol}\n")
                f.write("ts = 1.0\n")
                for mol in th_10_list:
                    f.write(f"{mol}\n")
        else:
            bad += 1

    print(f"good ratio: {float(good/(good + bad))}")