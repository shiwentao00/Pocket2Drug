"""
Analyse the Synthetic Accessibility scores of:
    1. Chembl28 dataset
    2. Molecules sampled from fold 0 pockets by Pocket2Drug, (char-level tokenization)
    3. Molecules sampled from fold 0 pockets by Pocket2Drug, (selfies tokenization)
"""
import os
from os import listdir
from os.path import isfile, join
from rdkit.Chem.QED import qed
from rdkit import Chem
from sascorer import calculateScore
import yaml
from tqdm import tqdm


def read_smiles_file(path):
    """Read .smi files. Need to exclude first line which is not SMILES"""
    with open(path, "r") as f:
         smiles = [line.strip("\n") for line in f.readlines()]
    
    return smiles


def load_sampled_smiles(folder_dir):
    """load a folder of yaml files, each file is a list of SMILES."""
    pocket_files = [f for f in listdir(
        folder_dir) if isfile(join(folder_dir, f))]
    
    smiles_list = []

    for pocket_file in pocket_files:
        sampled_mol_path = join(folder_dir, pocket_file)
        with open(sampled_mol_path, 'r') as f:
            sampled_smiles = yaml.full_load(f)
        smiles_list.extend(sampled_smiles)
    
    return smiles_list


def compute_scores(mol):
    """
    Returns QED score and SA score
    """
    # using default weights
    # MW=0.66, 
    # ALOGP=0.46, 
    # HBA=0.05, 
    # HBD=0.61, 
    # PSA=0.06, 
    # ROTB=0.65, 
    # AROM=0.48, 
    # ALERTS=0.95
    try:
        qed_score = qed(mol)
    except:
        print("Something went wrong when computing QED.")
        qed_score = None

    try:
        sa_score = calculateScore(mol)
    except:
        print("Something went wrong when computing SA.")
        sa_score = None

    return qed_score, sa_score


def compute_all_scores(smiles_list, qed_path, sa_path):
    qed_scores = []
    sa_scores = []
    num_invalid_qed = 0
    num_invalid_sa = 0
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            qed_score, sa_score = compute_scores(mol)
            if qed_score is None:
                num_invalid_qed += 1
            if sa_score is None:
                num_invalid_sa += 1
            if qed_score is not None and sa_score is not None:
                qed_scores.append(qed_score)
                sa_scores.append(sa_score)

    # compute valid ratio
    num_smiles = len(smiles_list)
    print("total number of smiles: {}".format(num_smiles))
    print("valid smiles ratio: {}".format(len(qed_scores)/num_smiles))
    print("number invalid QED: {}".format(num_invalid_qed))
    print("number of invalid SA: {}".format(num_invalid_sa))
    print("-------------------------------------------------")

    # save the scores
    with open(qed_path, "w") as f:
        yaml.dump(qed_scores, f)
    with open(sa_path, "w") as f:
        yaml.dump(sa_scores, f)


if __name__ == "__main__":
    chembl28_dataset_path = "/home/wentao/Desktop/local-workspace/molecule-generator-project/Molecule-RNN/dataset/chembl28-cleaned.smi"
    p2d_char_dataset_path = "/home/wentao/Desktop/local-workspace/pocket2drug-project/p2d_results/cv_results/cross_val_fold_0/val_pockets_sample_clustered/"
    p2d_selfies_dataset_path = "/home/wentao/Desktop/local-workspace/pocket2drug-project/p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered/"

    # chembl28_dataset = read_smiles_file(chembl28_dataset_path)
    # compute_all_scores(
        # chembl28_dataset, 
        # "./chembl28_qed_scores.yaml", 
        # "./chembl28_sa_scores.yaml"
    # )

    p2d_char_dataset = load_sampled_smiles(p2d_char_dataset_path)
    compute_all_scores(
        p2d_char_dataset,
        "./p2d_char_qed_scores.yaml",
        "./p2d_char_sa_scores.yaml"
    )

    p2d_selfies_dataset = load_sampled_smiles(p2d_selfies_dataset_path)
    compute_all_scores(
        p2d_selfies_dataset,
        "./p2d_selfies_qed_scores.yaml",
        "./p2d_selfies_sa_scores.yaml"
    )
