"""
Filter out the SMILES above the threshold and visualize.
"""
import random
from compute_scores import load_sampled_smiles, read_smiles_file
from sascorer import calculateScore
from rdkit import Chem


def compute_sa_score(mol):
    try:
        sa_score = calculateScore(mol)
    except:
        print("Something went wrong when computing SA.")
        sa_score = None

    return sa_score


def show_smiles_list(smiles_list):
    for i, smiles in enumerate(smiles_list):
        print(smiles)
        if i >= 6:
            return
    print("------------------------------------")


if __name__ == "__main__":
    # chembl28_dataset_path = "/home/wentao/Desktop/local-workspace/molecule-generator-project/Molecule-RNN/dataset/chembl28-cleaned.smi"
    # smiles_list = read_smiles_file(chembl28_dataset_path)

    # p2d_char_dataset_path = "/home/wentao/Desktop/local-workspace/pocket2drug-project/p2d_results/cv_results/cross_val_fold_0/val_pockets_sample_clustered/"
    # smiles_list = load_sampled_smiles(p2d_char_dataset_path)
    
    p2d_selfies_dataset_path = "/home/wentao/Desktop/local-workspace/pocket2drug-project/p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered/"
    smiles_list = load_sampled_smiles(p2d_selfies_dataset_path)
    

    random.shuffle(smiles_list)

    sa_9_10 = []
    sa_8_9 = []
    sa_7_8 = []
    sa_6_7 = []
    sa_5_6 = []
    sa_4_5 = []
    sa_small = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:
            sa_score = compute_sa_score(mol)
            if sa_score:
                if 9 < sa_score <= 10:
                    sa_9_10.append(smiles)
                elif 8 < sa_score <= 9:
                    sa_8_9.append(smiles)
                elif 7 < sa_score <= 8:
                    sa_7_8.append(smiles)
                elif 6 < sa_score <= 7:
                    sa_6_7.append(smiles)
                elif 5 < sa_score <= 6:
                    sa_5_6.append(smiles)
                elif 4 < sa_score <= 5:
                    sa_4_5.append(smiles)       
                else:
                    sa_small.append(smiles)

        if len(sa_9_10) > 6 and \
           len(sa_8_9) > 6 and \
           len(sa_7_8) > 6 and \
           len(sa_6_7) > 6 and \
           len(sa_5_6) > 6 and \
           len(sa_4_5) > 6 and \
           len(sa_small) > 6:
            break 

    print("sa score 9-10: ")
    show_smiles_list(sa_9_10)

    print("sa score 8-9: ")
    show_smiles_list(sa_8_9)

    print("sa score 7-8: ")
    show_smiles_list(sa_7_8)

    print("sa score 6-7: ")
    show_smiles_list(sa_6_7)

    print("sa score 5-6: ")
    show_smiles_list(sa_5_6)

    print("sa score 4-5: ")
    show_smiles_list(sa_4_5)

    print("sa score below 4: ")
    show_smiles_list(sa_small)

    