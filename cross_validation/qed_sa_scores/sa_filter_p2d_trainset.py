"""
Filter out the SMILES above the threshold and visualize.
"""
from rdkit_contrib.sascorer import calculateScore
from rdkit import Chem
import yaml


def read_smiles_yaml(path):
    with open(path, "r") as f:
        smiles_dict = yaml.full_load(f)
    return smiles_dict


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
        if i >= 20:
            return


if __name__ == "__main__":
    p2d_dataset_path = "../../data/pocket-smiles.yaml"
    smiles_dict = read_smiles_yaml(p2d_dataset_path)
    pockets = list(smiles_dict.keys())

    sa_9_10 = []
    sa_8_9 = []
    sa_7_8 = []
    sa_6_7 = []
    sa_5_6 = []
    sa_4_5 = []
    sa_small = []

    for pocket in pockets:
        smiles = smiles_dict[pocket]
        mol = Chem.MolFromSmiles(smiles)

        if mol:
            sa_score = compute_sa_score(mol)
            if sa_score:
                if 9 < sa_score <= 10:
                    sa_9_10.append((pocket, smiles))
                elif 8 < sa_score <= 9:
                    sa_8_9.append((pocket, smiles))
                elif 7 < sa_score <= 8:
                    sa_7_8.append((pocket, smiles))
                elif 6 < sa_score <= 7:
                    sa_6_7.append((pocket, smiles))
                elif 5 < sa_score <= 6:
                    sa_5_6.append((pocket, smiles))
                elif 4 < sa_score <= 5:
                    sa_4_5.append((pocket, smiles))
                else:
                    sa_small.append((pocket, smiles))

        if len(sa_9_10) > 20 and \
           len(sa_8_9) > 20 and \
           len(sa_7_8) > 20 and \
           len(sa_6_7) > 20 and \
           len(sa_5_6) > 20 and \
           len(sa_4_5) > 20 and \
           len(sa_small) > 20:
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
