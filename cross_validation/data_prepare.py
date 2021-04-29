"""
This script prepare the data for a 10-fold cross validation. 
 1. The proteins from the DUD-E dataset is excluded.
 2. The pockets used in previous case studies are excluded.
 3. The pockets are shuffled and divided into 10 folds.
 4. Each fold is saved
"""
import yaml
import random

if __name__ == "__main__":
    # all pockets
    smiles_dir = "../data/pocket-smiles.yaml"
    with open(smiles_dir, 'r') as f:
        smiles_dict = yaml.full_load(f)
    pockets = list(smiles_dict.keys())

    # load DUD-E dataset and case study pockets
    with open('./dude/dude_proteins.yaml', 'r') as f:
        dude_proteins = yaml.full_load(f)
    dude_proteins = set(dude_proteins)
    case_study_pockets = set(['4jdwA00', '6bwhB00', '1t5cA01', '5yhzA00', '2hs0A00'])  
    to_exclude = dude_proteins.union(case_study_pockets)

    # exclude DUD-E dataset and case study pockets
    new_pockets = []
    for pocket in pockets:
        if pocket[0:4] not in to_exclude:
            new_pockets.append(pocket)
    # print(len(pockets))
    # print(len(new_pockets))

    # shuffle
    random.shuffle(new_pockets)

    # divide into 10 folds
    num_pockets = len(new_pockets)
    fold_len = int(num_pockets / 10)

    # for each fold
    for f in range(10):
        start = f * fold_len
        if f == 9:
            end = max([(f + 1) * fold_len, num_pockets])
        else:
            end = (f + 1) * fold_len

        # the dictionary to save
        # key is pocket name, value is the target SMILES
        pocket_smiles_dict = {}
        for pocket in new_pockets[start:end]:
            pocket_smiles_dict[pocket] = smiles_dict[pocket]
        
        # save the pocket name - SMILES pairs of this fold
        with open('./folds/pockets_fold{}.yaml'.format(f), 'w') as f:
            yaml.dump(pocket_smiles_dict, f)

