"""
Generate smi files from the .yaml files containing SMILES.
"""
import os
from os import listdir
from os.path import isfile, join
import yaml
from tqdm import tqdm


if __name__ == "__main__":
    in_dir = "../../../p2d_results/cv_results/cross_val_fold_0/val_pockets_sample_clustered/"
    out_dir = "../../../p2d_results/cv_results/cross_val_fold_0/val_pockets_sample_clustered_smi/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    pocket_subset_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]

    for pocket_subset_file in tqdm(pocket_subset_files):
        pocket_name = pocket_subset_file.split('_')[0]
        pocket_out_path = join(out_dir, pocket_name + ".smi")

        with open(join(in_dir, pocket_subset_file), 'r') as f:
            smiles_list = yaml.full_load(f)

        with open(pocket_out_path, 'w') as f:
            for smiles in smiles_list:
                f.write(smiles + '\n')
            