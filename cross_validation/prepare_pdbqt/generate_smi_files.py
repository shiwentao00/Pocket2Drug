"""
Generate smi files from the .yaml files containing SMILES.
"""
import argparse
import os
from os import listdir
from os.path import isfile, join
import yaml
from tqdm import tqdm

def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-in_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered/",
                        help="directory of yaml files where each yaml is a list of sampled SMILES")

    parser.add_argument("-out_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered_smi/",
                        help="directory of yaml files after being selected by clustering algorithm")

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    in_dir = args.in_dir
    out_dir = args.out_dir
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
            