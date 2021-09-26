"""
Compute the docking boxes of the 14 label ligands.
"""
import os
from os import listdir
from os.path import isfile, join
import argparse
import subprocess
import yaml


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-in_dir",
                        required=False,
                        default="../../../../p2d_results_selfie/cv_results/cross_val_fold_0/zinc_sampled_pdbqt/",
                        help="input directory of the pdbqt files")

    parser.add_argument("-out_dir",
                        required=False,
                        default="../../../../p2d_results_selfie/cv_results/cross_val_fold_0/zinc_sampled_pdbqt_dockbox/",
                        help="output directory of the docking boxes")

    return parser.parse_args()


if __name__ == "__main__":
    # the tool used to compute the docking box size
    eboxsize = './eboxsize.pl'

    args = get_args()
    in_dir = args.in_dir
    out_dir = args.out_dir

    # list of pocket folders, where each folder contains multiple pdbqt files
    subfolders = [os.path.join(in_dir, pocket) for pocket in os.listdir(in_dir) 
                    if os.path.isdir(os.path.join(in_dir,pocket))]

    for pocket_folder in subfolders:
        pdbqt_files = [f for f in listdir(pocket_folder) if isfile(join(pocket_folder, f))]
        if len(pdbqt_files) != 100:
            print(pocket_folder, len(pdbqt_files))


    
