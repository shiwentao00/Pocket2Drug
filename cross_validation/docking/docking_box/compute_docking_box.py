"""
Compute the docking boxes of the 14 label ligands.
"""
import os
from os import listdir
from os.path import isfile, join
import argparse
import subprocess
import yaml
from tqdm import tqdm


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

    # input and output directories
    args = get_args()
    in_dir = args.in_dir
    out_dir = args.out_dir

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print("Input directory: ", in_dir)
    print("Output directory: ", out_dir)

    # list of pocket folders, where each folder contains multiple pdbqt files
    pockets = [pocket for pocket in os.listdir(in_dir) if os.path.isdir(os.path.join(in_dir, pocket))]

    for pocket in tqdm(pockets):
        docking_boxes = {}

        pocket_dir = os.path.join(in_dir, pocket) 

        # lsit of sampled molecules
        mols = [f for f in listdir(pocket_dir) if isfile(join(pocket_dir, f))]

        # compute size of docking box of each molecule
        for mol in mols:
            mol_path = os.path.join(pocket_dir, mol)
            p = subprocess.run(eboxsize + ' {}'.format(mol_path),
                               shell=True,
                               stdout=subprocess.PIPE,
                               text=True)

            # when there is no error
            if p.returncode == 0:
                result = p.stdout
                # add to result dict
                docking_boxes[mol.split()[0]] = float(result.strip())
            else:
                print('Something went wrong, ligand: {}'.format(mol_path))

        # save docking box sizes of this pocket
        with open(out_dir + pocket + "_docking_box.yaml", "w") as f:
            yaml.dump(docking_boxes, f)