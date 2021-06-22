"""
Compute the docking boxes of the 14 label ligands.
"""
from os import listdir
from os.path import isfile, join
import argparse

import os
import subprocess
import yaml
from tqdm import tqdm


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-fold",
                        required=False,
                        default=0,
                        help="which fold used for validation")

    parser.add_argument("-work_dir",
                        required=False,
                        default="../../../../p2d_results_selfie/cv_results/",
                        help="which fold used for validation")

    return parser.parse_args()


if __name__ == "__main__":
    # the tool used to compute the docking box size
    eboxsize = './eboxsize.pl'

    # input and output directories
    args = get_args()
    fold = args.fold
    work_dir = args.work_dir
    in_dir = work_dir + \
        "cross_val_fold_{}/val_pockets_sample_clustered_pdbqt/".format(fold)
    out_dir = work_dir + \
        "cross_val_fold_{}/val_pockets_sample_clustered_pdbqt_dock_box/".format(
            fold)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print("Input directory: ", in_dir)
    print("Output directory: ", out_dir)

    pocket_fold_dir = yaml.full_load(
        "../../folds/pockets_fold{}.yaml".format(fold))
    with open(pocket_fold_dir, "r") as f:
        pockets = yaml.full_load(f)
    pockets = pockets.keys()

    for pocket in tqdm(pockets):
        # for each pocket
        pocket_dir = in_dir + pocket + "/"
        docking_boxes = {}

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
