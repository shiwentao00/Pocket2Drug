"""
Docking script for ligands sampled by the P2D model.
"""
import argparse
import pickle
import yaml
import os
from os import listdir
from os.path import isfile, join
from dock import smina_dock


def get_args():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-fold',
                        type=int,
                        default=0,
                        help='which cross-validation fold to perform docking')

    parser.add_argument('-start',
                        type=int,
                        default=0,
                        help='start index of pocket')

    parser.add_argument('-end',
                        type=int,
                        default=1,
                        help='end index of pocket')

    parser.add_argument("-ligand_dir",
                        required=False,
                        default="../../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered_pdbqt/",
                        help="pdbqt files of ligands for docking")

    parser.add_argument("-dock_box_dir",
                        required=False,
                        default="../../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered_pdbqt_dock_box/",
                        help="pre-computed docking box sizes for docking")

    parser.add_argument("-docking_center_path",
                        required=False,
                        default="../docking_center/pocket_center.pickle",
                        help="pre-computed docking ceters for docking")

    parser.add_argument("-out_dir",
                        required=False,
                        default="../../../../p2d_results_selfie/cv_results/cross_val_fold_0/docking_results/",
                        help="output directory for docking results")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    fold = args.fold

    # list of pockets of determined order
    fold_dir = "../../../cross_validation/folds/pockets_fold{}.yaml".format(
        fold)
    with open(fold_dir, 'r') as f:
        pockets = yaml.full_load(f)
    pockets = list(pockets.keys())
    pockets.sort()

    # smina installed via conda
    smina_path = 'smina'

    # directory of proteins for docking
    protein_dir = '../../../../smina/pdbqt_googlenet/'

    # root directory of ligands
    ligand_dir = args.ligand_dir

    # directory to store output files
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # pocket cneters
    docking_center_path = args.docking_center_path
    with open(docking_center_path, 'rb') as f:
        pocket_centers = pickle.load(f)

    # directory of docking box files
    dock_box_dir = args.dock_box_dir

    num_errors = 0
    start, end = int(args.start), int(args.end)
    for pocket in pockets[start:end + 1]:
        if pocket in pocket_centers:
            # path of protein
            protein = protein_dir + pocket[0:-2] + '.pdbqt'

            # center of docking search space
            pocket_center = pocket_centers[pocket]

            # docking scores this pocket
            scores = {}

            # list of ligands to dock
            pocket_ligand_dir = ligand_dir + pocket + '/'
            ligands = [f for f in listdir(
                pocket_ligand_dir) if isfile(join(pocket_ligand_dir, f))]

            # dictionary of docking boxes
            docking_box_path = dock_box_dir + \
                "{}_docking_box.yaml".format(pocket)
            with open(docking_box_path, "r") as f:
                docking_boxes = yaml.full_load(f)

            # for each ligand candidate of this pocket
            for mol in ligands[0:10]:
                # current ligand path
                mol_path = pocket_ligand_dir + mol

                # docking box (cubic) size
                docking_box = docking_boxes[mol]

                score, _ = smina_dock(
                    'smina',
                    protein,
                    mol_path,
                    pocket_center,
                    docking_box
                )

                if score is None:
                    num_errors += 1
                    continue

                scores[mol] = score

            # save docking results
            docking_score_path = out_dir + \
                "{}_docking_score.yaml".format(pocket)
            with open(docking_score_path, "w") as f:
                yaml.dump(scores, f)

    # number of erroneous dockings
    print('total number of erroneous during docking: ', num_errors)
