"""
Take the mol files and convert them into pdbqt files using openbabel.
"""
import argparse
import yaml
import os
from os import listdir
from os.path import isfile, join
from tqdm import tqdm
from subprocess import call, STDOUT


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-result_dir",
                        required=True,
                        help="directory of files of sampled SMILES")

    parser.add_argument("-temperature",
                        required=False,
                        default=4.0,
                        help="the temperature paramter used during sampling")

    return parser.parse_args()


if __name__ == "__main__":
    FNULL = open(os.devnull, 'w')
    args = get_args()
    result_dir = args.result_dir
    temperature = str(args.temperature)
    pockets = ['4jdwA00', '6bwhB00', '1t5cA01', '5yhzA00', '2hs0A00']
    for pocket in pockets:
        print('generating pdbqt files for pocket {}...'.format(pocket))
        mol_folder = result_dir + "mol_sampled_" + pocket + temperature

        # make a directory for pdbqt files
        out_dir = result_dir + 'pdbqt_sampled_' + pocket + temperature + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # get list of mol files
        mol_files = [join(mol_folder, f) for f in listdir(mol_folder) if isfile(join(mol_folder, f))]

        # use openbabel to convert them into pdbqt
        for i, mol_file in tqdm(enumerate(mol_files)):
            pdbqt_file = out_dir + str(i) +  '.pdbqt'
            try:
                call(
                    ('obabel imol {0} -O {1}'.format(mol_file, pdbqt_file)).split(),
                    stdout=FNULL, 
                    stderr=STDOUT
                )
            except:
                print('Something went wrong for pocket {}, mol file {}'.format(pocket, mol_file))
                print('---------------------------------------------------------')