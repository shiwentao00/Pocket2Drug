"""
Run docking with the generated pdbqt files
"""
import subprocess
import argparse
import os
from os import listdir
from os.path import isfile, join
from tqdm import tqdm
import yaml


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


def parse_search_config(path):
    """Parse the earch configuation files.
       Returns (x, y, z centers, x, y, z dimensions)
    """
    f = open(path, 'r')
    info = f.readlines()
    f.close()
    centers = info[1]
    dimensions = info[2]
    centers = centers.split()
    dimensions = dimensions.split()
    out = []
    out.extend(centers[-3:])
    out.extend(dimensions[-3:])
    out = [float(x) for x in out]
    return out


def dock(vina_path, protein_path, ligand_path, config, out_path, exhaustiveness=8, num_modes=3, energy_range=99):
    """
    Call Autodock VINA program to perform docking.
    Arguments:
        vina_path - path of autodock vina's executable
        protein_path - pdbqt file of protein
        ligand_path - pdbqt file of ligand
        config - tuple containing 6 numbers (x, y, z centers, x, y, z dimensions).
    Return:
        docking_score: free energy left, the lower the more robust the docking.
    """
    p = subprocess.run(vina_path + ' --receptor {} --ligand {}'.format(protein_path, ligand_path) +
                                   ' --center_x {} --center_y {} --center_z {}'.format(
                                       config[0], config[1], config[2]) +
                                   ' --size_x {} --size_y {} --size_z {}'.format(
                                       config[3], config[4], config[5]) +
                                   ' --out {}'.format(out_path) +
                                   ' --exhaustiveness {}'.format(exhaustiveness) +
                                   ' --num_modes {}'.format(num_modes) +
                                   ' --energy_range {}'.format(energy_range),
                       shell=True,
                       stdout=subprocess.PIPE,
                       text=True)

    # check=True,
    # cwd='/home/wentao/Desktop/local-workspace/siamese-monet-project/glosa/glosa_v2.2/')  # working directory

    # when there is no error
    if p.returncode == 0:
        result = p.stdout
        return float(result.split()[-15])
    else:
        global error_cnt
        error_cnt += 1
        return 0

if __name__ == "__main__":
    # number of errors occured during docking
    global error_cnt
    error_cnt = 0

    args = get_args()
    result_dir = args.result_dir
    temperature = str(args.temperature)
    pockets = ['4jdwA00', '6bwhB00', '1t5cA01', '5yhzA00', '2hs0A00']

    # path of autodock vina's executable
    vina_path = '/home/wentao/Desktop/local-workspace/siamese-monet-project/vina/software/autodock_vina_1_1_2_linux_x86/bin/vina'

    # search space file directory
    cog_dir = '/home/wentao/Desktop/local-workspace/siamese-monet-project/vina/data/COG/'

    # protein file directory
    protein_dir = '/home/wentao/Desktop/local-workspace/siamese-monet-project/vina/data/proteins/'

    # for each pocket
    for pocket in pockets:
        print('docking on pocket {}...'.format(pocket))
        # output directory to store results
        out_dir = result_dir + 'docking_result_' + pocket + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        # path to store output file
        out_path = out_dir + pocket + '.out'

        # protein file
        protein = protein_dir + pocket[0:-2] + '.pdbqt'

        # search space configuration
        search_config = parse_search_config(cog_dir + pocket + '.out')

        # ligand files
        ligand_dir = result_dir + 'pdbqt_sampled_' + pocket + temperature + '/'
        ligand_files = [join(ligand_dir, f) for f in listdir(
            ligand_dir) if isfile(join(ligand_dir, f))]
        
        docking_scores = []
        # dock each ligand and save the docking score 
        for ligand_file in tqdm(ligand_files):
            score = dock(vina_path, protein, ligand_file,
                         search_config, out_path)
            if score != 0:
                docking_scores.append(score)

        # docking scores path
        docking_scores_path = out_dir + pocket + 'docking_socres.yaml'
        with open(docking_scores_path, 'w') as f:
            yaml.dump(docking_scores, f)

        
