"""
Select SMILES that have higher Tanimoto Smilarities,
then convert them to pdbqt files for docking
"""
import argparse
import yaml
import os
from tqdm import tqdm
import rdkit.Chem as Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import AllChem
from subprocess import call, STDOUT
import subprocess


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-result_dir",
                        required=True,
                        help="directory of files of sampled SMILES")

    parser.add_argument("-temperature",
                        required=False,
                        default=4.0,
                        help="the temperature paramter used during sampling")

    parser.add_argument("-th",
                        required=False,
                        default=0.9,
                        help="the threshold of tanimoto similarity to select ligands for docking")

    return parser.parse_args()


def compute_similarity(target_maccs, mol):
    """compute the similarity between two smiles"""
    mol_maccs = MACCSkeys.GenMACCSKeys(mol)
    return FingerprintSimilarity(target_maccs, mol_maccs)


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
    FNULL = open(os.devnull, 'w')
    args = get_args()
    result_dir = args.result_dir
    temperature = str(args.temperature)
    th = float(args.th)
    assert(0 <= th <= 1.0)
    pockets = ['5yhzA00', '4jdwA00', '6bwhB00', '1t5cA01', '2hs0A00']

    # dict of target SMILES
    smiles_dir = "../data/pocket-smiles.yaml"
    with open(smiles_dir, 'r') as f:
        smiles_dict = yaml.full_load(f)

    for pocket in pockets:
        print('pocket: ', pocket)

        # target SMILES
        target_smiles = smiles_dict[pocket]
        target_mol = Chem.MolFromSmiles(target_smiles)
        target_maccs = MACCSkeys.GenMACCSKeys(target_mol)

        # sampled smiles file
        sampled_smiles_file = result_dir + pocket + '_sampled_temp' + temperature + '.yaml'

        # read all sampled smiles of the pocket
        with open(sampled_smiles_file, 'r') as f:
            sampled_smiles_dict = yaml.full_load(f)

        # make a directory for ligand files
        out_dir = result_dir + 'selective_docking_' + pocket + temperature + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # path to store docking output file
        out_path = out_dir + pocket + '.out'

        # path of autodock vina's executable
        vina_path = '/home/wentao/Desktop/local-workspace/siamese-monet-project/vina/software/autodock_vina_1_1_2_linux_x86/bin/vina'

        # search space file directory
        cog_dir = '/home/wentao/Desktop/local-workspace/siamese-monet-project/vina/data/COG/'

        # protein file directory
        protein_dir = '/home/wentao/Desktop/local-workspace/siamese-monet-project/vina/data/proteins/'

        # for each smile
        docking_scores = {}
        for smiles in tqdm(list(sampled_smiles_dict.keys())):
            # read the smiles with rdkit
            mol = Chem.MolFromSmiles(smiles)

            # save as molecule file
            if mol is not None:
                try:
                    # compute similarity
                    s = compute_similarity(target_maccs, mol)

                    if s >= th:
                        mol = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(mol)
                        # save the ligand as a temporary mol file
                        Chem.MolToMolFile(mol, './temp.mol')

                        # convert to pdbqt file
                        call(
                            ('obabel imol {0} -O {1}'.format('./temp.mol',
                                                             './temp.pdbqt')).split(),
                            stdout=FNULL,
                            stderr=STDOUT
                        )
                    
                        # protein file
                        protein = protein_dir + pocket[0:-2] + '.pdbqt'

                        # search space configuration
                        search_config = parse_search_config(cog_dir + pocket + '.out')

                        # run docking
                        score = dock(vina_path, protein, './temp.pdbqt',
                                     search_config, out_path)
                        if score != 0:
                            docking_scores[smiles] = score
                except:
                    print('Something went wrong for pocket {}.'.format(pocket))
                    print('SMILES: ', smiles)
                    print('-----------------------------------------')

        # save the docking scores of selected molecules
        with open(out_dir + pocket + '_docking_score_th{}.yaml'.format(th), 'w') as f:
            yaml.dump(docking_scores, f)
