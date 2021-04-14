"""
Take the SMILES strings as input, convert and save them as 
molecule files. 
"""
import argparse
import yaml
import os
from tqdm import tqdm
import rdkit.Chem as Chem
from rdkit.Chem import AllChem


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
    args = get_args()
    result_dir = args.result_dir
    temperature = str(args.temperature)
    pockets = ['4jdwA00', '6bwhB00', '1t5cA01', '5yhzA00', '2hs0A00']
    for pocket in pockets:
        print('pocket: ', pocket)
        # smiles file
        smiles_file = result_dir + pocket + '_sampled_temp' + temperature + '.yaml'

        # read all sampled smiles of the pocket
        with open(smiles_file, 'r') as f:
            smiles_dict = yaml.full_load(f)

        # make a directory for ligand files
        out_dir = result_dir + 'mol_sampled_' + pocket + temperature + '/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # for each smile
        for i, smiles in tqdm(enumerate(list(smiles_dict.keys()))):
            # read the smiles with rdkit
            mol = Chem.MolFromSmiles(smiles)

            # save as molecule file
            if mol is not None:
                try:
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                    mol_path = out_dir + str(i) + '.mol'
                    Chem.MolToMolFile(mol, mol_path)
                except:
                    print(
                        'something went wrong for pocket {}, {}th SMILES.'.format(
                            pocket, i)
                    )
                    print('SMILES: ', smiles)
                    print('-----------------------------------------')
