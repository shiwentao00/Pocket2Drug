"""
Analyse the sampled SMILES of the case study pockets, 
draw figures to visualize the distributions of similarity
scores and frequencies of the pockets
"""
import yaml
import argparse
import random
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import FingerprintSimilarity
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np


# suppress rdkit error
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument("-result_dir",
                        required=True,
                        help="directory of result files including configuration, \
                            loss, trained model, and sampled molecules")

    parser.add_argument("-pocket_list",
                        required=False,
                        default="./data/excluded_pockets.yaml",
                        help="the pocket name used to condition the RNN. If None. \
                            a random vector will be used.")

    parser.add_argument("-pocket_dir",
                        required=False,
                        default="../data/googlenet-dataset/",
                        help="the directory of pocket files")

    parser.add_argument("-random_molecule_path",
                        required=False,
                        default="../mol-rnn-sampled/run_34/sampled_molecules.out",
                        help="the directory of popsa files associated with the pockets")

    return parser.parse_args()


def read_random_pockets(mol_path, num):
    """read num pockets randomly from the input file"""
    with open(mol_path, 'r') as f:
        mols = f.readlines()
    mols = [x.strip('\n') for x in mols]
    mols = list(set(mols))
    assert (num <= len(mols))
    sampled = random.sample(mols, num)
    return sampled


def compute_similarity(target_maccs, mol):
    """compute the similarity between two smiles"""
    mol_maccs = MACCSkeys.GenMACCSKeys(mol)
    return FingerprintSimilarity(target_maccs, mol_maccs)


def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    # ax.set_xlabel('Sample name')


def compute_all_similarities(sampled_smiles, target_maccs):
    similarities = []
    unique_similarities = []
    for smiles in tqdm(sampled_smiles):
        # compute the simlarity
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            s = compute_similarity(target_maccs, mol)
            # if there are multiple occurences
            if sampled_smiles[smiles] > 0:
                similarities.extend([s] * sampled_smiles[smiles])
            else:
                similarities.append(s)
            unique_similarities.append(s)
    return similarities, unique_similarities


if __name__ == "__main__":
    args = get_args()
    result_dir = args.result_dir
    pocket_list = args.pocket_list
    pocket_dir = args.pocket_dir
    random_molecule_path = args.random_molecule_path

    # load the configuartion file in output
    config_dir = result_dir + "config.yaml"
    with open(config_dir, 'r') as f:
        config = yaml.full_load(f)

    # load vocab
    vocab, vocab_path = config["vocab"], config["vocab_path"]

    # load the list of pockets
    with open(pocket_list, 'r') as f:
        test_pockets = yaml.full_load(f)

    # load the dictionary of SMILES
    smiles_dir = config['smiles_dir']
    with open(smiles_dir, 'r') as f:
        smiles_dict = yaml.full_load(f)

    # for each case study pocket
    for pocket in test_pockets:
        print("plotting figures for pocket {}...".format(pocket))
        # get the target smiles
        target_smiles = smiles_dict[pocket]
        target_mol = Chem.MolFromSmiles(target_smiles)
        target_maccs = MACCSkeys.GenMACCSKeys(target_mol)

        # load the sampled smiles of the pocket
        sampled_smiles_path = result_dir + pocket + '_sampled.yaml'
        with open(sampled_smiles_path, 'r') as f:
            sampled_smiles = yaml.full_load(f)
        print("{} unique SMILES sampled".format(len(sampled_smiles.keys())))

        # load same number of unique random smiles
        random_molecules = read_random_pockets(
            random_molecule_path,
            len(sampled_smiles.keys())
        )

        # similarities of sampled molecules
        print("computing similarities of sampled molecules with target ligand...")
        similarities, unique_similarities = compute_all_similarities(
            sampled_smiles, target_maccs)

        # similarities of random SMILES
        print('computing similarities of random moldecules with target ligand...')
        random_similarity = []
        for smiles in tqdm(random_molecules):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol_maccs = MACCSkeys.GenMACCSKeys(mol)
                random_similarity.append(
                    FingerprintSimilarity(target_maccs, mol_maccs)
                )

        # similarites sampled at temperature 4.0
        sampled_smiles_path = result_dir + pocket + '_sampled_temp4.0.yaml'
        with open(sampled_smiles_path, 'r') as f:
            sampled_smiles = yaml.full_load(f)
        print("{} unique SMILES sampled".format(len(sampled_smiles.keys())))
        print(
            "computing similarities of sampled molecules with target ligand at temp 4.0...")
        similarities_4, unique_similarities_4 = compute_all_similarities(
            sampled_smiles, target_maccs)

        # plot the distributions of the similarites
        fig, ax = plt.subplots(1, 1)
        # ax.hist(
        # similarities,
        # bins=100,
        # label='SMILES based on pocket',
        # histtype='step',
        # stacked=True,
        # fill=False
        # )
        ax.hist(
            random_similarity,
            bins=100,
            label='Random',
            histtype='step',
            stacked=True,
            fill=False
        )
        ax.hist(
            unique_similarities,
            bins=100,
            label='Sampled & Unique (temp=1.0)',
            histtype='step',
            stacked=True,
            fill=False
        )
        ax.hist(
            unique_similarities_4,
            bins=100,
            label='Sampled & Unique (temp=4.0)',
            histtype='step',
            stacked=True,
            fill=False
        )
        ax.legend(loc='upper right')
        ax.set_xlabel('Tanimoto Similarity with Target Ligand')
        ax.set_ylabel('Ligand Count')
        ax.set_title('Pocket: {}'.format(pocket))
        plt.savefig(
            result_dir + "{}_similarity_temp.png".format(pocket),
            dpi=150,
            bbox_inches='tight'
        )

        # violin plot of the distributions of the similarites
        fig, ax = plt.subplots(1, 1)
        ax.violinplot(random_similarity, positions=[1])
        ax.violinplot(unique_similarities, positions=[2])
        ax.violinplot(unique_similarities_4, positions=[3])
        ax.violinplot(similarities, positions=[4])
        ax.violinplot(similarities_4, positions=[5])
        # ax.legend(loc='upper left')
        # ax.set_xlabel('Similarities')
        labels = ['Random',
                  'Unique-1.0', 'Unique-4.0',
                  'All-1.0', 'All-4.0']
        set_axis_style(ax, labels)
        ax.set_ylabel('Tanimoto Similarity')
        ax.set_title('Pocket: {}'.format(pocket))
        plt.savefig(
            result_dir + "{}_similarity_violin_temp.png".format(pocket),
            dpi=150,
            bbox_inches='tight'
        )

        # box plot of the distributions of the similarites
        fig, ax = plt.subplots(1, 1)
        ax.boxplot(random_similarity, sym='g.', positions=[1])
        ax.boxplot(unique_similarities, sym='g.', positions=[2])
        ax.boxplot(unique_similarities_4, sym='g.', positions=[3])
        ax.boxplot(similarities, sym='g.', positions=[4])
        ax.boxplot(similarities_4, sym='g.', positions=[5])
        # ax.legend(loc='upper left')
        # ax.set_xlabel('Similarities')
        labels = ['Random',
                  'Unique-1.0', 'Unique-4.0',
                  'All-1.0', 'All-4.0']
        set_axis_style(ax, labels)
        ax.set_ylabel('Tanimoto Similarity')
        ax.set_title('Pocket: {}'.format(pocket))
        plt.savefig(
            result_dir + "{}_similarity_box_temp.png".format(pocket),
            dpi=150,
            bbox_inches='tight'
        )
