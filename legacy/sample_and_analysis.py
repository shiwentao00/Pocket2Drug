# Copyright: Wentao Shi, 2021
"""
Sample from the RNN of the model.
"""
from dataloader import PocketDataset
from model import Pocket2Drug
import yaml
import argparse
import random
import torch
from torch_geometric.data import DataLoader
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import FingerprintSimilarity
from tqdm import tqdm
import matplotlib.pyplot as plt

# suppress rdkit error
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument("-result_dir",
                        required=True,
                        help="directory of result files including configuration, \
                            loss, trained model, and sampled molecules"
                        )

    parser.add_argument("-batch_size",
                        required=False,
                        default=2048,
                        help="number of samples to generate per mini-batch"
                        )

    parser.add_argument("-num_batches",
                        required=False,
                        default=200,
                        help="number of batches to generate"
                        )

    parser.add_argument("-pocket_list",
                        required=False,
                        default="./data/excluded_pockets.yaml",
                        help="the pocket name used to condition the RNN. If None. \
                            a random vector will be used."
                        )

    parser.add_argument("-pocket_dir",
                        required=False,
                        default="../data/googlenet-dataset/",
                        help="the directory of pocket files"
                        )

    parser.add_argument("-popsa_dir",
                        required=False,
                        default="../data/pops-googlenet/",
                        help="the directory of popsa files associated with the pockets"
                        )

    parser.add_argument("-random_molecule_path",
                        required=False,
                        default="../mol-rnn-sampled/run_34/sampled_molecules.out",
                        help="the directory of popsa files associated with the pockets"
                        )

    return parser.parse_args()


def read_random_pockets(mol_path, num):
    """read num pockets randomly from the input file"""
    with open(mol_path, 'r') as f:
        mols = f.readlines()
    mols = [x.strip('\n') for x in mols]
    # print('number of molecules from Chembl28 RNN: ', len(mols))
    mols = list(set(mols))
    # print('number of molecules after removing duplicates: ', len(mols))
    assert (num <= len(mols))
    sampled = random.sample(mols, num)
    return sampled


def compute_similarity(target_maccs, mol):
    """compute the similarity between two smiles"""
    mol_maccs = MACCSkeys.GenMACCSKeys(mol)
    return FingerprintSimilarity(target_maccs, mol_maccs)


if __name__ == "__main__":
    args = get_args()
    result_dir = args.result_dir
    batch_size = int(args.batch_size)
    num_batches = int(args.num_batches)
    pocket_list = args.pocket_list
    pocket_dir = args.pocket_dir
    popsa_dir = args.popsa_dir
    random_molecule_path = args.random_molecule_path

    # load the configuartion file in output
    config_dir = result_dir + "config.yaml"
    with open(config_dir, 'r') as f:
        config = yaml.full_load(f)

    # detect cpu or gpu
    device = torch.device(
        'cuda' if torch.cuda.is_available() else 'cpu'
    )
    print('device: ', device)

    # load vocab
    vocab, vocab_path = config["vocab"], config["vocab_path"]

    # load the model
    encoder_config = config['encoder_config']
    decoder_config = config['decoder_config']
    model = Pocket2Drug(encoder_config, decoder_config).to(device)
    model.load_state_dict(
        torch.load(
            config['out_dir'] + 'trained_model.pt',
            map_location=torch.device(device)
        )
    )
    model.eval()

    # load the list of pockets
    with open(pocket_list, 'r') as f:
        test_pockets = yaml.full_load(f)

    # load the dictionary of SMILES
    smiles_dir = config['smiles_dir']
    with open(smiles_dir, 'r') as f:
        smiles_dict = yaml.full_load(f)

    # create a dataset of the case study pockets
    features_to_use = config['features_to_use']
    dataset = PocketDataset(
        pockets=test_pockets,
        pocket_dir=pocket_dir,
        pop_dir=popsa_dir,
        smiles_dict=smiles_dict,
        features_to_use=features_to_use,
        vocab=vocab,
        vocab_path=vocab_path
    )

    # wrap the dataset with the dataloader
    trainloader = DataLoader(
        dataset,
        batch_size=1,
        shuffle=False,
        num_workers=1,
        drop_last=False
    )

    # sample SMILES for each pocket
    # dataloader has batch_size of 1
    for data in trainloader:
        data = data.to(device)
        pocket_name = data.pocket_name[0]
        print('sampling SMILES for pocket {}...'.format(pocket_name))
        target_smiles_ints = data.y[0]
        target_smiles = []
        for x in target_smiles_ints:
            if dataset.vocab.int2tocken[x] == '<eos>':
                break
            else:
                target_smiles.append(dataset.vocab.int2tocken[x])
        target_smiles = "".join(target_smiles[1:])
        print('target SMILES: ', target_smiles)

        # sample
        molecules = model.sample_from_pocket(
            data, num_batches,
            batch_size, dataset.vocab,
            device
        )

        # remove duplicate samples
        print('removing duplicates...')
        original_num = len(molecules)
        molecules = list(set(molecules))
        print('unique ratio after sampling: {}%'.format(
            len(molecules) * 100.0 / original_num
        ))

        # filter out invalid SMILES
        print('filtering out invalid sampled SMILES and computing similarities...')
        num_valid, num_invalid = 0, 0
        valid_molecules = []
        sampled_similarity = []
        target_mol = Chem.MolFromSmiles(target_smiles)
        target_maccs = MACCSkeys.GenMACCSKeys(target_mol)
        for smiles in tqdm(molecules):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                num_invalid += 1
            else:
                num_valid += 1
                valid_molecules.append(smiles)
                sampled_similarity.append(
                    compute_similarity(target_maccs, mol)
                )
        print('sampled {} unique SMILES, {}% of which are valid'.format(
            num_valid + num_invalid,
            num_valid * 100.0 / (num_valid + num_invalid)
        ))

        # save sampled SMILES
        out_path = result_dir + pocket_name + '_sampled.yaml'
        with open(out_path, 'w') as f:
            yaml.dump(valid_molecules, f)

        # get same number of random molecules as control group
        random_molecules = read_random_pockets(
            random_molecule_path,
            len(valid_molecules)
        )

        # compute the Tanimoto coefficients of
        # random molecules with the target ligand
        print('computing similarities of random moldecules with target ligand...')
        random_similarity = []
        for smiles in tqdm(random_molecules):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol_maccs = MACCSkeys.GenMACCSKeys(mol)
                random_similarity.append(
                    FingerprintSimilarity(target_maccs, mol_maccs)
                )

        # plot the distributions of the Tanimoto coefficients
        fig, ax = plt.subplots(1, 1)
        ax.hist(
            sampled_similarity,
            bins=20,
            label='model based on pocket',
            histtype='step',
            stacked=True,
            fill=False
        )
        ax.hist(
            random_similarity,
            bins=20,
            label='model based on chembl28',
            histtype='step',
            stacked=True,
            fill=False
        )
        ax.legend(loc='upper left')
        ax.set_xlabel('Tanimoto Similarity with Target Ligand')
        ax.set_ylabel('Ligand Count')
        ax.set_title('Pocket name: {}'.format(pocket_name))
        plt.savefig(
            result_dir + "{}_similarity.png".format(pocket_name),
            dpi=150
        )
