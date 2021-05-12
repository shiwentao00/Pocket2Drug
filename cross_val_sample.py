# Copyright: Wentao Shi, 2021
"""
Sample from a trained model and a given pocket, 
remove the duplicate pockets and save the 
freuency of each pocket. The output file is
a dictionary in yaml format.
"""
from dataloader import PocketDataset
from model import Pocket2Drug
import os
import yaml
import argparse
import torch
from torch_geometric.data import DataLoader
from rdkit import Chem


def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument("-result_dir",
                        required=True,
                        help="directory of result files including configuration, \
                            loss, trained model, and sampled molecules")

    parser.add_argument("-batch_size",
                        required=False,
                        default=2048,
                        help="number of samples to generate per mini-batch")

    parser.add_argument("-num_batches",
                        required=False,
                        default=10,
                        help="number of batches to generate")

    parser.add_argument("-temperature",
                        required=False,
                        default=1.0,
                        help="the temperature paramter used to reshape softmax")
    
    parser.add_argument("-fold",
                        required=False,
                        default=0,
                        help="which fold to sample")
    
    parser.add_argument("-pocket_folds_dir",
                        required=False,
                        default="./cross_validation/folds/",
                        help="The directory of yaml files contains a dictionary of a fold of pockets \
                              to sample. The keys are the pocket names and the values are \
                              the target SMILES.")

    parser.add_argument("-pocket_dir",
                        required=False,
                        default="../data/googlenet-dataset/",
                        help="the directory of pocket files")

    parser.add_argument("-popsa_dir",
                        required=False,
                        default="../data/pops-googlenet/",
                        help="the directory of popsa files associated with \
                              the pockets")

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    result_dir = args.result_dir
    batch_size = int(args.batch_size)
    num_batches = int(args.num_batches)
    fold = int(args.fold)
    temperature = float(args.temperature)
    pocket_folds_dir = args.pocket_folds_dir
    pocket_dir = args.pocket_dir
    popsa_dir = args.popsa_dir

    # directory to store the sampled molecules
    sampled_mols_dir = result_dir + '/val_pockets_sample/'
    if not os.path.exists(sampled_mols_dir):
        os.makedirs(sampled_mols_dir)

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

    # load the list of pockets in the validation fold
    val_pockets = pocket_folds_dir + 'pockets_fold{}.yaml'.format(fold)
    with open(val_pockets, 'r') as f:
        val_pockets = yaml.full_load(f)
    val_pockets = list(val_pockets.keys())
    num_val_pockets = len(val_pockets)

    # load the dictionary of SMILES
    smiles_dir = config['smiles_dir']
    with open(smiles_dir, 'r') as f:
        smiles_dict = yaml.full_load(f)

    # create a dataset of the case study pockets
    features_to_use = config['features_to_use']
    dataset = PocketDataset(
        pockets=val_pockets,
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
    for i, data in enumerate(trainloader):
        data = data.to(device)
        pocket_name = data.pocket_name[0]
        print('sampling SMILES for pocket {}...'.format(pocket_name))
        # target_smiles_ints = data.y[0]
        # target_smiles = []
        # for x in target_smiles_ints:
            # if dataset.vocab.int2tocken[x] == '<eos>':
                # break
            # else:
                # target_smiles.append(dataset.vocab.int2tocken[x])
        # target_smiles = "".join(target_smiles[1:])
        # print('target SMILES: ', target_smiles)

        # sample
        molecules = model.sample_from_pocket(
            data,
            num_batches,
            batch_size,
            temperature,
            dataset.vocab,
            device
        )

        # a dictionary that stores the frequency of each valid SMILES
        mol_dict = {}

        # filter out invalid SMILES
        print('filtering out invalid sampled SMILES...')
        num_valid, num_invalid = 0, 0
        for smiles in molecules:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                num_invalid += 1
            elif smiles in mol_dict:
                mol_dict[smiles] += 1
                num_valid += 1
            else:
                mol_dict[smiles] = 1
                num_valid += 1

        valid_rate = float(num_valid / (num_valid + num_invalid))
        print("valid rate: {}%".format(valid_rate * 100))

        num_unique = len(list(mol_dict.keys()))
        unique_rate = float(num_unique / (num_valid + num_invalid))
        print("unique valid rate {}%".format(unique_rate * 100))

        print('{}/{} pockets sampled.'.format(i + 1, num_val_pockets))

        # save sampled SMILES
        out_path = sampled_mols_dir + pocket_name + \
            '_sampled_temp{}.yaml'.format(temperature)
        with open(out_path, 'w') as f:
            yaml.dump(mol_dict, f)