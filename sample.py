"""
Sample from the RNN of the model.
"""
import yaml
import argparse
import torch
from tqdm import tqdm
from torch_geometric.data import DataLoader

from model import Pocket2Drug
from dataloader import PocketDataset


def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument("-result_dir",
                        required=True,
                        help="directory of result files including configuration, \
                            loss, trained model, and sampled molecules"
                        )

    parser.add_argument("-batch_size",
                        required=False,
                        default=8,
                        help="number of samples to generate per mini-batch"
                        )

    parser.add_argument("-num_batches",
                        required=False,
                        default=2,
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

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    result_dir = args.result_dir
    batch_size = args.batch_size
    num_batches = args.num_batches
    pocket_list = args.pocket_list
    pocket_dir = args.pocket_dir
    popsa_dir = args.popsa_dir

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

    # sample
    for data in trainloader:
        data = data.to(device)
        pocket_name = data.pocket_name[0]
        print(
            'sampe SMILES for pocket {}...'.format(pocket_name)
        )
        smiles = data.y

        # sample
        molecules = model.sample_from_pocket(
            data, num_batches,
            batch_size, dataset.vocab,
            device
        )

        # save sampled SMILES
        out_path = result_dir + pocket_name + '_sampled.yaml'
        with open(out_path, 'w') as f:
            yaml.dump(molecules, f)
        break
