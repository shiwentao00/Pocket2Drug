"""
Sample from a trained model for specified pockets, 
remove the duplicate pockets and save the freuency of each pocket. The output file is
a dictionary in yaml file.
"""
from dataloader import PocketDataset
from model import Pocket2Drug
import os
import yaml
import argparse
import torch
from torch_geometric.loader import DataLoader
from rdkit import Chem
from rdkit_contrib.sascorer import calculateScore

def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument("-result_dir",
                        required=True,
                        help="directory of result files including configuration, \
                            loss, trained model, and sampled molecules")

    parser.add_argument("-batch_size",
                        required=False,
                        default=32,
                        help="number of samples to generate per mini-batch")

    parser.add_argument("-num_batches",
                        required=False,
                        default=2,
                        help="number of batches to generate")

    parser.add_argument("-temperature",
                        required=False,
                        default=1.0,
                        help="the temperature paramter used to reshape softmax")
    
    parser.add_argument("-fold",
                        required=False,
                        default=42,
                        help="which fold to sample")
    
    parser.add_argument("-pocket_folds_dir",
                        required=False,
                        default="./data/folds/",
                        help="The directory of yaml files contains a dictionary of a fold of pockets \
                              to sample. The keys are the pocket names and the values are \
                              the target SMILES.")

    parser.add_argument("-pocket_dir",
                        required=False,
                        default="../test_data/pocket-data/",
                        help="the directory of pocket files")

    parser.add_argument("-popsa_dir",
                        required=False,
                        default="../test_data/protein-data/",
                        help="the directory of popsa files associated with \
                              the pockets")

    parser.add_argument("-profile_dir",
                        required=False,
                        default="../test_data/protein-data/",
                        help="the directory of profile files associated with \
                              the pockets")

    return parser.parse_args()


def compute_sa_score(mol):
    try:
        sa_score = calculateScore(mol)
    except:
        print("Something went wrong when computing SA.")
        sa_score = None

    return sa_score


if __name__ == "__main__":
    args = get_args()
    result_dir = args.result_dir
    batch_size = int(args.batch_size)
    num_batches = int(args.num_batches)
    fold = int(args.fold)
    pocket_folds_dir = args.pocket_folds_dir    
    temperature = float(args.temperature)
    pocket_dir = args.pocket_dir
    popsa_dir = args.popsa_dir
    profile_dir = args.profile_dir

    # directory to store the sampled molecules
    sampled_mols_dir = result_dir + f"/val_pockets_sample_{batch_size * num_batches}/"
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
    model_path = result_dir + "trained_model.pt"
    encoder_config = config['encoder_config']
    decoder_config = config['decoder_config']
    model = Pocket2Drug(encoder_config, decoder_config).to(device)
    model.load_state_dict(
        torch.load(
            model_path,
            map_location=torch.device(device)
        )
    )
    model.eval()

    # load the list of input pockets 
    if 0 <= fold <= 9:
        # input pockets are from the validation fold
        val_pockets = pocket_folds_dir + 'pockets_fold{}.yaml'.format(fold)
        with open(val_pockets, 'r') as f:
            val_pockets = yaml.full_load(f)
        val_pockets = list(val_pockets.keys())
    else:
        # input pockets are from user-defined directory
        val_pockets = [pocket.name for pocket in os.scandir(pocket_dir) if pocket.is_dir()]
    num_val_pockets = len(val_pockets)

    # load the dictionary of SMILES
    smiles_dir = config['smiles_dir']
    with open(smiles_dir, 'r') as f:
        smiles_dict = yaml.full_load(f)

    # create a valset of the pockets
    features_to_use = config['features_to_use']
    valset = PocketDataset(
        pockets=val_pockets,
        pocket_dir=pocket_dir,
        pop_dir=popsa_dir,
        profile_dir=profile_dir,
        smiles_dict=smiles_dict,
        features_to_use=features_to_use,
        vocab=vocab,
        vocab_path=vocab_path
    )

    # wrap the valset with the dataloader
    val_loader = DataLoader(
        valset,
        batch_size=1,
        shuffle=False,
        num_workers=1,
        drop_last=False
    )

    # sample SMILES for each pocket
    # dataloader has batch_size of 1
    for i, data in enumerate(val_loader):
        data = data.to(device)
        pocket_name = data.pocket_name[0]
        print('sampling SMILES for pocket {}...'.format(pocket_name))
        # target_smiles_ints = data.y[0]
        # target_smiles = []
        # for x in target_smiles_ints:
            # if valset.vocab.int2tocken[x] == '<eos>':
                # break
            # else:
                # target_smiles.append(valset.vocab.int2tocken[x])
        # target_smiles = "".join(target_smiles[1:])
        # print('target SMILES: ', target_smiles)

        # sample
        molecules = model.sample_from_pocket(
            data,
            num_batches,
            batch_size,
            temperature,
            valset.vocab,
            device
        )

        # a dictionary that stores the frequency of each valid SMILES
        mol_dict = {}

        # filter out invalid SMILES
        print('filtering out invalid sampled SMILES...')
        num_valid, num_invalid = 0, 0
        for smiles in molecules:
            if smiles is None:
                print('SMILES of None value in sample')
                num_invalid += 1
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                num_invalid += 1
            else:
                sa_score = compute_sa_score(mol)
                if sa_score and 1 <= sa_score <= 6:
                    if smiles in mol_dict:
                        mol_dict[smiles] += 1
                        num_valid += 1
                    else:
                        mol_dict[smiles] = 1
                        num_valid += 1
                else:
                    num_invalid += 1

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
