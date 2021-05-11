"""
Compute the SELFIES vocabulary of the dataset. 
"""
import yaml
import selfies as sf
from tqdm import tqdm


def read_dataset(dataset_path):
    with open(dataset_path, "r") as f:
        smiles_dict = yaml.full_load(f)
    return smiles_dict


if __name__ == "__main__":
    smiles_dict = read_dataset("../data/pocket-smiles.yaml")
    smiles_list = smiles_dict.values()

    dataset = []
    for smiles in tqdm(smiles_list):
        selfies = sf.encoder(smiles)
        if selfies is not None:
            dataset.append(selfies)
    print('{} smiles converted to selfies.'.format(len(dataset)))

    vocab = list(sf.get_alphabet_from_selfies(dataset))
    print(vocab)

    vocab_dict = {}
    for i, token in enumerate(vocab):
        vocab_dict[token] = i

    i += 1
    vocab_dict['<eos>'] = i
    i += 1
    vocab_dict['<sos>'] = i
    i += 1
    vocab_dict['<pad>'] = i

    with open("./selfies_vocab.yaml", 'w') as f:
        yaml.dump(vocab_dict, f)
