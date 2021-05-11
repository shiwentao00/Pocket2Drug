"""
Merge the SELFIES vocab of two datasets.
"""
import yaml


if __name__ == "__main__":
    dataset_0_path = "./selfies_vocab.yaml"
    dataset_1_path = "../../../molecule-generator-project/Molecule-RNN/vocab/chembl_selfies_vocab.yaml"

    with open(dataset_0_path, 'r') as f:
        dataset_0 = yaml.full_load(f)

    with open(dataset_1_path, 'r') as f:
        dataset_1 = yaml.full_load(f)

    vocab_0 = set(list(dataset_0.keys()))
    vocab_1 = set(list(dataset_1.keys()))
    print('size of vocab 0', len(vocab_0))
    print('size of vocab 1', len(vocab_1))

    vocab = vocab_0.union(vocab_1)
    vocab.remove('<eos>')
    vocab.remove('<sos>')
    vocab.remove('<pad>')

    vocab_dict = {}
    for i, token in enumerate(list(vocab)):
        vocab_dict[token] = i

    i += 1
    vocab_dict['<eos>'] = i
    i += 1
    vocab_dict['<sos>'] = i
    i += 1
    vocab_dict['<pad>'] = i

    with open("./selfies_merged_vocab.yaml", 'w') as f:
        yaml.dump(vocab_dict, f)


    
