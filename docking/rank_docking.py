"""
Rank the SMILES according to docking scores
"""
import argparse
import yaml


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-docking_result",
                        required=True,
                        help="directory of files of sampled SMILES")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    docking_result = args.docking_result
    with open(docking_result, 'r') as f:
        docking_scores = yaml.full_load(f)

    sorted_docking_scores = sorted(docking_scores.items(), key=lambda item: item[1])
    for smiles, score in sorted_docking_scores:    
        print('{}: {}'.format(smiles, score))

