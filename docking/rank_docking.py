"""
Rank the SMILES according to docking scores
"""
import argparse
import yaml


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-docking_result",
                        required=False,
                        default="../../p2d_results/run_9/selective_docking_2hs0A004.0/2hs0A00_docking_score_th0.94.yaml",
                        help="directory sampled SMILES and their docking scores")

    parser.add_argument("-smiles_dict",
                        required=False,
                        default="../../p2d_results/run_9/2hs0A00_sampled_temp4.0.yaml",
                        help="directory sampled SMILES and their frequencies")
                        
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    docking_result = args.docking_result
    smiles_dict = args.smiles_dict

    with open(docking_result, 'r') as f:
        docking_scores = yaml.full_load(f)

    with open(smiles_dict, 'r') as f:
        smiles_freq = yaml.full_load(f)

    sorted_docking_scores = sorted(docking_scores.items(), key=lambda item: item[1])
    for smiles, score in sorted_docking_scores:    
        print('{}: {}, {}'.format(smiles, score, smiles_freq[smiles]))

