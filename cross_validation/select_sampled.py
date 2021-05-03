"""
Using Rdkit to select a diverse subset of sampled molecules of 
a pocket in the validation set. 
"""
import argparse
import yaml

def get_args():
    parser = argparse.ArgumentParser("python")
    parser.add_argument("-mol_dir",
                        required=False,
                        default="../../p2d_results/cross_val_fold_0/val_pockets_sample/",
                        help="directory of yaml files where each yaml is a list of sampled SMILES")

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    mol_dir = args.mol_dir