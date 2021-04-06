"""
Sample from the RNN of the model.
"""
from Pocket2Drug.model import Pocket2Drug
import argparse
import torch

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
                        default=1024,
                        help="number of samples to generate per mini-batch"
                        )

    parser.add_argument("-num_batches",
                        required=False,
                        default=10,
                        help="number of batches to generate"
                        )

    parser.add_argument("-pocket",
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
    pocket = args.pocket
    pocket_dir = args.pocket_dir
    popsa_dir = args.popsa_dir

    # load the model

    # load the pockets as graphs

    # sample
