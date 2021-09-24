"""
Generate pdbqt files from smi files
"""
import argparse
import os
from os import listdir
from os.path import isfile, join
import subprocess
from multiprocessing import Pool


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-num_workers",
                        required=False,
                        default=4,
                        help="number of workders to generate the pdbqt files")

    parser.add_argument("-in_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/zinc_sampled_smi/",
                        help="directory of smi files after being selected by clustering algorithm")

    parser.add_argument("-out_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/zinc_sampled_pdbqt/",
                        help="directory of pdbqt files of sampled SMILES")

    return parser.parse_args()


def worker_func(in_tuples):
    """
    The worker function for the multi-processing module.
    Input: a tuple containing 
        1. input directory
        2. output directory
        3. pocket_subset_file where each line is a SMILES
    """
    for in_tuple in in_tuples:
        # upack input
        in_dir, out_dir, pocket_subset_file = in_tuple

        pocket_name = pocket_subset_file.split('.')[0]

        # make a folder for each pocket
        pocket_out_dir = join(out_dir, pocket_name)
        if not os.path.exists(pocket_out_dir):
            os.makedirs(pocket_out_dir)

        # input .smi file name
        in_file = join(in_dir, pocket_subset_file)

        # run openbabel to convert files
        p = subprocess.run(
            'obabel' + 
            ' -ismi {}'.format(in_file) +
            ' -opdbqt' +
            ' -O {}'.format(join(pocket_out_dir, pocket_name + '-.pdbqt')) +
            ' -h -m --gen3D',
            shell=True,
            capture_output=True
        )

        if p.returncode != 0:
            print('Something went wrong! pocket: ', pocket_name)


if __name__ == "__main__":
    args = get_args()
    num_workers = int(args.num_workers)
    in_dir = args.in_dir
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    pocket_subset_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]

    # a list of input tuples to be processed by the worker function
    in_tuples = [(in_dir, out_dir, pocket_subset_file) for pocket_subset_file in pocket_subset_files]

    # divide the input tuples according to the number of workers
    in_tuple_list = []
    jobs_per_worker = int(len(in_tuples) / num_workers)
    for start in range(0, len(in_tuples), jobs_per_worker):        
        end = start + jobs_per_worker
        in_tuple_list.append(in_tuples[start:end])

    # parallel processing
    with Pool(num_workers) as p:
        p.map(worker_func, in_tuple_list)