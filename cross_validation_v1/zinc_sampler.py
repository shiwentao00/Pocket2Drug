"""
For the subset of drugs of each pocket, sample the same number of molecules
from the Zinc dataset, such that the molecules from the Zinc dataset have
the same molecular weight distribution.
"""
import os
from os import listdir
from os.path import isdir, isfile, join
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
import yaml
import pickle
import random
import argparse
from tqdm import tqdm 


class ZincDistributor:
    """
    Read SMILES from the Zinc dataset, then distribute the smiles into 9
    bins according to their molecular weights. Those bins are the same
    as the Zinc tranches: https://zinc.docking.org/tranches/home/#. The
    bins will then be saved as a yaml file
    """
    def __init__(self, zinc_path):
        # bins to store smiles according to molecular weights
        self.bins = [[] for _ in range(9)]

        # the tranch folders
        zinc_folders = []
        for child in listdir(zinc_path):
            d = join(zinc_path, child)
            if isdir(d):
                zinc_folders.append(d)
        
        # load smiles into memory
        print('loading Zinc SMILES into memory...')
        for zinc_folder in tqdm(zinc_folders):
            smiles_files = [f for f in listdir(zinc_folder) if isfile(join(zinc_folder, f))]
            for smiles_file in smiles_files:
                self.read_mols(join(zinc_folder, smiles_file))

        # test correctness of bin 0
        for smiles in tqdm(self.bins[0]):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol_weight = Descriptors.MolWt(mol)
                if mol_weight >= 250:
                    print("something went wrong:")
                    print(smiles, mol_weight)

        print('all the SMILES are now distributed in bins')
        print('sizes of the bins: ', [len(b) for b in self.bins])
        print('total number of SMILES: ', sum([len(b) for b in self.bins]))

    def read_mols(self, smi_path):
        """
        Read a .smi file, get the SMILES, compute the molecular weights, and
        distribute the SMILES into the bins
        """
        with open(smi_path) as fp:
            lines = fp.readlines()
        
        smiles_list = [line.split(' ')[0] for line in lines]

        # The molecules in the same file are from the same tranch,
        # so they should fall into the same bin. Therefore, we
        # only need to calculate the weight of only one molecule.
        for smiles in smiles_list[1:]:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol_weight = Descriptors.MolWt(mol)

                # the range of bins are the same as Zinc tranches:
                # https://zinc.docking.org/tranches/home/#
                if 200 <= mol_weight < 250:
                    self.bins[0].append(smiles)
                elif 250 <= mol_weight < 300:
                    self.bins[1].append(smiles)
                elif 300 <= mol_weight < 325:
                    self.bins[2].append(smiles)
                elif 325 <= mol_weight < 350:
                    self.bins[3].append(smiles)
                elif 350 <= mol_weight < 375:
                    self.bins[4].append(smiles)
                elif 375 <= mol_weight < 400:
                    self.bins[5].append(smiles)
                elif 400 <= mol_weight < 425:
                    self.bins[6].append(smiles)
                elif 425 <= mol_weight < 450:
                    self.bins[7].append(smiles)
                elif 450 <= mol_weight < 500:
                    self.bins[8].append(smiles)

    def save_bins(self, out_path):
        """
        Save the bins of SMILES into a yaml file.
        """
        print(f'saving the bins to {out_path}...')
        with open(out_path, "wb") as f:
            #yaml.dump(self.bins, f)
            pickle.dump(self.bins, f)

class ZincSampler:
    """
    Zinc Sampler is used to sample a set of molecules from the Zinc dataset.
    The sampling function takes a list of molecular weights as input, and 
    returns a list of smiles
    """
    def __init__(self, zinc_path):
        # load bins of molecules
        with open(zinc_path, 'rb') as f:
            #bins = yaml.full_load(f)
            self.bins = pickle.load(f)

    def sample_mols(self, smiles_list):
        """
        Sample from the bins according to the input molecular weights.
        mol_weights - list of integers
        """
        # compute the molecular weights
        mol_weights = smiles_to_weights(smiles_list)

        # calculate the count of each bin
        counts = count_mol_weights(mol_weights)

        # sample from each bin according to the counts
        sampled = []
        for i, cnt in enumerate(counts):
            sampled.extend(random.sample(self.bins[i], cnt))

        return sampled
        

def count_mol_weights(mol_weights):
    """
    Count the number of weights in each bin, and 
    return a list of counts.
    """
    counts = [0 for _ in range(9)]

    for mol_weight in mol_weights:
        if mol_weight < 250:
            counts[0] += 1    
        elif 250 <= mol_weight < 300:
            counts[1] += 1
        elif 300 <= mol_weight < 325:
            counts[2] += 1
        elif 325 <= mol_weight < 350:
            counts[3] += 1
        elif 350 <= mol_weight < 375:
            counts[4] += 1
        elif 375 <= mol_weight < 400:
            counts[5] += 1
        elif 400 <= mol_weight < 425:
            counts[6] += 1
        elif 425 <= mol_weight < 450:
            counts[7] += 1
        else:
            counts[8] += 1

    return counts

def smiles_to_weights(smiles_list):
    """
    Compute the molecular weights given a list of SMILES.
    """
    molecular_weights = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_weight = Descriptors.MolWt(mol)
            molecular_weights.append(mol_weight)

    return molecular_weights


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-zinc_dataset",
                        required=False,
                        default="../../../zinc-sampler-project/zinc-dataset-bins.pickle",
                        help="the SMILES from Zinc in bins")

    parser.add_argument("-in_dir",
                        required=False,
                        default="../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample/",
                        help="input molecules for molecular weight distribution")

    parser.add_argument("-out_dir",
                        required=False,
                        default="../../p2d_results_selfie/cv_results/cross_val_fold_0/zinc_rank_sampled/",
                        help="directory for sampled SMILES from Zinc")

    return parser.parse_args()


if __name__=="__main__":
    args = get_args()
    zinc_dataset = args.zinc_dataset
    in_dir = args.in_dir
    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load the Zinc molecules in bins to the sampler
    sampler = ZincSampler(zinc_dataset)

    # the yaml file list of the molecules generated by the GNN model for each pocket
    mol_files = [f for f in listdir(in_dir) if isfile(join(in_dir, f))]

    # read each file and sample according to the molecular weight distribution
    for molFile in tqdm(mol_files):
        with open(join(in_dir, molFile), 'r') as f:
            in_mols = yaml.full_load(f)

        # sort the SMILES according to their frequencies
        in_mols = [(freq, smiles) for smiles, freq in in_mols.items()]
        in_mols.sort(key=lambda x: x[0], reverse=True)
        in_mols = in_mols[0:100]
        in_mols = [x[1] for x in in_mols]
        
        # sample
        out_mols = sampler.sample_mols(in_mols)

        out_file = molFile.split('_')[0] + '_zinc.yaml'
        with open(join(out_dir, out_file), 'w') as f:
            yaml.dump(out_mols, f)