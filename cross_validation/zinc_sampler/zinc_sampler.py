import os
from os import listdir
from os.path import isdir, isfile, join
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
import yaml
import pickle
import random


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