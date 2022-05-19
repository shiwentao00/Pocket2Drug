# Copyright: Wentao Shi, 2021
import numpy as np
import pandas as pd
import statistics
import random

import torch
from torch_geometric.data import Data, Dataset
#from torch_geometric.data import DataLoader
from torch_geometric.loader import DataLoader
from biopandas.mol2 import PandasMol2
from scipy.spatial import distance

# self-defined utilities
import utils


def pocket_single_loader_gen(smiles_dict,
                             pocket_dir,
                             pop_dir,
                             profile_dir,
                             features_to_use,
                             vocab,
                             vocab_path,
                             batch_size,
                             shuffle=False,
                             num_workers=1):
    """
    Dataloader used to wrap PocketDataset. Generate only one dataloader

    Arguments:
        smiles_dict - a python dictionary of pocket-smiles pairs
        pocket_dir - root directory of the pockets
        pop_dir - root directory of the popsa files
        profile_dir - root directory of the profile files
        features_to_use - which node features to use
        vocab - which vocabulary to use
        vocab_path - path to load the vocabular
        batch_size - size of the mini-batch for training
        shuffle - whether to shuffle the dataset druing training
        test_split - ratio of test data of entire dataset
        num_workers - number of worker threads to load the data
    """
    # split pockets into train/test split
    pockets = list(smiles_dict.keys())
    random.shuffle(pockets)

    dataset = PocketDataset(
        pockets=pockets,
        pocket_dir=pocket_dir,
        pop_dir=pop_dir,
        profile_dir=profile_dir,
        smiles_dict=smiles_dict,
        features_to_use=features_to_use,
        vocab=vocab,
        vocab_path=vocab_path
    )

    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        drop_last=False
    )

    return dataloader, len(dataset)

def pocket_loader_gen(smiles_dict,
                      pocket_dir,
                      pop_dir,
                      profile_dir,
                      features_to_use,
                      vocab,
                      vocab_path,
                      batch_size,
                      shuffle=False,
                      test_split=0.1,
                      num_workers=1):
    """
    Dataloader used to wrap PocketDataset. Generate both train and validation
    dataloaders.

    Arguments:
        smiles_dict - a python dictionary of pocket-smiles pairs
        pocket_dir - root directory of the pockets
        pop_dir - root directory of the popsa files
        profile_dir - root directory of the profile files
        features_to_use - which node features to use
        vocab - which vocabulary to use
        vocab_path - path to load the vocabular
        batch_size - size of the mini-batch for training
        shuffle - whether to shuffle the dataset druing training
        test_split - ratio of test data of entire dataset
        num_workers - number of worker threads to load the data
    """
    # split pockets into train/test split
    pockets = list(smiles_dict.keys())
    random.shuffle(pockets)
    num_pockets = len(pockets)
    num_test_pockets = int(num_pockets * test_split)
    test_pockests = pockets[0:num_test_pockets]
    train_pockets = pockets[num_test_pockets:]

    trainset = PocketDataset(
        pockets=train_pockets,
        pocket_dir=pocket_dir,
        pop_dir=pop_dir,
        profile_dir=profile_dir,
        smiles_dict=smiles_dict,
        features_to_use=features_to_use,
        vocab=vocab,
        vocab_path=vocab_path
    )

    valset = PocketDataset(
        pockets=test_pockests,
        pocket_dir=pocket_dir,
        pop_dir=pop_dir,
        profile_dir=profile_dir,
        smiles_dict=smiles_dict,
        features_to_use=features_to_use,
        vocab=vocab,
        vocab_path=vocab_path
    )

    trainloader = DataLoader(
        trainset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        drop_last=False
    )

    valloader = DataLoader(
        valset,
        batch_size=1,
        shuffle=False,
        num_workers=num_workers,
        drop_last=False
    )
    return trainloader, valloader, len(trainset), len(valset)


class PocketDataset(Dataset):
    """Dataset to generate single pocket graphs for inference/testing."""

    def __init__(self,
                 pockets,
                 pocket_dir,
                 pop_dir,
                 profile_dir,
                 smiles_dict,
                 features_to_use,
                 vocab,
                 vocab_path):
        """
        pocket: a list of pockets to include.
        pockect_dir: directoy of pockets and the .profile files.
        pop_dir: directory of the sasa files.
        profile_dir - directory of the profile files
        smiles_dict: a python dictionary of pocket-smiles pairs
        features_to_use: which node features to use. 
        vocab: which vocabulary to use. options: ['char', 'selfies']
        vocab_path: path of the pre-computed vocab file
        """
        self.pockets = pockets
        self.pocket_dir = pocket_dir
        self.pop_dir = pop_dir
        self.profile_dir = profile_dir
        self.smiles_dict = smiles_dict

        # distance threshold to form an undirected edge between two atoms
        self.threshold = 4.5

        # hard coded info to generate 2 node features
        self.hydrophobicity = {'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5,
                               'CYS': 2.5, 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
                               'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9,
                               'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8,
                               'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2}
        self.binding_probability = {'ALA': 0.701, 'ARG': 0.916, 'ASN': 0.811, 'ASP': 1.015,
                                    'CYS': 1.650, 'GLN': 0.669, 'GLU': 0.956, 'GLY': 0.788,
                                    'HIS': 2.286, 'ILE': 1.006, 'LEU': 1.045, 'LYS': 0.468,
                                    'MET': 1.894, 'PHE': 1.952, 'PRO': 0.212, 'SER': 0.883,
                                    'THR': 0.730, 'TRP': 3.084, 'TYR': 1.672, 'VAL': 0.884}

        total_features = ['x', 'y', 'z', 'r', 'theta', 'phi', 'sasa', 'charge',
                          'hydrophobicity', 'binding_probability', 'sequence_entropy']

        # features to use should be subset of total_features
        assert(set(features_to_use).issubset(set(total_features)))
        self.features_to_use = features_to_use

        # initialize the vocabulary used to tokenize smiles
        if vocab == 'char':
            self.vocab = utils.CharVocab(vocab_path)
        elif vocab == 'selfies':
            self.vocab = utils.SELFIESVocab(vocab_path)
        elif vocab == 'regex':
            self.vocab = utils.RegExVocab(vocab_path)
        else:
            raise ValueError("invalid vocab value.")

    def __len__(self):
        return len(self.pockets)

    def __getitem__(self, idx):
        pocket = self.pockets[idx]

        # form the graph data
        pocket_dir = self.pocket_dir + pocket + '/' + pocket + '.mol2'
        profile_dir = f'{self.profile_dir}{pocket[0:-2]}/{pocket[0:-2]}.profile'
        pop_dir = f'{self.pop_dir}{pocket[0:-2]}/{pocket[0:-2]}.pops'

        x, edge_index, edge_attr = read_pocket(
            pocket_dir,
            profile_dir,
            pop_dir,
            self.hydrophobicity,
            self.binding_probability,
            self.features_to_use,
            self.threshold
        )
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

        if self.smiles_dict is not None:
            # read the smile data
            smile = self.smiles_dict[pocket]
    
            # convert the smiles to integers according to vocab
            smile = self.vocab.tokenize_smiles(smile)
            data.y = smile

        # save the pocket name in data
        data.pocket_name = pocket
        return data


def read_pocket(mol_path,
                profile_path,
                pop_path,
                hydrophobicity,
                binding_probability,
                features_to_use,
                threshold):
    """
    Read the mol2 file as a dataframe.
    """
    atoms = PandasMol2().read_mol2(mol_path)
    atoms = atoms.df[['atom_id', 'subst_name',
                      'atom_type', 'atom_name', 'x', 'y', 'z', 'charge']]
    atoms['residue'] = atoms['subst_name'].apply(lambda x: x[0:3])
    atoms['hydrophobicity'] = atoms['residue'].apply(
        lambda x: hydrophobicity[x])
    atoms['binding_probability'] = atoms['residue'].apply(
        lambda x: binding_probability[x])

    r, theta, phi = compute_spherical_coord(atoms[['x', 'y', 'z']].to_numpy())
    if 'r' in features_to_use:
        atoms['r'] = r
    if 'theta' in features_to_use:
        atoms['theta'] = theta
    if 'phi' in features_to_use:
        atoms['phi'] = phi

    siteresidue_list = atoms['subst_name'].tolist()

    if 'sasa' in features_to_use:
        qsasa_data = extract_sasa_data(siteresidue_list, pop_path)
        atoms['sasa'] = qsasa_data

    if 'sequence_entropy' in features_to_use:
        # sequence entropy data with subst_name as keys
        seq_entropy_data = extract_seq_entropy_data(
            siteresidue_list, profile_path)
        atoms['sequence_entropy'] = atoms['subst_name'].apply(
            lambda x: seq_entropy_data[x])

    if atoms.isnull().values.any():
        print('invalid input data (containing nan):')
        print(mol_path)

    bonds = bond_parser(mol_path)
    node_features, edge_index, edge_attr = form_graph(
        atoms, bonds, features_to_use, threshold)

    return node_features, edge_index, edge_attr


def bond_parser(pocket_path):
    f = open(pocket_path, 'r')
    f_text = f.read()
    f.close()
    bond_start = f_text.find('@<TRIPOS>BOND')
    bond_end = -1
    df_bonds = f_text[bond_start:bond_end].replace('@<TRIPOS>BOND\n', '')
    df_bonds = df_bonds.replace('am', '1')  # amide
    df_bonds = df_bonds.replace('ar', '1.5')  # aromatic
    df_bonds = df_bonds.replace('du', '1')  # dummy
    df_bonds = df_bonds.replace('un', '1')  # unknown
    df_bonds = df_bonds.replace('nc', '0')  # not connected
    df_bonds = df_bonds.replace('\n', ' ')
    df_bonds = np.array([np.float(x) for x in df_bonds.split()]).reshape(
        (-1, 4))  # convert the the elements to integer
    df_bonds = pd.DataFrame(
        df_bonds, columns=['bond_id', 'atom1', 'atom2', 'bond_type'])
    df_bonds.set_index(['bond_id'], inplace=True)
    return df_bonds


def compute_edge_attr(edge_index, bonds):
    """
    Compute the edge attributes according to the chemical bonds. 
    """
    sources = edge_index[0, :]
    targets = edge_index[1, :]
    edge_attr = np.zeros((edge_index.shape[1], 1))
    for index, row in bonds.iterrows():
        # find source == row[1], target == row[0]
        # minus one because in new setting atom id starts with 0
        source_locations = set(list(np.where(sources == (row[1]-1))[0]))
        target_locations = set(list(np.where(targets == (row[0]-1))[0]))
        edge_location = list(
            source_locations.intersection(target_locations))[0]
        edge_attr[edge_location] = row[2]

        # find source == row[0], target == row[1]
        source_locations = set(list(np.where(sources == (row[0]-1))[0]))
        target_locations = set(list(np.where(targets == (row[1]-1))[0]))
        edge_location = list(
            source_locations.intersection(target_locations))[0]
        edge_attr[edge_location] = row[2]
    return edge_attr


def form_graph(atoms, bonds, features_to_use, threshold):
    """
    Form a graph data structure (Pytorch geometric) according to the input data frame.
    Rule: Each atom represents a node. If the distance between two atoms are less than or 
    equal to 4.5 Angstrom (may become a tunable hyper-parameter in the future), then an 
    undirected edge is formed between these two atoms. 
    Input:
    atoms: dataframe containing the 3-d coordinates of atoms.
    bonds: dataframe of bond info.
    threshold: distance threshold to form the edge (chemical bond).
    Output:
    A Pytorch-gemometric graph data with following contents:
        - node_attr (Pytorch Tensor): Node feature matrix with shape [num_nodes, num_node_features]. e.g.,
          x = torch.tensor([[-1], [0], [1]], dtype=torch.float)
        - edge_index (Pytorch LongTensor): Graph connectivity in COO format with shape [2, num_edges*2]. e.g.,
          edge_index = torch.tensor([[0, 1, 1, 2],
                                     [1, 0, 2, 1]], dtype=torch.long)
    Forming the final output graph:
        data = Data(x=x, edge_index=edge_index)
    """
    A = atoms.loc[:, 'x':'z']  # sample matrix
    A_dist = distance.cdist(A, A, 'euclidean')  # the distance matrix

    # set the element whose value is larger than threshold to 0
    threshold_condition = A_dist > threshold

    # set the element whose value is larger than threshold to 0
    A_dist[threshold_condition] = 0

    result = np.where(A_dist > 0)
    result = np.vstack((result[0], result[1]))
    edge_attr = compute_edge_attr(result, bonds)
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)
    edge_index = torch.tensor(result, dtype=torch.long)

    # normalize large features
    atoms['x'] = atoms['x']/300
    atoms['y'] = atoms['y']/300
    atoms['z'] = atoms['z']/300

    node_features = torch.tensor(
        atoms[features_to_use].to_numpy(), dtype=torch.float32)
    return node_features, edge_index, edge_attr


def compute_spherical_coord(data):
    """
    Shift the geometric center of the pocket to origin, then compute its spherical coordinates.             
    """
    center = np.mean(data, axis=0)
    shifted_data = data - center  # center the data around origin

    r, theta, phi = cartesian_to_spherical(shifted_data)
    return r, theta, phi


def cartesian_to_spherical(data):
    """Convert cartesian coordinates to spherical coordinates.
    Arguments:   
    data - numpy array with shape (n, 3) which is the 
    cartesian coordinates (x, y, z) of n points.
    Returns:   
    numpy array with shape (n, 3) which is the spherical 
    coordinates (r, theta, phi) of n points.
    """
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    # distances to origin
    r = np.sqrt(x**2 + y**2 + z**2)

    # angle between x-y plane and z
    theta = np.arccos(z/r)/np.pi

    # angle on x-y plane
    phi = np.arctan2(y, x)/np.pi

    #spherical_coord = np.vstack([r, theta, phi])
    #spherical_coord = np.transpose(spherical_coord)
    return r, theta, phi


def extract_sasa_data(siteresidue_list, pop):
    '''extracts accessible surface area data from .pops file generated by POPSlegacy.
        then matches the data in the .pops file to the binding site in the mol2 file.
        Used POPSlegacy https://github.com/Fraternalilab/POPSlegacy '''
    # Extracting sasa data from .pops file
    residue_list = []
    qsasa_list = []
    with open(pop) as popsa:
        for line in popsa:
            line_list = line.split()
            if len(line_list) == 12:
                residue_type = line_list[2] + line_list[4]
                if residue_type in siteresidue_list:
                    qsasa = line_list[7]
                    residue_list.append(residue_type)
                    qsasa_list.append(qsasa)
    qsasa_list = [float(x) for x in qsasa_list]
    median = statistics.median(qsasa_list)
    qsasa_new = [median if x == '-nan' else x for x in qsasa_list]

    # Matching amino acids from .mol2 and .out files and creating dictionary
    qsasa_data = []
    fullprotein_data = list(zip(residue_list, qsasa_new))
    for i in range(len(fullprotein_data)):
        if fullprotein_data[i][0] in siteresidue_list:
            qsasa_data.append(float(fullprotein_data[i][1]))
    return qsasa_data


def extract_seq_entropy_data(siteresidue_list, profile):
    '''extracts sequence entropy data from .profile'''
    # Opening and formatting lists of the probabilities and residues
    with open(profile) as profile:
        ressingle_list = []
        probdata_list = []

        # extracting relevant information
        for line in profile:
            line_list = line.split()
            residue_type = line_list[0]
            prob_data = line_list[1:]
            prob_data = list(map(float, prob_data))
            ressingle_list.append(residue_type)
            probdata_list.append(prob_data)
    ressingle_list = ressingle_list[1:]
    probdata_list = probdata_list[1:]

    # Changing single letter amino acid to triple letter with its corresponding number
    count = 0
    restriple_list = []
    for res in ressingle_list:
        newres = res.replace(res, amino_single_to_triple(res))
        count += 1
        restriple_list.append(newres + str(count))

    # Calculating information entropy
    with np.errstate(divide='ignore'):      # suppress warning
        prob_array = np.asarray(probdata_list)
        log_array = np.log2(prob_array)

        # change all infinite values to 0
        log_array[~np.isfinite(log_array)] = 0
        entropy_array = log_array * prob_array
        entropydata_array = np.sum(a=entropy_array, axis=1) * -1
        entropydata_list = entropydata_array.tolist()

    # Matching amino acids from .mol2 and .profile files and creating dictionary
    fullprotein_data = dict(zip(restriple_list, entropydata_list))
    seq_entropy_data = {k: float(
        fullprotein_data[k]) for k in siteresidue_list if k in fullprotein_data}
    return seq_entropy_data


def amino_single_to_triple(single):
    """
    converts the single letter amino acid abbreviation to the triple letter abbreviation
    """

    single_to_triple_dict = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
                             'G': 'GLY', 'Q': 'GLN', 'E': 'GLU', 'H': 'HIS', 'I': 'ILE',
                             'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                             'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    for i in single_to_triple_dict.keys():
        if i == single:
            triple = single_to_triple_dict[i]
    return triple
