"""
Compute the geometric centers of each pocket by averaging the coordinates of
the atoms
"""
import os
import numpy as np
from biopandas.mol2 import PandasMol2
from tqdm import tqdm
import pickle


def compute_geo_center(mol_path):
    """Compute the geometric center of a pocket"""
    atoms = PandasMol2().read_mol2(mol_path)
    atoms = atoms.df[['x', 'y', 'z']]
    center = np.mean(atoms.to_numpy(), axis=0)
    return list(center)


if __name__ == "__main__":
    pocket_dir = '../../../../data/googlenet-dataset/'
    out_path = './pocket_center.pickle'

    pocket_center = {}
    for x in tqdm(os.listdir(pocket_dir)):
        mol_path = os.path.join(pocket_dir, x)
        if os.path.isdir(mol_path):
            mol_path = mol_path + '/' + x + '.mol2'
            center = compute_geo_center(mol_path)
            pocket_center[x] = center
        else:
            print("something went wrong: {}".format(mol_path))

    # cannot serialize dictionary of lists/numpy array,
    # so use pickle
    with open(out_path, "wb") as f:
        pickle.dump(pocket_center, f, protocol=pickle.HIGHEST_PROTOCOL)
