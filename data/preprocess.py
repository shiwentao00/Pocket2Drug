"""
Pre-process the dataset, and generate a yaml file that contains valid
pocket-SMILES pairs.
"""
import os
import glob
from rdkit import Chem, RDLogger
from rdkit.Chem import MolStandardize
import yaml
from tqdm import tqdm
RDLogger.DisableLog('rdApp.*')


if __name__ == "__main__":
    # list of pocket names
    pocket_dir = "/home/wentao/Desktop/local-workspace/pocket2drug-project/data/googlenet-dataset/"
    pockets = [d for d in os.listdir(
        pocket_dir) if os.path.isdir(os.path.join(pocket_dir, d))]
    print("total number of pocket files: ", len(pockets))

    # protein hashset according to popsa files
    pop_dir = '../../data/pops-googlenet/'
    proteins = [d for d in os.listdir(
        pop_dir) if os.path.isfile(os.path.join(pop_dir, d))]
    proteins = set([x[0:-5] for x in proteins])
    print("total number of corresponding proteins: ", len(proteins))

    # filter out pockets without popsa files
    pockets_with_pops = []
    for pocket in pockets:
        # only append the pockets that are in the protein set
        if pocket[0:-2] in proteins:
            pockets_with_pops.append(pocket)
    pockets = pockets_with_pops
    print("number of pockets with popsa files: ", len(pockets))

    # Count number of ligands. it turns out each pocket file
    # is associated with exactly 1 ligand file.
    num_ligands = 0
    for pocket in pockets:
        pocket_path = pocket_dir + pocket + '/'
        os.chdir(pocket_path)
        sdf_files = glob.glob('*.sdf')
        num_ligands += len(sdf_files)
    print("number of ligands: ", num_ligands)

    # filter out data points that have invalid SMILES
    # which can not be parsed by rdkit
    # hashmap to store the SMILE of ligands of pockets
    good_smiles = {}
    problematic_data = []
    multi_ligand_pockets = []
    normarizer = MolStandardize.normalize.Normalizer()
    choose_frag = MolStandardize.fragment.LargestFragmentChooser()
    for pocket in tqdm(pockets):
        pocket_path = pocket_dir + pocket + '/'
        os.chdir(pocket_path)
        sdf_files = glob.glob('*.sdf')
        assert(len(sdf_files) == 1)
        sdf_file = sdf_files[0]
        sdf_file = pocket_path + sdf_file
        supp = Chem.SDMolSupplier(sdf_file)
        for i, mol in enumerate(supp):
            if i > 0:
                multi_ligand_pockets.append(pocket)
            if mol is None:
                problematic_data.append(pocket)
            else:
                mol = normarizer.normalize(mol)
                mol = choose_frag.choose(mol)
                mol = Chem.MolToSmiles(
                    mol, isomericSmiles=False, canonical=True
                )
                if mol is not None:
                    good_smiles[pocket] = mol
                else:
                    problematic_data.append(pocket)
            # print(mol.GetNumAtoms())

    print("number of problematic pockets: ", len(problematic_data))
    print("number of multi-ligand pockets: ",
          len(set(multi_ligand_pockets)))
    print("number of good pockets: ", len(good_smiles))

    os.chdir(
        "/home/wentao/Desktop/local-workspace/pocket2drug-project/Pocket2Drug/data/"
    )
    out = open('./pocket-smiles.yaml', 'w')
    yaml.dump(good_smiles, out)
