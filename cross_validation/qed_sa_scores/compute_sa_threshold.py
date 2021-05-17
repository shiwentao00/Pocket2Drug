"""
Compute the threshold of SA score that can cover a certain ratio 
of molecules in chembl.
"""
from compute_scores import read_smiles_file

if __name__ == "__main__":
    chembl28_dataset_path = "/home/wentao/Desktop/local-workspace/molecule-generator-project/Molecule-RNN/dataset/chembl28-cleaned.smi"

    chembl28_dataset = read_smiles_file(chembl28_dataset_path)

    