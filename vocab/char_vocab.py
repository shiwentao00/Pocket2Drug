"""
Read the SMILES of good pockets and verify the char-level vocab.
"""
import yaml
import tqdm


def read_smiles(pocket_list_dir, smile_dir):
    """read all the used SMILES as a list of strings"""


if __name__ == "__main__":
    pocket_list_dir = "../../smiles/pocket_that_has_popsa_list.yaml"
    smiles_dir = "../../smiles/smile_files/"

    smiles = read_smiles(pocket_list_dir, smiles_dir)
