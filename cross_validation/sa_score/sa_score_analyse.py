"""
Analyse the Synthetic Accessibility scores of:
    1. Chembl28 dataset
    2. Molecules sampled from fold 0 pockets by Pocket2Drug, (char-level tokenization)
    3. Molecules sampled from fold 0 pockets by Pocket2Drug, (selfies tokenization)
"""
from sascorer import calculateScore


if __name__ == "__main__":
    chemb28_dataset = "~/Desktop/local-workspace/molecule-generator-project/Molecule-RNN/dataset/chembl28-cleaned.smi"
    p2d_char_dataset = "/home/wentao/Desktop/local-workspace/pocket2drug-project/p2d_results/cv_results/cross_val_fold_0/val_pockets_sample_clustered/"
    p2d_selfies_dataset = ""
