"""
Plot the best docking scores of the molecules generated by the model
and the random molecules sampled from Zinc.
"""
import os
import argparse
from os import listdir
from os.path import isfile, join
from tqdm import tqdm
import yaml
import matplotlib.pyplot as plt


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-model_mols_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_docking_results_ligand_efficiency/",
                        help="directory of docking scores of molecules generated by model")

    parser.add_argument("-zinc_mols_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/zinc_docking_results_ligand_efficiency/",
                        help="directory of docking scores of molecules sampled from Zinc")

    parser.add_argument("-output_images_dir",
                        required=False,
                        default="../../../p2d_results_selfie/cv_results/cross_val_fold_0/ligand_efficiency_plots/",
                        help="directory of docking scores of molecules sampled from Zinc")        

    return parser.parse_args()

if __name__=="__main__":
    args = get_args()
    model_mols_dir = args.model_mols_dir
    zinc_mols_dir = args.zinc_mols_dir
    output_images_dir = args.output_images_dir

    if not os.path.exists(output_images_dir):
        os.makedirs(output_images_dir)

    model_docking_score_files = [f for f in listdir(
                model_mols_dir) if isfile(join(model_mols_dir, f))]

    # tuples of docking scores
    best_docking_scores = []
    num_problematic_pocket = 0
    print("loading docking scores...")
    for docking_file in tqdm(model_docking_score_files):
        # docking scores of model molecules 
        with open(join(model_mols_dir, docking_file), "r") as f:
            model_docking_scores = yaml.full_load(f)
        
        # docking scores of random zinc molecules
        with open(join(zinc_mols_dir, docking_file), "r") as f:
            zinc_docking_scores = yaml.full_load(f)
        
        if zinc_docking_scores and model_docking_scores:
            fig, ax = plt.subplots()
            ax.set_title('Pocket: ' + docking_file.split('_')[0])
            #ax.set_xlabel('Docking Score')
            ax.set_ylabel('Ligand Efficiency')
            ax.set_xticklabels(['model', 'zinc'])
            ax.boxplot([list(model_docking_scores.values()), list(zinc_docking_scores.values())], widths=[0.5, 0.5])
            ax.legend()
            plt.savefig(join(output_images_dir, docking_file.split('_')[0] + ".png"), dpi=300)
            plt.close()
        else:
            num_problematic_pocket += 1

    print("finished, number of problematic pocket: ", num_problematic_pocket)        