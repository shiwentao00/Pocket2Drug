"""
Plot normalized histograms of the QED and SA scores.
"""
import yaml
import matplotlib.pyplot as plt


if __name__ == "__main__":
    chembel28 = {
        "title": "Chembl28",
        "qed_file": "./scores/chembl28_qed_scores.yaml",
        "sa_file": "./scores/chembl28_sa_scores.yaml"
    }

    p2d_dataset = {
        "title": "P2D - dataset",
        "qed_file": "./scores/p2d_dataset_qed_scores.yaml",
        "sa_file": "./scores/p2d_dataset_sa_scores.yaml"
    }

    p2d_char = {
        "title": "P2D - Character",
        "qed_file": "./scores/p2d_char_qed_scores.yaml",
        "sa_file": "./scores/p2d_char_sa_scores.yaml"
    }

    p2d_selfies = {
        "title": "P2D - SELFIES",
        "qed_file": "./scores/p2d_selfies_qed_scores.yaml",
        "sa_file": "./scores/p2d_selfies_sa_scores.yaml"
    }
    
    configs = [chembel28, p2d_dataset, p2d_char, p2d_selfies]

    # ax0 for qed, ax1 for sa
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    ax0.set_xlabel("QED")
    ax0.set_ylabel("Frequency")
    ax1.set_xlabel("SA")
    ax1.set_ylabel("Frequency")

    for i in range(4):
        config = configs[i]
        
        with open(config["qed_file"], "r") as f:
            qed_scores = yaml.full_load(f)

        if i == 0:
            ax0.hist(
                qed_scores,
                bins=20,
                label=config['title'],
                histtype='bar',
                facecolor="grey",
                alpha=0.5,
                stacked=True,
                density=True,
                fill=True
            )
        else:
            ax0.hist(
                qed_scores,
                bins=20,
                label=config['title'],
                histtype='step',
                stacked=True,
                density=True,
                fill=False
            )
        
        with open(config["sa_file"], "r") as f:
            sa_scores = yaml.full_load(f)

        if i == 0:
            ax1.hist(
                sa_scores,
                bins=20,
                label=config['title'],
                histtype='bar',
                facecolor="grey",
                alpha=0.5,
                stacked=True,
                density=True,
                fill=True
            )        
        else:    
            ax1.hist(
                sa_scores,
                bins=20,
                label=config['title'],
                histtype='step',
                stacked=True,
                density=True,
                fill=False
            )
    
    ax0.legend(loc='upper left')
    ax1.legend(loc='upper right')
    # fig.suptitle('QED vs SA')
    plt.savefig("./figures/scores_hist_norm.png", bbox_inches="tight", dpi=600)
    
