"""
Plot QED score vs SA score
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
        "title": "P2D - Char",
        "qed_file": "./scores/p2d_char_qed_scores.yaml",
        "sa_file": "./scores/p2d_char_sa_scores.yaml"
    }

    p2d_selfies = {
        "title": "P2D - SELFIES",
        "qed_file": "./scores/p2d_selfies_qed_scores.yaml",
        "sa_file": "./scores/p2d_selfies_sa_scores.yaml"
    }
    
    configs = [chembel28, p2d_dataset, p2d_char, p2d_selfies]

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(
        nrows=1, ncols=4, figsize=(24, 6), sharey=True)

    for i, ax in enumerate([ax0, ax1, ax2, ax3]):
        config = configs[i]

        with open(config["qed_file"], "r") as f:
            qed_scores = yaml.full_load(f)

        with open(config["sa_file"], "r") as f:
            sa_scores = yaml.full_load(f)

        ax.scatter(qed_scores, sa_scores, c="g", s=0.01, marker=".")
        ax.set_title(config["title"])
        ax.set_xlabel("QED")
        ax.set_ylabel("SA")

    # fig.suptitle('QED vs SA')
    plt.savefig("./figures/scores.png", bbox_inches="tight", dpi=600)
