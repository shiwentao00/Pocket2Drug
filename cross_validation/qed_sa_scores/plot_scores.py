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

    p2d_char = {
        "title": "Pocket2Drug - character tokenization",
        "qed_file": "./scores/p2d_char_qed_scores.yaml",
        "sa_file": "./scores/p2d_char_sa_scores.yaml"
    }

    p2d_selfies = {
        "title": "Pocket2Drug - selfies tokenization",
        "qed_file": "./scores/p2d_selfies_qed_scores.yaml",
        "sa_file": "./scores/p2d_selfies_sa_scores.yaml"
    }
    
    configs = [chembel28, p2d_char, p2d_selfies]

    fig, (ax0, ax1, ax2) = plt.subplots(
        nrows=1, ncols=3, figsize=(18, 6), sharey=True)

    for i, ax in enumerate([ax0, ax1, ax2]):
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
    plt.savefig("./scores.png", bbox_inches="tight", dpi=600)
