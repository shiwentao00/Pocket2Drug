"""Rank the pockets according to their mean sequence identity with training data"""
import pickle

if __name__ =="__main__":
    with open("../../../sequence_identity/mean_seq_id.pickle", "rb") as f:
        mean_seq_identity = pickle.load(f)

    pockets = []
    for pocket in mean_seq_identity:
        pockets.append((mean_seq_identity[pocket], pocket))

    pockets.sort()
    with open("./mean_seq_identities.txt", "w") as f:
        for score, pocket in pockets:
            f.write(f"{score} {pocket}\n")