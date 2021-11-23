"""Compute the max of each test pocket's sequence identity with the train pockets."""
from tqdm import tqdm
import pickle

if __name__ == "__main__":
    # read the file 
    print("loading the file into memory...")
    with open("../../../sequence_identity/sequence_identity_fold0.txt", "r") as f:
        lines = f.readlines()

    # update the result line by line 
    max_seq_identity = {}
    curr_pocket = None
    curr_scores = []
    for line in tqdm(lines):
        line = line.rstrip()
        line = line.split()
        test_pocket = line[0]
        score = line[2]
        if test_pocket != curr_pocket:
            if curr_pocket and curr_scores:
                max_seq_identity[curr_pocket] = max(curr_scores)
            curr_scores = []
            curr_pocket = test_pocket
        curr_scores.append(float(score))

    with open("../../../sequence_identity/max_seq_id.pickle", "wb") as f:
        pickle.dump(max_seq_identity, f)