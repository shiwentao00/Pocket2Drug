"""Compute the mean of each test pocket's sequence identity with the train pockets."""
from tqdm import tqdm
from statistics import mean, median

if __name__ == "__main__":
    selection_function = "median"
    if selection_function=="mean":
        select_func = mean
    elif selection_function=="median":
        select_func = median
    elif selection_function=="max":
        select_func = max

    # read the file 
    print("loading the file into memory...")
    with open("../../../sequence_identity/sequence_identity_fold0.txt", "r") as f:
        lines = f.readlines()

    # update the result line by line 
    seq_identity = {}
    curr_pocket = None
    curr_scores = []
    for line in tqdm(lines):
        line = line.rstrip()
        line = line.split()
        test_pocket = line[0]
        score = line[2]
        if test_pocket != curr_pocket:
            if curr_pocket and curr_scores:
                seq_identity[curr_pocket] = select_func(curr_scores)
            curr_scores = []
            curr_pocket = test_pocket
        curr_scores.append(float(score))

    pockets = []
    for pocket in seq_identity:
        pockets.append((seq_identity[pocket], pocket))

    pockets.sort()
    with open(f"./{selection_function}_seq_identities.txt", "w") as f:
        for score, pocket in pockets:
            f.write(f"{score} {pocket}\n")