"""
This script prepare the data for a 10-fold cross validation. 
 1. The proteins from the DUD-E dataset is excluded.
 2. The pockets used in previous case studies are excluded.
 3. The pockets are shuffled and divided into 10 folds.
 4. Each fold is saved
"""


if __name__ == "__main__":