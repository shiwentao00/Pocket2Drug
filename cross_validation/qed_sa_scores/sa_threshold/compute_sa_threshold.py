"""
Compute the threshold of SA score that can cover a certain ratio 
of molecules in chembl.
"""
import yaml

def load_yaml(path):
    with open(path, 'r') as f:
        ans = yaml.full_load(f)
    return ans


def find_percentage(scores, th):
    """percentage of data below threshold th, scores is already sorted"""
    num_data = len(scores)
    for i in range(num_data):
        if scores[i] >= th:
            return float(i) / num_data 
    return 1.0


if __name__ == "__main__":
    p2d_dataset_sa_score_path = "../scores/p2d_dataset_sa_scores.yaml"
    p2d_dataset_sa_scores = load_yaml(p2d_dataset_sa_score_path)
    p2d_dataset_sa_scores.sort()

    p2d_char_sa_score_path = "../scores/p2d_char_sa_scores.yaml"
    p2d_char_sa_scores = load_yaml(p2d_char_sa_score_path)
    p2d_char_sa_scores.sort()
    
    p2d_selfies_sa_score_path = "../scores/p2d_selfies_sa_scores.yaml"
    p2d_selfies_sa_scores = load_yaml(p2d_selfies_sa_score_path)
    p2d_selfies_sa_scores.sort()

    sa_score_path = "../scores/chembl28_sa_scores.yaml"
    with open(sa_score_path, 'r') as f:
        sa_scores = yaml.full_load(f)

    ratio_list = [0.9999, 0.999, 0.99, 0.98, 0.95, 0.90, 0.85, 0.8]
    
    sa_scores.sort()
    num_mols = len(sa_scores)
    
    for r in ratio_list:
        th = sa_scores[int(num_mols * r)]
        print("percentage of chembl28 data: {}, threshold: {}".format(r, th))
        print("percentage of p2d dataset: {}".format(
            find_percentage(p2d_dataset_sa_scores, th)))
        print("percentage of p2d char model: {}".format(
            find_percentage(p2d_char_sa_scores, th)))
        print("percentage of p2d selfies model: {}".format(
            find_percentage(p2d_selfies_sa_scores, th)))
        print('-----------------------------------------------------')

