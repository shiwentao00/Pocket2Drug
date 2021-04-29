"""
Parse the DUD-E dataset txt file to get the pdb ids.
"""
import yaml

if __name__ == "__main__":
    proteins = []
    with open('./dude.txt', 'r') as f:
        lines = f.readlines()
    for line in lines:
        proteins.append(line.split()[2])
    
    with open('./dude_proteins.yaml', 'w') as f:
        yaml.dump(proteins, f)