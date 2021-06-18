## Cross Validation
Docking is the main method used to evaluate the results. The dataset is divided into 10 folds and evaluated by a 10-fold cross validation. In this way, the model is evaluated against all the pockets in the dataset. The cross validation has the following procedures:
	
1. Exclude the data points in DUDE dataset. 
2. Exclude the 5 pockets for case study.
3. Filter the dataset with Synthetic Accessibility (SA) score. We keep the pockets that have binding molecules of SA scores between 1 to 6.
4. Divide the dataset into 10 folds.
5. Train the model 10 times for the 10 folds.
6. Sample 20480 molecules for each pocket in the validation foldes. Filter the molecules with SA score range 1 to 6.
7. Compute a representative subset from the sampled molecules of each pocket. Right now the Maxmin pickingalgorithm is used. It is implemented in RdKit. 
```
cd ./subset
python compute_subset.py -mol_dir <path/to/sampled/molecules> -out_dir <path/of/subset/molecules> 
```
8. Compute molecular weight for each of the selected molecules.
9. Prepare pdbqt files for all the pockets in the dataset.
10. Prepare a large dataset of random drugs from Zinc/Chembl.
11. For fold 0, develop a set of procedures to perform docking and collect results:
For each pocket in validation fold 0,   
	i. Dock the sampled molecules    
	ii. Get the same number of random molecules from Zinc, dock them.   
	iii. Collect the docking scores.    
	iv. Normalize the docking scores with molecular weight.   
	v. Save both normalized and unnormalized docking scores.      
	iv. Draw the distributions of random and our model.   

The folder structure of each fold looks like this:
```
    .
    ├── p2d_results      
        ├── cv_results
            ├── cross_val_fold_0
	        ├── val_pockets_sample
                ├── val_pockets_sample_clustered
                ├── val_pockets_sample_clustered_smi
                    ├── 1a3eB01-.smi
                    ├── 2hr7A02-.smi
                    ├── ...
                    └── 8icsC00-.smi
                ├── val_pockets_sample_clustered_pdbqt
                    ├── 1a3eB01
                        ├── 1a3eB01-1.pdbqt
                        ├── 1a3eB01-2.pdbqt
                        ├── ...
                        └── 1a3eB01-100.pdbqt
                    ├── 2hr7A02
                    ├── ...
                    └── 8icsC00
            ├── cross_val_fold_1
            ├── ...
            └── cross_val_fold_9
	└── random_drugs
	    ├── smiles
	    └── pdbqt
```

11.Repeat step 10 for the rest validation folds.

### Docking
The docking script has the following steps:
	1. Convert SMILES to pdbqt file.
	2. Compute the docking center and the docking box sizes, store the docking center and box dimension in a temporary file.
    3. Use Vina to dock.
