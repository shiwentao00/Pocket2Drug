## Cross Validation
Docking is the main method used to evaluate the results. The dataset is divided into 10 folds and evaluated by a 10-fold cross validation. In this way, the model is evaluated against all the pockets in the dataset. The cross validation has the following procedures:
	
    1. Exclude the data points in DUDE dataset. 
	2. Exclude the 5 pockets for case study.
	3. Divide the dataset into 10 folds.
	4. Train the model 10 times for the 10 folds.
	5. Sample 20480 molecules for each pocket in the validation foldes.
    6. Compute a representative subset from the sampled molecules of each pocket. Right now the Maxmin picking algorithm is used. It is implemented in RdKit. 
	<!-- 6. Compute the binary string of sampled molecules for each pocket in the validation folds.
	7. For each pocket, cluster the molecules using the subset algorithm, select representative molecules and save them. -->
	7. Compute molecular weight for each of the selected molecules.
	8. Prepare pdbqt files for all the pockets in the dataset.
	9. Prepare a large dataset of random drugs from Zinc/Chembl.
	10. For fold 0, develop a set of procedures to perform docking and collect results:
	For each pocket in validation fold 0,
		i. Dock the sampled molecules
		ii. Get the same number of random molecules from Zinc, dock them.
		iii. Collect the docking scores.
		iv. Normalize the docking scores with molecular weight.
		v. Save both normalized and unnormalized docking scores.
	Draw the distributions of random and our model.
Repeat step 7 for the rest validation folds.

### Docking
The docking script has the following steps:
	1. Convert SMILES to pdbqt file.
	2. Compute the docking center and the docking box sizes, store the docking center and box dimension in a temporary file.
    3. Use Vina to dock.