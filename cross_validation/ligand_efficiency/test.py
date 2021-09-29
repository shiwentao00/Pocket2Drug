from openbabel import openbabel

pdbqt_file = "../../../p2d_results_selfie/cv_results/cross_val_fold_0/val_pockets_sample_clustered_pdbqt/1a3bB00/1a3bB00-1.pdbqt"


obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdbqt", "smi")

mol = openbabel.OBMol()
obConversion.ReadFile(mol, pdbqt_file)   # Open Babel will uncompress automatically
print(mol.GetMolWt())
obConversion.WriteFile(mol, './test.smi')

# C(=O)(N([C][C]C1=[C][C]=[C][C]=C1[C]=[C])[C][C]/[C]=C(/C#C[C][C]O[C])\[C]=O)/[C]=[C]/[C]=[C]/[C]	=
