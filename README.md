# Pocket2Drug
Pocket2Drug is an encoder-decoder deep neural network that predicts binding drugs given protein binding sites (pockets). The pocket graphs are generated using [Graphsite](https://github.com/shiwentao00/Graphsite). The encoder is a graph neural network, and the decoder is a recurrent neural network. The [SELFIES](https://github.com/aspuru-guzik-group/selfies) molecule representation is used as the tokenization scheme instead of SMILES. The pipeline of Pocket2Drug is illustrated below:
<p align="center">
<img width="820" height="260" src="doc/pipeline.png">
</p>

If you find Pocket2Drug helpful, please cite our paper in your work :)  
Shi, Wentao, et al. "Pocket2Drug: An encoder-decoder deep neural network for the target-based drug design." Frontiers in Pharmacology: 587.

## Usage
### Dependency installation
1. Install [Pytorch](https://pytorch.org/get-started/locally/):
```
conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
```

2. Install [Pytorch-geometric](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html):
```
conda install pyg -c pyg
```

3. Install [BioPandas](http://rasbt.github.io/biopandas/):
```
conda install biopandas -c conda-forge
```

4. Install [selfies](https://github.com/aspuru-guzik-group/selfies):
```
pip install selfies
```

5. Install [Rdkit](https://www.rdkit.org/docs/Install.html):
```
conda install rdkit -c conda-forge
```

### Dataset
All the related data can be downloaded [here](https://osf.io/qacwj/). After extraction, there will be two folders:
1. pocket-data: files that contain information of the pockets. We will use the ```.mol2``` files.
2. protein-data: files that contain information of the proteins. We wiil use the ```.pops``` and ```.profile``` files.

### Train
The configurations for training can be updated in ```train.yaml```. Set the ```pocket_dir``` to the path of ```pocket-data```, then set ```pop_dir``` and ```profile_dir``` to the path of ```protein-data```. Set the ```out_dir``` the folder where you want to save the output results. The other configurations are for hyper-parameter tuning and they are self-explanatory according to their names. The script ```train.py``` trains the model on a 90%-10% split of the dataset, and you can specify which fold is used for validation:  
```
python train.py -val_fold 0
```
In addition, you can use a pretrained RNN to initialize the decoder, the pretrained model can be found [here](https://osf.io/qacwj/). The pretrained RNN is trained on the chembl dataset and can improve the performance of the model. I have wrote an exmaple for pretraining RNN [here](https://github.com/shiwentao00/Molecule-RNN)).

<!---
Note: the results presented in our research paper were produced with Pytorch v1.7.1 and selfies v1.0, which are far behind current version. There has been a major update for selfies, so the generated token vocabulary of the molecules has changed and the pre-trained RNN can not be used for training (The pre-trained RNN learns the distribution of the chembl database, which improves the performance of the model. I have wrote an exmaple for pretraining RNN [here](https://github.com/shiwentao00/Molecule-RNN)). As a result, the performance of the model will be affected. I will update the configuration once I get a chance to re-do the pretraining and hyper-parameter tuning.
-->

### Sample molecules
After training, the trained model will be saved at ```out_dir```, and we can use it to sample molecules for the pockets in the validation fold:  
```
python sample.py -batch_size 1024 -num_batches 2 -pocket_dir path_to_dataset_folder -popsa_dir path_to_pops_folder -profile_dir path_to_profile_folder -result_dir path_to_training_output_folder -fold 0
```
Of course, the model can be used to sample molecules for the unseen pockets defined by user. Simply omit the ```-fold``` option, the code will run on the specified input directories:
```
python sample.py -batch_size 1024 -num_batches 2 -pocket_dir path_to_dataset_folder -popsa_dir path_to_pops_folder -profile_dir path_to_profile_folder -result_dir path_to_training_output_folder
```
