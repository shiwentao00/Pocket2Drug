
## Usage
### Dependencies
Pocket2Drug is implemented based on the folloing Python packages:
1. [Pytorch](https://pytorch.org/get-started/locally/)
2. [Pytorch-geometric](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html)
3. [Rdkit](https://www.rdkit.org/docs/Install.html)
4. [SELFIES](https://github.com/aspuru-guzik-group/selfies)
5. Pandas 
6. [BioPandas](http://rasbt.github.io/biopandas/)
7. Numpy
8. Scipy

### Installation
1. Install Pytorch
```
conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
```

2. Install Pytorch geometric:

3. Install Biopandas:
```
conda install biopandas -c conda-forge
```

4. Install SELFIES:
```
pip install selfies
```

### Dataset
All the related data can be downloaded [here](). There are two dataset files:
1. dataset.tar.gz: contains all binding site data in this project.
2. pops.tar.gz: contains information of node feature contact surface area.

### Configuration

### Train
The configurations for training can be updated in ```train.yaml```. Modify the ```pocket_dir``` and ```pop_dir``` entries to the paths of the extracted dataset. Modify the ```out_dir``` entry to the folder where you want to save the output results. Then,
```
python train.py
```

### Inference
After training, the trained model will be saved at ```out_dir```, and we can use it to sample predicted molecules:
```
python sample.py -batch_size 1024 -num_batches 1 -pocket_dir path_to_dataset_folder -popsa_dir path_to_pops_folder
```

