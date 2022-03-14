# Pocket2Drug
Pocket2Drug is an encoder-decoder deep neural network that predicts binding drugs given protein binding sites (pockets). The pocket graphs are generated using [Graphsite](https://github.com/shiwentao00/Graphsite). The encoder is a graph neural network, and the decoder is a recurrent neural network. The [SELFIES](https://github.com/aspuru-guzik-group/selfies) molecule representation is used as the tokenization scheme instead of SMILES. The pipeline of Pocket2Drug is illustrated below:
<p align="center">
<img width="820" height="260" src="doc/pipeline.png">
</p>

## Usage
Please refer the instructions [here]() to run Pocket2Drug.

## Citation
Shi, Wentao, et al. "Pocket2Drug: An encoder-decoder deep neural network for the target-based drug design." Frontiers in Pharmacology: 587.