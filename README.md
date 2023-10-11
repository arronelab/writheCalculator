# SWRITHE
iPython notebook and supporting data to accompany "The SKMT Algorithm: A method for assessing and comparing
underlying protein entanglement"

[![DOI](https://zenodo.org/badge/617996455.svg)](https://zenodo.org/badge/latestdoi/617996455)

This notebook allows you to compute writhe and average crossing number fingerprints of smoothed representations of protein backbones. The code downloads PDB files, and smoothes the alpha-carbon backbone to a minimal representation reserving the underlying topological structure. With these writhe fingerprints, once can analyse the backbone for large scale helical structures, and compare it to a database of smoothed structures to identify similarities.



## Installation 
We recommend using the SWRITHE package locally with Linux or MacOS. To install, run the below in the terminal.
```shell
git clone https://github.com/arronelab/SWRITHE.git
cd SWRITHE
make
```

This notebook is maintained by Arron Bale, email me at arron.n.bale@durham.ac.uk if you need any help with installation or potential new features
