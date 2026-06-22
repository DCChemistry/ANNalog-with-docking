# ANNalog with docking
This project aims to bias the generation of the [ANNalog](https://github.com/DVNecromancer/ANNalog/tree/main) Seq2Seq model towards a given protein target by means of selecting a subset of generated molecules with good docking scores, as docked with Autodock Vina. This subset of molecules will be used for a subsequent round of generation with ANNalog, and the process will repeat for as many times as the user specifies.

# Installation
Using the [pixi package manager](https://github.com/prefix-dev/pixi), simply type the following to install the project environment:
```
pixi install
```